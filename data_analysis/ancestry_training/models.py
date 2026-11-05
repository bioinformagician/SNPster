"""Reference-only base fits, training-only heads, and held-out inference."""
import gzip
import pickle
import re
import shutil
import sys
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.mixture import GaussianMixture

from core import digest, java_heap, probability, read, rounded_flare_probability, run, verify_tools, write, write_vcf


def dump(path, value):
    with Path(path).open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load(path):
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def base_policy(model):
    params = dict(model["params"])
    params.pop("crf_c2", None)
    if model["family"] == "flare":
        params = {key: value for key, value in params.items()
                  if key in ("reference_cap", "reference_source_cap", "model_file", "components", "component_map", "reference_label")}
    elif model["family"] in ("pca_ridge", "pca_rf", "pca_xgboost"):
        params = {key: value for key, value in params.items() if key in ("reference_cap", "n_components")}
    return {"family": model["family"], "params": params}


def plan(config):
    import hashlib
    import json
    result = {"models": [], "bases": {}}
    for model in config["models"]:
        base = base_policy(model)
        identity = hashlib.sha256(json.dumps(base, sort_keys=True).encode()).hexdigest()[:16]
        result["bases"][identity] = base
        result["models"].append(dict(model, base_id=identity))
    return result


def reference_subset(meta, params):
    reference = meta["reference"]
    cap = params.get("reference_cap")
    source_cap = params.get("reference_source_cap")
    if source_cap is not None:
        if cap is not None or type(source_cap) is not int or source_cap < 1:
            raise ValueError("Source-balanced capping requires a positive integer and no fine-population cap")
        selected = []
        for source in meta["sources"]:
            populations = {}
            for index, row in enumerate(reference):
                if row["source"] == source:
                    if not row["fine"]:
                        raise ValueError("Source-balanced reference requires fine-population labels")
                    populations.setdefault(row["fine"], []).append(index)
            for population in populations.values():
                population.sort(key=lambda index: reference[index]["sample_id"])
            count = 0
            while any(populations.values()) and count < source_cap:
                for fine in sorted(populations):
                    if populations[fine] and count < source_cap:
                        selected.append(populations[fine].pop(0))
                        count += 1
            if count == 0:
                raise ValueError("Source-balanced reference loses a declared source")
        return sorted(selected)
    if cap is None:
        return list(range(len(reference)))
    if not isinstance(cap, int) or cap < 1:
        raise ValueError("reference_cap must be a positive per-fine-population count")
    groups = {}
    for index, row in enumerate(reference):
        if not row["fine"]:
            raise ValueError("Reference capping requires authoritative fine-population labels")
        groups.setdefault(row["fine"], []).append(index)
    chosen = [index for fine in sorted(groups)
              for index in sorted(groups[fine], key=lambda index: reference[index]["sample_id"])[:cap]]
    if {reference[index]["source"] for index in chosen} != set(meta["sources"]):
        raise ValueError("Reference cap loses source coverage")
    return chosen


def flare_reference_labels(meta, params, selected):
    rows = [meta["reference"][index] for index in selected]
    if "reference_label" not in params:
        return [row["source"] for row in rows]
    if params["reference_label"] != "fine" or "model_file" in params:
        raise ValueError("Fine reference labels cannot be combined with a supplied native model")
    populations = {}
    for row in rows:
        fine, source = row["fine"], row["source"]
        if (not re.fullmatch(r"[A-Za-z0-9_]+", fine) or source not in meta["sources"]
                or (fine in populations and populations[fine] != source)):
            raise ValueError("Fine reference labels are missing, unsafe or map to multiple sources")
        populations[fine] = source
    components = params.get("components", [])
    if len(components) != len(populations) or set(components) != set(populations):
        raise ValueError("Fine components must exactly cover the selected reference populations")
    expected = [[float(populations[name] == source) for source in meta["sources"]] for name in components]
    mapping = probability(params.get("component_map", []), len(meta["sources"]))
    if not np.array_equal(mapping, expected):
        raise ValueError("Fine component map must be the exact reference-label one-hot mapping")
    if set(populations.values()) != set(meta["sources"]):
        raise ValueError("Fine reference labels lose source coverage")
    return [row["fine"] for row in rows]


def layout(meta, params):
    positions = np.asarray([key[1] for key in meta["keys"]])
    points = np.loadtxt(meta["config"]["maps"][str(meta["chromosome"])], dtype=str, ndmin=2)
    points = points[[str(row[0]).removeprefix("chr") == str(meta["chromosome"]) for row in points]]
    bp, cm = points[:, 3].astype(float), points[:, 2].astype(float)
    if (len(bp) < 2 or (np.diff(bp) <= 0).any() or (np.diff(cm) < 0).any()
            or bp[0] > positions[0] or bp[-1] < positions[-1] or cm[-1] <= 0):
        raise ValueError("Genetic map must monotonically bracket the exact SNP panel")
    size = int(round(params.get("window_size_cM", .5) * len(positions) / cm[-1]))
    count = len(positions) // max(size, 1)
    if size < 2 or count < 2:
        raise ValueError("Insufficient marker density for the requested Gnomix window")
    if int(size * params.get("context_ratio", .5)) < 1:
        raise ValueError("Gnomix context must contain at least one SNP")
    starts = np.arange(count) * size
    ends = np.r_[starts[1:], len(positions)]
    from core import physical_weights
    cells = physical_weights(positions, meta["span"])
    return {"M": size, "W": count, "starts": starts, "ends": ends,
            "weights": np.asarray([cells[a:b].sum() for a, b in zip(starts, ends)]),
            "positions": positions}


def gnomix_classes(config):
    root = Path(config["tools"]["gnomix_root"])
    expected = "cd15f65eadfca8e2bbf209eb69a82978a200c5fd"
    actual = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    dirty = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain"], text=True)
    if actual != expected or dirty:
        raise ValueError(f"Gnomix requires the clean verified upstream {expected}")
    if "academic research use only" not in (root / "LICENSE.md").read_text():
        raise ValueError("Gnomix licensing notice is missing")
    sys.path.insert(0, str(root))
    from src.model import Gnomix
    from src.preprocess import WindowedNpyData
    return Gnomix, WindowedNpyData


def pca_features(reference, params, seed):
    dosage = reference.sum(axis=1).astype(float)
    frequency = dosage.mean(axis=0) / 2
    scale = np.sqrt(np.maximum(2 * frequency * (1 - frequency), .01))
    count = params.get("n_components", 16)
    if not isinstance(count, int) or not 1 <= count < min(dosage.shape):
        raise ValueError("PCA n_components must be below both reference sample and marker counts")
    model = PCA(n_components=count, svd_solver="randomized", random_state=seed)
    model.fit((dosage - 2 * frequency) / scale)
    return model, frequency, scale


def transform(base, calls):
    return base["pca"].transform((calls.sum(axis=1) - 2 * base["frequency"]) / base["scale"])


def fit_base(panel, base_id, output):
    output.mkdir()
    meta = read(panel / "panel.json")
    verify_tools(meta)
    spec = plan(meta["config"])["bases"][base_id]
    params, family = spec["params"], spec["family"]
    data = np.load(panel / "data.npz")
    selected = reference_subset(meta, params)
    reference = data["reference"][selected]
    labels = np.asarray([meta["sources"].index(meta["reference"][index]["source"]) for index in selected])
    base = {"family": family, "params": params, "sources": meta["sources"], "base_id": base_id,
            "panel_sha256": digest(panel / "data.npz"), "reference_indices": selected}
    started = time.monotonic()
    if family == "flare":
        reference_labels = flare_reference_labels(meta, params, selected)
        write_vcf(output / "reference.vcf.gz", [tuple(key) for key in meta["keys"]],
                  [meta["reference"][index]["sample_id"] for index in selected], reference)
        (output / "reference.panel").write_text("".join(
            f"{meta['reference'][index]['sample_id']}\t{label}\n"
            for index, label in zip(selected, reference_labels)))
        base["fitting_semantics"] = "Frozen labelled reference; target-dependent FLARE parameters are inferred later"
    elif family == "gnomix":
        cls, windowed = gnomix_classes(meta["config"])
        grid = layout(meta, params)
        haplotypes = reference.reshape(-1, reference.shape[-1])
        x, y = output / "X.npy", output / "y.npy"
        np.save(x, haplotypes)
        np.save(y, np.repeat(np.repeat(labels, 2)[:, None], len(meta["keys"]), axis=1).astype(np.int8))
        model = cls(C=len(meta["keys"]), M=grid["M"], A=len(meta["sources"]), S=1,
                    mode="fast", snp_pos=grid["positions"],
                    snp_ref=np.asarray([key[2] for key in meta["keys"]]),
                    snp_alt=np.asarray([key[3] for key in meta["keys"]]),
                    population_order=np.asarray(meta["sources"]), missing_encoding=2,
                    n_jobs=1, calibrate=False, context_ratio=params.get("context_ratio", .5),
                    mode_filter=False, seed=meta["config"]["seed"], verbose=False)
        model.base.base_multithread = False
        for estimator in model.base.models:
            estimator.set_params(C=params.get("logistic_c", 1), random_state=meta["config"]["seed"])
        model.base.train_windowed(windowed([str(x)], [str(y)], len(meta["keys"]), grid["M"]), verbose=False)
        base.update(model=model, layout=grid)
        x.unlink()
        y.unlink()
    elif family == "nnls":
        base["frequencies"] = np.stack(
            [reference[labels == index].mean(axis=(0, 1)) for index in range(len(meta["sources"]))], axis=1)
    else:
        pca, frequency, scale = pca_features(reference, params, meta["config"]["seed"])
        base.update(pca=pca, frequency=frequency, scale=scale)
        if family == "pca_gmm":
            count = params.get("k", len(meta["sources"]))
            if not isinstance(count, int) or not 2 <= count <= len(reference):
                raise ValueError("GMM k must be an integer between 2 and the reference sample count")
            gmm = GaussianMixture(n_components=count, covariance_type="full", n_init=5,
                                  reg_covar=params.get("reg_covar", 1e-5),
                                  random_state=meta["config"]["seed"])
            features = transform(base, reference)
            responsibilities = gmm.fit_predict(features)
            if not gmm.converged_:
                raise ValueError("Reference-only GMM did not converge")
            soft = gmm.predict_proba(features)
            mapping = soft.T @ np.eye(len(meta["sources"]))[labels]
            if (mapping.sum(axis=1) == 0).any() or len(set(responsibilities)) < count:
                raise ValueError("GMM produced an unsupported/empty cluster")
            base.update(gmm=gmm, component_map=probability(mapping / mapping.sum(axis=1, keepdims=True)),
                        component_names=[f"cluster_{index}" for index in range(count)])
    dump(output / "base.pkl", base)
    write(output / "base.json", {
        "base_id": base_id, "family": family, "fitting_role": "reference_only",
        "sources": meta["sources"], "reference_people": len(selected),
        "panel_sha256": base["panel_sha256"], "seconds": time.monotonic() - started,
        "sha256": digest(output / "base.pkl"),
        "artifact_sha256": {path.name: digest(path) for path in output.iterdir()
                           if path.name in ("base.pkl", "reference.vcf.gz", "reference.vcf.gz.csi", "reference.panel")},
    })


def regression(kind, params, seed):
    if kind in ("pca_ridge", "ridge"):
        return Ridge(alpha=params.get("alpha", 1))
    if kind in ("pca_rf", "rf"):
        return RandomForestRegressor(n_estimators=params.get("n_estimators", 200),
                                     min_samples_leaf=params.get("min_samples_leaf", 5),
                                     max_depth=params.get("max_depth"), random_state=seed, n_jobs=1)
    if kind in ("pca_xgboost", "xgboost"):
        from sklearn.multioutput import MultiOutputRegressor
        from xgboost import XGBRegressor
        return MultiOutputRegressor(XGBRegressor(
            n_estimators=params.get("n_estimators", 200), max_depth=params.get("max_depth", 3),
            learning_rate=params.get("learning_rate", .05), random_state=seed, n_jobs=1), n_jobs=1)
    raise ValueError(f"Unknown regression head: {kind}")


def fit_model(panel, base_folder, model_id, output):
    output.mkdir()
    meta = read(panel / "panel.json")
    spec = next(item for item in plan(meta["config"])["models"] if item["id"] == model_id)
    verify_tools(meta)
    if spec["family"] == "gnomix":
        gnomix_classes(meta["config"])
    base = load(base_folder / "base.pkl")
    if (base["base_id"] != spec["base_id"] or base["panel_sha256"] != digest(panel / "data.npz")
            or digest(base_folder / "base.pkl") != read(base_folder / "base.json")["sha256"]):
        raise ValueError("Base model/panel provenance differs")
    data = np.load(panel / "data.npz")
    if spec["family"] == "gnomix":
        gnomix_classes(meta["config"])
        model, grid = base["model"], base["layout"]
        model.smooth.model.CRF.set_params(c2=spec["params"].get("crf_c2", 1))
        labels = data["train_labels"].reshape(-1, len(meta["keys"]))
        window_labels = np.asarray([
            [np.bincount(row[a:b], minlength=len(meta["sources"])).argmax()
             for a, b in zip(grid["starts"], grid["ends"])] for row in labels])
        if set(window_labels.flat) != set(range(len(meta["sources"]))):
            raise ValueError("Training-only CRF labels do not cover every source")
        model.smooth.train(model.base.predict_proba(data["train"].reshape(-1, len(meta["keys"]))), window_labels)
        base["model"] = model
    elif spec["family"] in ("pca_ridge", "pca_rf", "pca_xgboost"):
        head = regression(spec["family"], spec["params"], meta["config"]["seed"])
        head.fit(transform(base, data["train"]), data["truth"][:len(meta["training"])])
        base["head"] = head
    base["spec"] = spec
    dump(output / "model.pkl", base)
    write(output / "model.json", {
        "id": model_id, "family": spec["family"], "base_id": spec["base_id"],
        "fitting_roles": ["reference", "train"] if spec["family"] in ("gnomix", "pca_ridge", "pca_rf", "pca_xgboost") else ["reference"],
        "license_track": "research_only" if spec["family"] == "gnomix" else "permissive_software",
        "sha256": digest(output / "model.pkl"), "panel_sha256": base["panel_sha256"],
        "effective_crf_c2": base["model"].smooth.model.CRF.get_params()["c2"] if spec["family"] == "gnomix" else None,
    })


def flare_posterior(path, expected_keys, components, span, sample_id=None):
    from core import physical_weights
    weights = physical_weights([key[1] for key in expected_keys], span)
    value = np.zeros(len(components))
    count, maximum_correction, header_seen = 0, 0., False
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM") and sample_id is not None:
                if line.rstrip().split("\t")[9:] != [sample_id]:
                    raise ValueError("FLARE changed the singleton sample identity")
                header_seen = True
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if count >= len(expected_keys) or key != tuple(expected_keys[count]) or len(fields) != 10:
                raise ValueError("FLARE changed marker coverage/order or singleton sample axis")
            names, calls = fields[8].split(":"), fields[9].split(":")
            haplotypes = []
            for column in ("ANP1", "ANP2"):
                raw = np.asarray(calls[names.index(column)].split(","), dtype=float)
                rounded_flare_probability(raw, len(components))
                maximum_correction = max(maximum_correction, abs(raw.sum() - 1))
                haplotypes.append(raw / raw.sum())
            value += weights[count] * np.mean(haplotypes, axis=0)
            count += 1
    if count != len(expected_keys) or (sample_id is not None and not header_seen):
        raise ValueError("FLARE posterior is missing required markers")
    return probability(value / (span[1] - span[0])), maximum_correction


def flare_command(meta, params, base_folder, target, prefix):
    config = meta["config"]
    command = [
        "java", java_heap(), "-jar", config["tools"]["flare_jar"], f"gt={target}",
        f"ref={base_folder / 'reference.vcf.gz'}", f"ref-panel={base_folder / 'reference.panel'}",
        f"map={config['maps'][str(meta['chromosome'])]}", f"out={prefix}",
        "probs=true", "nthreads=1", f"seed={config['seed']}", "min-maf=0", "min-mac=0"]
    for name in ("array", "em", "gen", "update-p", "states", "panel-weight"):
        if name in params:
            command.append(f"{name}={str(params[name]).lower()}")
    if "model_file" in params:
        command.append(f"model={params['model_file']}")
    return command


def infer(panel, phased, base_folder, fitted, output, *, query_start=0, query_end=None, reuse_manifest=None):
    from core import project_simplex
    output.mkdir()
    meta = read(panel / "panel.json")
    spec = read(fitted / "model.json")
    chunked = query_end is not None
    if query_end is None:
        if query_start != 0:
            raise ValueError("A query chunk requires both range endpoints")
        query_end = len(meta["queries"])
    if (type(query_start) is not int or type(query_end) is not int
            or not 0 <= query_start < query_end <= len(meta["queries"])):
        raise ValueError("Inference query range must be nonempty and inside the declared person axis")
    if (chunked or query_start != 0) and spec["family"] != "flare":
        raise ValueError("Query chunks are supported only for singleton FLARE inference")
    samples = meta["queries"][query_start:query_end]
    verify_tools(meta)
    if spec["family"] == "gnomix":
        gnomix_classes(meta["config"])
    if spec["sha256"] != digest(fitted / "model.pkl"):
        raise ValueError("Saved fitted model changed")
    if digest(phased / "query.npy") != read(phased / "phase.json")["query_sha256"]:
        raise ValueError("Phased query changed")
    model = load(fitted / "model.pkl")
    base_record = read(base_folder / "base.json")
    if base_record["base_id"] != model["base_id"] or base_record["panel_sha256"] != model["panel_sha256"]:
        raise ValueError("Inference reference base differs from fitted model")
    for filename, expected in base_record["artifact_sha256"].items():
        if digest(base_folder / filename) != expected:
            raise ValueError(f"Saved reference-base artifact changed: {filename}")
    if model["panel_sha256"] != digest(panel / "data.npz"):
        raise ValueError("Inference panel differs from fit panel")
    query = np.load(phased / "query.npy", mmap_mode="r")[query_start:query_end]
    config, family = meta["config"], spec["family"]
    params = model["spec"]["params"]
    started = time.monotonic()
    components = meta["sources"]
    corrections = []
    reused = []
    if family == "flare":
        from reuse import copy_native, inference_inputs, native_job
        reusable = native_job(reuse_manifest, meta, spec, inference_inputs(panel, phased, base_folder, fitted),
                              panel=panel, params=params)
        components = params.get("components", meta["sources"])
        values = []
        for index, sample in enumerate(samples):
            target = output / f"{sample['sample_id']}.vcf.gz"
            prefix = output / sample["sample_id"]
            if reusable is not None and sample["sample_id"] in reusable["samples"]:
                copy_native(reusable["samples"][sample["sample_id"]], output, sample["sample_id"], meta, query[index])
                reused.append(sample["sample_id"])
            else:
                write_vcf(target, meta["keys"], [sample["sample_id"]], query[index:index + 1])
                run(flare_command(meta, params, base_folder, target, prefix), output / f"{sample['sample_id']}.flare")
            names = [line.split() for line in Path(str(prefix) + ".model").read_text().splitlines()
                     if line.strip() and not line.lstrip().startswith("#")][0]
            if len(names) != len(components) or set(names) != set(components):
                raise ValueError(f"FLARE ancestry names differ from the declared component axis: {names}")
            value, correction = flare_posterior(str(prefix) + ".anc.vcf.gz", meta["keys"], components,
                                                meta["span"], sample["sample_id"])
            values.append(value[[names.index(component) for component in components]])
            corrections.append(correction)
        values = np.asarray(values)
    elif family == "gnomix":
        gnomix_classes(config)
        grid, learner = model["layout"], model["model"]
        local = np.asarray(learner.smooth.predict_proba(
            learner.base.predict_proba(query.reshape(-1, len(meta["keys"])))), dtype=float)
        classes = list(map(int, learner.smooth.model.classes_))
        if sorted(classes) != list(range(len(components))):
            raise ValueError("Gnomix changed its source class axis")
        local = probability(local[:, :, np.argsort(classes)], len(components), tolerance=2e-6)
        local = local / local.sum(axis=-1, keepdims=True)
        np.savez_compressed(output / "local_posteriors.npz", probabilities=local, weights=grid["weights"],
                            starts=grid["starts"], ends=grid["ends"])
        values = np.einsum("nhwa,w->na", local.reshape(-1, 2, grid["W"], len(components)),
                           grid["weights"]) / (2 * (meta["span"][1] - meta["span"][0]))
    elif family == "nnls":
        frequencies = 2 * model["frequencies"]
        values = []
        for calls in query:
            dosage = calls.sum(axis=0)
            result = minimize(lambda p: float(np.mean((frequencies @ p - dosage)**2)),
                              np.ones(len(components)) / len(components), method="SLSQP",
                              bounds=[(0, 1)] * len(components),
                              constraints={"type": "eq", "fun": lambda p: p.sum() - 1},
                              options={"ftol": 1e-12, "maxiter": 1000})
            if not result.success:
                raise ValueError(f"Constrained frequency control did not converge: {result.message}")
            values.append(probability(result.x, len(components)))
        values = np.asarray(values)
    elif family == "pca_gmm":
        values = probability(model["gmm"].predict_proba(transform(model, query)) @ model["component_map"])
    else:
        values = project_simplex(model["head"].predict(transform(model, query)))
    probability(values, len(components))
    chunk_record = {}
    if chunked:
        from reuse import inference_inputs
        chunk_record = {
            "inference_chunk": {"start": query_start, "end": query_end, "total": len(meta["queries"])},
            "inference_input_sha256": inference_inputs(panel, phased, base_folder, fitted),
            "artifact_sha256": {path.name: digest(path) for path in sorted(output.iterdir()) if path.is_file()},
        }
    np.savez_compressed(output / "predictions.npz", values=values,
                        truth=np.load(panel / "data.npz")["truth"][
                            len(meta["training"]) + query_start:len(meta["training"]) + query_end])
    write(output / "prediction.json", {
        **spec, "chromosome": meta["chromosome"], "route": meta["route"], "span": meta["span"],
        "components": components, "sources": meta["sources"], "samples": samples,
        "seconds": time.monotonic() - started, "max_native_rounding_correction": max(corrections, default=0),
        "native_probability_policy": (
            {"export_component_half_width": .005, "require_feasible_unit_mass": True,
             "integration": "explicit_mass_renormalization_of_rounded_export",
             "exact_unrounded_probabilities_recoverable": False} if family == "flare" else None),
        "sha256": digest(output / "predictions.npz"), "config": config, **chunk_record,
        "reused_singletons": reused, "reuse_manifest_sha256": digest(reuse_manifest) if reuse_manifest else None,
    })


def collect_inference(panel, folders, output):
    output.mkdir()
    if not folders:
        raise ValueError("Require completed FLARE inference chunks")
    meta = read(panel / "panel.json")
    records = sorted(((folder, read(folder / "prediction.json")) for folder in folders),
                     key=lambda item: item[1]["inference_chunk"]["start"])
    first = records[0][1]
    if any(first[key] != meta[key] for key in ("chromosome", "route", "span", "sources", "config")):
        raise ValueError("Inference chunk identity differs from the original panel")
    inputs = first["inference_input_sha256"]
    if inputs["panel.json"] != digest(panel / "panel.json") or inputs["data.npz"] != digest(panel / "data.npz"):
        raise ValueError("Inference chunks belong to a different panel")
    truth = np.load(panel / "data.npz")["truth"][len(meta["training"]):]
    cursor, values, artifacts, chunks = 0, [], {}, []
    shared = ("id", "family", "license_track", "chromosome", "route", "span", "components", "sources",
              "config", "inference_input_sha256", "native_probability_policy")
    for folder, record in records:
        bounds = record["inference_chunk"]
        start, end = bounds["start"], bounds["end"]
        if (record["family"] != "flare" or type(start) is not int or type(end) is not int
                or start != cursor or not start < end <= len(meta["queries"])
                or bounds["total"] != len(meta["queries"])):
            raise ValueError("Incomplete, overlapping or duplicate inference query chunks")
        if any(record[key] != first[key] for key in shared):
            raise ValueError("Inference chunk model, input or feature-axis mismatch")
        if record.get("reuse_manifest_sha256") != first.get("reuse_manifest_sha256"):
            raise ValueError("Inference chunks use different reuse manifests")
        if record["samples"] != meta["queries"][start:end]:
            raise ValueError("Inference chunk person identity/order mismatch")
        if record["sha256"] != digest(folder / "predictions.npz"):
            raise ValueError("Inference chunk predictions changed")
        data = np.load(folder / "predictions.npz")
        value = probability(data["values"], len(record["components"]))
        if value.shape != (end - start, len(record["components"])):
            raise ValueError("Inference chunk prediction person axis differs")
        if not np.array_equal(data["truth"], truth[start:end]):
            raise ValueError("Inference chunk truth differs from the original panel")
        required = {f"{sample['sample_id']}{suffix}" for sample in record["samples"]
                    for suffix in (".model", ".anc.vcf.gz")}
        if not required <= set(record["artifact_sha256"]):
            raise ValueError("Inference chunk is missing singleton FLARE evidence")
        for name, expected in record["artifact_sha256"].items():
            if Path(name).name != name or name in artifacts or name in ("prediction.json", "predictions.npz"):
                raise ValueError("Duplicate or invalid inference chunk artifact name")
            if digest(folder / name) != expected:
                raise ValueError(f"Inference chunk native artifact changed: {name}")
            shutil.copyfile(folder / name, output / name)
            artifacts[name] = expected
        chunks.append({**bounds, "metadata_sha256": digest(folder / "prediction.json"),
                       "predictions_sha256": record["sha256"]})
        values.append(value)
        cursor = end
    if cursor != len(meta["queries"]):
        raise ValueError("Incomplete inference query chunks")
    np.savez_compressed(output / "predictions.npz", values=np.concatenate(values), truth=truth)
    write(output / "prediction.json", {
        **{key: value for key, value in first.items() if key != "inference_chunk"},
        "samples": meta["queries"], "sha256": digest(output / "predictions.npz"),
        "seconds": sum(record["seconds"] for _, record in records),
        "seconds_semantics": "sum_of_chunk_inference_seconds_not_elapsed_wall_time",
        "max_native_rounding_correction": max(record["max_native_rounding_correction"] for _, record in records),
        "artifact_sha256": artifacts, "inference_chunks": chunks,
        "reused_singletons": [sample for _, record in records for sample in record.get("reused_singletons", [])],
    })
