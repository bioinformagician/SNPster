"""Input contracts and common physical-span estimands for ancestry experiments."""
import csv
from collections import defaultdict
import gzip
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np


ROLES = ("reference", "train", "calibration", "validation", "test")
FAMILIES = {
    "flare": {"array", "em", "gen", "update-p", "states", "panel-weight", "model_file", "components", "component_map", "reference_label"},
    "gnomix": {"window_size_cM", "context_ratio", "logistic_c", "crf_c2"},
    "nnls": set(),
    "pca_ridge": {"n_components", "alpha"},
    "pca_rf": {"n_components", "n_estimators", "min_samples_leaf", "max_depth"},
    "pca_xgboost": {"n_components", "n_estimators", "max_depth", "learning_rate"},
    "pca_gmm": {"n_components", "k", "reg_covar"},
}
for _parameters in FAMILIES.values():
    _parameters.add("reference_cap")
FAMILIES["flare"].add("reference_source_cap")
INTERPRETERS = {
    "native": set(), "simplex": {"l2"}, "ridge": {"alpha"},
    "rf": {"n_estimators", "min_samples_leaf", "max_depth"},
    "xgboost": {"n_estimators", "max_depth", "learning_rate"},
}


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    pending = path.with_suffix(path.suffix + ".pending")
    pending.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    pending.replace(path)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            value.update(block)
    return value.hexdigest()


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command, log):
    command = [str(part) for part in command]
    write(str(log) + ".command.json", command)
    with Path(str(log) + ".log").open("w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def probability(value, size=None, tolerance=1e-8):
    value = np.asarray(value, dtype=float)
    if (value.ndim < 1 or (size is not None and value.shape[-1] != size)
            or not np.isfinite(value).all() or (value < 0).any()
            or not np.allclose(value.sum(axis=-1), 1, rtol=0, atol=tolerance)):
        raise ValueError("Invalid ancestry probability axis, bounds or mass")
    return value


def ancestry_tv_pp(value, true, sources, *, combine_eas_sea=False):
    value, true = probability(value, len(sources)), probability(true, len(sources))
    if value.shape != (len(sources),) or true.shape != value.shape:
        raise ValueError("TV requires one matching named source vector per person")
    errors = value - true
    if combine_eas_sea:
        if "EAS" not in sources or "SEA" not in sources:
            raise ValueError("Combined metrics require separate EAS and SEA source axes")
        eas, sea = sources.index("EAS"), sources.index("SEA")
        return float(50 * (np.abs(np.delete(errors, [eas, sea])).sum() + abs(errors[eas] + errors[sea])))
    return float(50 * np.abs(errors).sum())


def rounded_flare_probability(value, size):
    value = np.asarray(value, dtype=float)
    if (value.shape != (size,) or not np.isfinite(value).all()
            or (value < 0).any() or (value > 1).any() or value.sum() <= 0):
        raise ValueError("Invalid rounded FLARE posterior axis or bounds")
    # The pinned writer exports each component to two decimal places.
    lower = np.maximum(value - .005, 0).sum()
    upper = np.minimum(value + .005, 1).sum()
    if lower > 1 + 1e-12 or upper < 1 - 1e-12:
        raise ValueError("FLARE posterior mass is incompatible with per-component export rounding")
    return value


def project_simplex(value):
    """Explicit decoder link; never used to repair native posterior errors."""
    value = np.asarray(value, dtype=float)
    if value.ndim != 2 or not np.isfinite(value).all():
        raise ValueError("Decoder returned invalid regression output")
    ordered = np.sort(value, axis=1)[:, ::-1]
    cumulative = np.cumsum(ordered, axis=1) - 1
    rank = np.arange(1, value.shape[1] + 1)
    rho = (ordered - cumulative / rank > 0).sum(axis=1) - 1
    threshold = cumulative[np.arange(len(value)), rho] / (rho + 1)
    return probability(np.maximum(value - threshold[:, None], 0))


def physical_weights(positions, span):
    positions = np.asarray(positions, dtype=np.int64)
    if (len(positions) < 2 or np.any(np.diff(positions) <= 0)
            or not span[0] <= positions[0] or not positions[-1] < span[1]):
        raise ValueError("Markers must be unique, ordered and inside a half-open physical span")
    weights = np.diff(np.r_[span[0], (positions[:-1].astype(float) + positions[1:]) / 2, span[1]])
    if (weights <= 0).any() or weights.sum() != span[1] - span[0]:
        raise ValueError("Physical posterior cells do not cover the chromosome span")
    return weights


def validate_roles(samples, sources, *, development_only=False):
    required = {"sample_id", "role", "family_id", "donors", "source", "track", "fine"}
    if not samples or any(not required <= set(row) for row in samples):
        raise ValueError(f"Sample manifest requires {sorted(required)}")
    ids = [row["sample_id"] for row in samples]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate sample identities")
    ownership = {}
    for row in samples:
        if row["role"] not in ROLES or not row["family_id"] or not row["donors"]:
            raise ValueError("Every sample needs an explicit role, family and founder donors")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", row["sample_id"]):
            raise ValueError("Unsupported sample ID")
        if row["role"] == "reference" and row["source"] not in sources:
            raise ValueError("Reference source label is unknown")
        for token in [("family", row["family_id"]), ("donor", row["sample_id"])] + [
                ("donor", donor) for donor in row["donors"].split(",")]:
            if not token[1] or (token in ownership and ownership[token] != row["role"]):
                raise ValueError(f"Cross-role founder/family leakage: {token}")
            ownership[token] = row["role"]
    expected_roles = set(ROLES) - ({"test"} if development_only else set())
    if set(row["role"] for row in samples) != expected_roles:
        raise ValueError(f"Require exactly the disjoint roles {sorted(expected_roles)}")
    if {row["source"] for row in samples if row["role"] == "reference"} != set(sources):
        raise ValueError("Reference must cover every declared source")


def validate_config(config):
    required = {
        "version", "seed", "mode", "sources", "chromosomes", "spans", "samples", "segments",
        "reference", "training", "maps", "routes", "models", "interpretations", "tools",
        "reference_filters", "allow_research_only", "independent_test",
    }
    if set(config) - {"combine_eas_sea_for_metrics"} != required or config["version"] not in (1, 2):
        raise ValueError(f"Config keys must be exactly {sorted(required)}; version must be 1 or 2")
    if type(config["seed"]) is not int or not 0 <= config["seed"] < 2**31:
        raise ValueError("seed must be an unsigned 31-bit integer")
    if any(type(config[key]) is not bool for key in ("allow_research_only", "independent_test")):
        raise ValueError("Research/independent-test declarations must be booleans")
    if (len(config["sources"]) < 2 or len(set(config["sources"])) != len(config["sources"])
            or any(not re.fullmatch(r"[A-Za-z0-9_]+", item) for item in config["sources"])):
        raise ValueError("Require unique, explicitly named ancestry sources")
    combined = config.get("combine_eas_sea_for_metrics", False)
    if type(combined) is not bool or (combined and not {"EAS", "SEA"} <= set(config["sources"])):
        raise ValueError("combine_eas_sea_for_metrics must be boolean and requires separate EAS/SEA axes")
    chromosomes = config["chromosomes"]
    if (not chromosomes or any(type(chrom) is not int for chrom in chromosomes)
            or len(set(chromosomes)) != len(chromosomes) or set(chromosomes) - set(range(1, 23))):
        raise ValueError("Require unique autosomes")
    modes = ("smoke", "full") if config["version"] == 1 else ("development_screen", "development_full")
    if config["mode"] not in modes:
        raise ValueError(f"mode must be one of {modes} for this config version")
    if config["mode"] in ("full", "development_full") and chromosomes != list(range(1, 23)):
        raise ValueError("Partial-autosome full scores are forbidden")
    if config["mode"] == "smoke" and config["independent_test"]:
        raise ValueError("A smoke cannot certify independent accuracy")
    if config["version"] == 2 and config["independent_test"]:
        raise ValueError("Development excludes the test role and cannot certify independent accuracy")
    for chrom in chromosomes:
        for key in ("reference", "training", "maps", "spans"):
            if str(chrom) not in config[key]:
                raise ValueError(f"Missing {key}/chr{chrom}")
        span = config["spans"][str(chrom)]
        if len(span) != 2 or any(type(value) is not int for value in span) or not 0 <= span[0] < span[1]:
            raise ValueError("Invalid half-open chromosome span")
    for group in ("models", "interpretations", "routes"):
        if not config[group] or len({item["id"] for item in config[group]}) != len(config[group]):
            raise ValueError(f"Empty or duplicate {group}")
        if any(not re.fullmatch(r"[A-Za-z0-9_]+", item["id"]) for item in config[group]):
            raise ValueError(f"Unsafe {group} identifier")
    for model in config["models"]:
        if set(model) != {"id", "family", "params"} or model["family"] not in FAMILIES:
            raise ValueError("Unknown model family or model fields")
        if set(model["params"]) - FAMILIES[model["family"]]:
            raise ValueError(f"Unknown or unsupported model hyperparameters: {model['id']}")
        if model["family"] == "gnomix" and not config["allow_research_only"]:
            raise ValueError("Gnomix requires explicit academic/research-only opt-in")
        validate_parameters(model["params"])
        fine_reference = model["params"].get("reference_label") == "fine"
        if "reference_label" in model["params"] and not fine_reference:
            raise ValueError("reference_label must be 'fine' when explicitly supplied")
        if fine_reference and ("model_file" in model["params"] or not {
                "components", "component_map"} <= set(model["params"])):
            raise ValueError("Fine reference labels require explicit components/mapping and no supplied model file")
        if "components" in model["params"]:
            components = model["params"]["components"]
            if not components or len(set(components)) != len(components):
                raise ValueError("Native components must have a unique, stable declared order")
            if "model_file" not in model["params"] and not fine_reference:
                raise ValueError("Custom FLARE components require a supplied native model file")
        if "component_map" in model["params"]:
            mapping = probability(model["params"]["component_map"], len(config["sources"]))
            if mapping.shape != (len(model["params"].get("components", config["sources"])), len(config["sources"])):
                raise ValueError("Explicit component map shape differs from declared components/sources")
    signatures = [
        json.dumps({"family": model["family"], "params": {
            key: value for key, value in model["params"].items()
            if not (model["family"] == "flare" and key == "array")}}, sort_keys=True)
        for model in config["models"]]
    if len(signatures) != len(set(signatures)):
        raise ValueError("Duplicate model settings (FLARE array is an alias with the common explicit MAC/MAF policy)")
    for item in config["interpretations"]:
        if set(item) != {"id", "kind", "params"} or item["kind"] not in INTERPRETERS:
            raise ValueError("Unknown interpretation")
        if set(item["params"]) - INTERPRETERS[item["kind"]]:
            raise ValueError("Unknown interpretation hyperparameters")
        validate_parameters(item["params"])
    if not any(item["kind"] == "native" for item in config["interpretations"]):
        raise ValueError("The uncalibrated native track must always be retained")
    for route in config["routes"]:
        if set(route) != {"id", "stage", "vcfs"} or route["stage"] not in ("harmonized", "imputed"):
            raise ValueError("Route stage must be harmonized or imputed")
        if set(route["vcfs"]) != set(map(str, chromosomes)):
            raise ValueError("Every route must cover the exact declared chromosome scope")
    filters = config["reference_filters"]
    if (set(filters) != {"min_maf", "min_mac", "minimum_markers"}
            or not 0 <= filters["min_maf"] < .5 or filters["min_mac"] < 0
            or filters["minimum_markers"] < 2):
        raise ValueError("Invalid common reference-only marker filters")
    if set(config["tools"]) - {"flare_jar", "beagle_jar", "gnomix_root"}:
        raise ValueError("Unknown tool; declare pinned FLARE/Beagle jars or Gnomix checkout")
    if any(row["family"] == "flare" for row in config["models"]) and "flare_jar" not in config["tools"]:
        raise ValueError("FLARE models require tools.flare_jar")
    if any(row["family"] == "gnomix" for row in config["models"]) and "gnomix_root" not in config["tools"]:
        raise ValueError("Gnomix models require tools.gnomix_root")


def validate_parameters(params):
    integer = {"reference_cap", "reference_source_cap", "n_components", "k", "n_estimators", "min_samples_leaf", "max_depth", "states"}
    nonnegative = {"crf_c2", "l2", "alpha", "reg_covar"}
    positive = {"gen", "window_size_cM", "logistic_c", "learning_rate"}
    for key, value in params.items():
        if key in integer and (type(value) is not int or value < 1):
            raise ValueError(f"{key} must be a positive integer")
        if key in nonnegative | positive and (
                type(value) not in (int, float) or not np.isfinite(value)
                or value < 0 or (key in positive and value == 0)):
            raise ValueError(f"Invalid numeric hyperparameter: {key}")
        if key in {"em", "array", "update-p"} and type(value) is not bool:
            raise ValueError(f"{key} must be boolean")
        if key in {"context_ratio", "panel-weight"} and (
                type(value) not in (float, int) or not np.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"{key} must lie in [0,1]")


def java_heap():
    value = int(os.environ.get("SNPSTER_JAVA_HEAP_MB", "6144"))
    if value < 256:
        raise ValueError("Task memory is too small for Java")
    return f"-Xmx{value}m"


def input_paths(config):
    paths = [config["samples"], config["segments"]]
    for name in ("reference", "training", "maps"):
        paths += list(config[name].values())
    for route in config["routes"]:
        paths += list(route["vcfs"].values())
    paths += [value for key, value in config["tools"].items() if key != "gnomix_root"]
    if any(model["family"] == "gnomix" for model in config["models"]):
        root = Path(config["tools"]["gnomix_root"])
        paths += [str(path) for path in sorted(root.rglob("*.py"))]
        paths += [str(root / "LICENSE.md")]
    for model in config["models"]:
        if "model_file" in model["params"]:
            paths.append(model["params"]["model_file"])
    return sorted(set(paths))


def freeze(config_path, output):
    config = read(config_path)
    validate_config(config)
    samples = rows(config["samples"])
    validate_roles(samples, config["sources"], development_only=config["version"] == 2)
    hashes = {}
    for path in input_paths(config):
        if not Path(path).is_absolute() or not Path(path).is_file():
            raise ValueError(f"Inputs must be existing absolute file paths: {path}")
        hashes[path] = digest(path)
    for chrom in config["chromosomes"]:
        truth(config, samples, chrom)
    output.mkdir()
    write(output / "config.json", config)
    write(output / "manifest.json", {
        "config_sha256": digest(config_path), "input_sha256": hashes,
        "scope": config["mode"], "declared_role_disjointness": True,
        "genetic_relatedness_independently_certified": False,
        "production_promotion": False, "pgs_changed": False,
    })


def decode_gt(text, sample_count, *, phased=False, allow_missing=False):
    raw = (text + "\t").encode("ascii")
    if len(raw) != 4 * sample_count and allow_missing:
        raw = ("\t".join("./." if token == "." else token for token in text.split("\t")) + "\t").encode("ascii")
    allowed = b"01." if allow_missing else b"01"
    alleles = raw[::2]
    if (len(raw) != 4 * sample_count or alleles.translate(None, allowed)
            or raw[1::4].translate(None, b"|/") or raw[3::4].translate(None, b"\t")):
        raise ValueError("Missing/non-biallelic GT; filling/DS substitution forbidden")
    if phased and raw[1::4].translate(None, b"|"):
        raise ValueError("Unphased GT supplied to a phase-dependent model")
    table = bytes.maketrans(b".01", bytes((255, 0, 1)))
    return np.frombuffer(alleles.translate(table), dtype=np.int8).reshape(sample_count, 2)


def vcf(path, chromosome, identities, phased):
    header = subprocess.check_output(["bcftools", "query", "-l", str(path)], text=True).splitlines()
    if len(set(header)) != len(header) or not set(identities) <= set(header):
        raise ValueError(f"VCF identity mismatch: {path}")
    order = [header.index(identity) for identity in identities]
    # Stream only the requested samples; bcftools preserves the requested sample order.
    command = ["bcftools", "query", "-s", ",".join(identities),
               "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", str(path)]
    keys, values = [], []
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t", 4)
            if fields[0] != str(chromosome):
                continue
            if len(fields[2]) != 1 or len(fields[3]) != 1 or fields[2] not in "ACGT" or fields[3] not in "ACGT":
                continue
            try:
                calls = decode_gt(fields[4], len(order), phased=phased)
            except ValueError as error:
                raise ValueError(f"{error}: {path}:{fields[1]}") from error
            keys.append((fields[0], int(fields[1]), fields[2], fields[3]))
            values.append(calls)
        if process.wait():
            raise RuntimeError(f"bcftools query failed: {path}")
    if not keys or len(set(keys)) != len(keys) or any(
            keys[index][1] <= keys[index - 1][1] for index in range(1, len(keys))):
        raise ValueError(f"Empty, ambiguous or unordered SNP panel: {path}")
    return keys, np.asarray(values, dtype=np.int8).transpose(1, 2, 0)


def write_vcf(path, keys, identities, calls, phased=True):
    if calls.shape != (len(identities), 2, len(keys)):
        raise ValueError("VCF writing axes differ")
    if not calls.size or calls.dtype.kind not in "biu" or calls.min() < 0 or calls.max() > 1:
        raise ValueError("Model VCF writing requires complete biallelic GT")
    sep = "|" if phased else "/"
    encoded = np.asarray([f"{a}{sep}{b}\t".encode("ascii") for a in (0, 1) for b in (0, 1)], dtype="S4")
    raw = path.with_suffix("")
    with raw.open("wb") as handle:
        handle.write((f"##fileformat=VCFv4.2\n##contig=<ID={keys[0][0]}>\n"
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(identities) + "\n").encode("utf-8"))
        for index, (chrom, pos, ref, alt) in enumerate(keys):
            gt = encoded[2 * calls[:, 0, index] + calls[:, 1, index]].tobytes()
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t".encode("ascii"))
            handle.write(gt[:-1] + b"\n")
    run(["bcftools", "view", "-Oz", "-o", path, raw], path.with_suffix(".compress"))
    run(["bcftools", "index", path], path.with_suffix(".index"))
    raw.unlink()


def truth(config, samples, chrom, positions=None):
    segments = defaultdict(list)
    with Path(config["segments"]).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["chromosome"] == str(chrom):
                segments[(row["sample_id"], row["haplotype"])].append(
                    (int(row["start"]), int(row["end"]), row["source"]))
    span = config["spans"][str(chrom)]
    fractions, labels = [], []
    for sample in samples:
        if sample["role"] == "reference":
            continue
        counts = np.zeros(len(config["sources"]))
        haplabels = []
        for hap in (0, 1):
            selected = sorted(segments[(sample["sample_id"], str(hap))])
            cursor = span[0]
            result = np.full(len(positions) if positions is not None else 0, -1, dtype=np.int16)
            for start, end, source in selected:
                if start != cursor or not start < end <= span[1] or source not in config["sources"]:
                    raise ValueError(f"Invalid/gapped/overlapping physical truth: {sample['sample_id']}/{chrom}/{hap}")
                index = config["sources"].index(source)
                counts[index] += end - start
                if positions is not None:
                    result[(positions >= start) & (positions < end)] = index
                cursor = end
            if cursor != span[1]:
                raise ValueError(f"Incomplete physical truth: {sample['sample_id']}/{chrom}/{hap}")
            haplabels.append(result)
        fractions.append(counts / (2 * (span[1] - span[0])))
        labels.append(haplabels)
    return probability(fractions, len(config["sources"])), np.asarray(labels)


def prepare(frozen, route_id, chrom, output):
    from models import plan
    config = read(frozen / "config.json")
    manifest = read(frozen / "manifest.json")
    route = next(row for row in config["routes"] if row["id"] == route_id)
    samples = rows(config["samples"])
    reference = [row for row in samples if row["role"] == "reference"]
    training = [row for row in samples if row["role"] == "train"]
    queries = [row for row in samples if row["role"] in ("calibration", "validation", "test")]
    for path in [config["reference"][str(chrom)], config["training"][str(chrom)],
                 route["vcfs"][str(chrom)], config["samples"], config["segments"], config["maps"][str(chrom)]]:
        if digest(path) != manifest["input_sha256"][path]:
            raise ValueError(f"Frozen scientific input changed: {path}")
    output.mkdir()
    refkeys, ref = vcf(config["reference"][str(chrom)], chrom, [row["sample_id"] for row in reference], True)
    trainkeys, train = vcf(config["training"][str(chrom)], chrom, [row["sample_id"] for row in training], True)
    querykeys, query = vcf(route["vcfs"][str(chrom)], chrom, [row["sample_id"] for row in queries], route["stage"] == "imputed")
    span = config["spans"][str(chrom)]
    common = set(trainkeys) & set(querykeys)
    counts = ref.sum(axis=(0, 1))
    mac = np.minimum(counts, 2 * len(reference) - counts)
    filters = config["reference_filters"]
    keys = [key for index, key in enumerate(refkeys)
            if key in common and span[0] <= key[1] < span[1]
            and mac[index] >= filters["min_mac"]
            and mac[index] / (2 * len(reference)) >= filters["min_maf"]]
    if len(keys) < filters["minimum_markers"]:
        raise ValueError(f"Too few exact common reference-qualified markers: {len(keys)}")
    def align(original, calls):
        offsets = {key: index for index, key in enumerate(original)}
        return calls[:, :, [offsets[key] for key in keys]]
    ref, train, query = align(refkeys, ref), align(trainkeys, train), align(querykeys, query)
    weights = physical_weights([key[1] for key in keys], span)
    targets, labels = truth(config, training + queries, chrom, np.asarray([key[1] for key in keys]))
    np.savez_compressed(output / "data.npz", reference=ref, train=train, query=query,
                        train_labels=labels[:len(training)], truth=targets,
                        weights=weights, keys=np.asarray(keys, dtype=str))
    write_vcf(output / "reference.vcf.gz", keys, [row["sample_id"] for row in reference], ref)
    (output / "reference.panel").write_text(
        "".join(f"{row['sample_id']}\t{row['source']}\n" for row in reference))
    write(output / "panel.json", {
        "route": route, "chromosome": chrom, "span": span, "reference": reference,
        "training": training, "queries": queries, "keys": keys, "sources": config["sources"],
        "config": config, "data_sha256": digest(output / "data.npz"),
        "available_reference_markers": len(refkeys), "retained_markers": len(keys),
        "common_quality_policy": filters,
        "base_plan": plan(config),
        "tool_sha256": {path: value for path, value in manifest["input_sha256"].items()
                        if path not in [config["samples"], config["segments"]]
                        and path not in list(config["reference"].values()) + list(config["training"].values())
                        + list(config["maps"].values()) + [value for item in config["routes"] for value in item["vcfs"].values()]},
    })


def phase(panel, output):
    output.mkdir()
    meta = read(panel / "panel.json")
    data = np.load(panel / "data.npz")
    query = data["query"]
    keys = [tuple(key) for key in meta["keys"]]
    config = meta["config"]
    requires_phase = any(model["family"] in ("flare", "gnomix") for model in config["models"])
    verify_tools(meta)
    if meta["route"]["stage"] == "harmonized" and requires_phase:
        if "beagle_jar" not in config["tools"]:
            raise ValueError("Harmonized FLARE/Gnomix input requires the declared phase-only Beagle jar")
        query = query.copy()
        for index, sample in enumerate(meta["queries"]):
            target = output / f"{sample['sample_id']}.vcf.gz"
            write_vcf(target, keys, [sample["sample_id"]], data["query"][index:index + 1], phased=False)
            prefix = output / f"{sample['sample_id']}.phase"
            run(["java", java_heap(), "-jar", config["tools"]["beagle_jar"], f"gt={target}",
                 f"ref={panel / 'reference.vcf.gz'}", f"map={config['maps'][str(meta['chromosome'])]}",
                 f"out={prefix}", "impute=false", "nthreads=1", f"seed={config['seed']}"], prefix)
            actualkeys, calls = vcf(str(prefix) + ".vcf.gz", meta["chromosome"], [sample["sample_id"]], True)
            if actualkeys != keys or not np.array_equal(calls.sum(axis=1), data["query"][index:index + 1].sum(axis=1)):
                raise ValueError("Phase-only adapter changed observed keys or diploid GT")
            query[index] = calls[0]
    np.save(output / "query.npy", query)
    write(output / "phase.json", {
        "stage": meta["route"]["stage"], "phase_only": meta["route"]["stage"] == "harmonized" and requires_phase,
        "added_variants": 0, "changed_diploid_calls": 0, "query_sha256": digest(output / "query.npy"),
    })


def verify_tools(meta):
    for path, expected in meta["tool_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"Frozen native software/model changed: {path}")
