#!/usr/bin/env python3
"""Known-haplotype, identity-disjoint singleton array/imputation transfer pilot."""

import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import model_family_gnomix as adapter
from realized_ancestry import load_bp
from screen_model_families import verify


CANDIDATE = "gnomix_fast_050"
PROJECT = Path(__file__).resolve().parents[2]
REFERENCE_PATTERN = "hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.SNV_biallelic.numericCHR.vcf.gz"


def checked(command, output, label):
    command = list(map(str, command))
    adapter.write_json(output / f"{label}.command.json", command)
    started = time.monotonic()
    with (output / f"{label}.log").open("w") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
    return time.monotonic() - started


def fingerprint_files(paths):
    return {str(path): adapter.digest(path) for path in sorted(set(map(Path, paths)))}


def verify_hashes(hashes):
    for path, expected in hashes.items():
        if adapter.digest(path) != expected:
            raise ValueError(f"Pilot provenance changed: {path}")


def ancestor_donors(plan, sample):
    nodes = {row["sample_id"]: row for row in plan["nodes"]}
    pending, seen, donors = [sample], set(), set()
    while pending:
        person = pending.pop()
        if person in ("", "0") or person in seen:
            continue
        seen.add(person)
        if person in plan["founders"]:
            donors.add(plan["founders"][person]["vcf_sample"])
        else:
            node = nodes[person]
            pending.extend([node["father_id"], node["mother_id"]])
    if not donors:
        raise ValueError("No original founder donors found")
    return donors


def pedigree_haplotypes(plan, bp, sample, donor_ids, donors, positions):
    hap_to_donor = {}
    for founder, info in plan["founders"].items():
        if info["vcf_sample"] not in donor_ids:
            continue
        for hap_index, name in enumerate(("h0", "h1")):
            values = {segment[2] for segment in bp[founder][name]}
            if len(values) != 1 or next(iter(values)) in hap_to_donor:
                raise ValueError("Ambiguous founder haplotype mapping")
            hap_to_donor[next(iter(values))] = (donor_ids.index(info["vcf_sample"]), hap_index)
    output = np.empty((2, len(positions)), dtype=np.int8)
    for hap, name in enumerate(("h0", "h1")):
        covered = np.zeros(len(positions), dtype=np.int8)
        for start, inclusive_end, identity in bp[sample][name]:
            if identity not in hap_to_donor:
                raise ValueError("An inherited donor is absent from the dense founder panel")
            first, end = np.searchsorted(positions, [start, inclusive_end + 1])
            donor, allele = hap_to_donor[identity]
            output[hap, first:end] = donors[2 * donor + allele, first:end]
            covered[first:end] += 1
        if not np.all(covered == 1):
            raise ValueError("Inherited haplotypes do not cover the requested physical span exactly once")
    return output


def phase_metrics(truth, inferred):
    if (truth.shape != inferred.shape or truth.ndim != 2 or truth.shape[0] != 2
            or np.any(truth < 0) or np.any(truth > 1)
            or np.any(inferred < 0) or np.any(inferred > 1)):
        raise ValueError("Phase comparison requires complete paired biallelic haplotypes")
    td, pdos = truth.sum(axis=0), inferred.sum(axis=0)
    called = (td == 1) & (pdos == 1)
    orientation = truth[0, called] != inferred[0, called]
    n = len(orientation)
    switches = int(np.count_nonzero(orientation[1:] != orientation[:-1]))
    alternatives = int(orientation.sum())
    return dict(
        n_sites=len(td), n_genotype_errors=int(np.count_nonzero(td != pdos)),
        n_truth_heterozygotes=int((td == 1).sum()),
        n_correct_heterozygotes=n, switch_pairs=max(n - 1, 0), switches=switches,
        phase_hamming_errors=min(alternatives, n - alternatives),
        dosage_squared_error=int(np.square(td.astype(int) - pdos.astype(int)).sum()))


def genotype_strata(truth, inferred, observed):
    if len(observed) != truth.shape[1]:
        raise ValueError("Observed mask differs from genotype marker axis")
    result = {}
    for label, mask in (("all", np.ones(len(observed), dtype=bool)),
                        ("observed", observed), ("imputed", ~observed)):
        result[label] = phase_metrics(truth[:, mask], inferred[:, mask])
    return result


def prepare(screen, output, reference, beagle):
    if (output / "pilot.json").exists():
        raise ValueError("Pilot preparation is one-shot; select a new output directory")
    output.mkdir(parents=True, exist_ok=True)
    cohort = verify(screen)
    transfer = json.loads((PROJECT / "results_v3_real_transfer_audit_03/audit.json").read_text())
    controls = transfer["matched_50_50_eur_asian_controls"]
    plan_root = PROJECT / "results_v3_pedigree_pilot_10/pedigree"
    samples, batches, donor_ids, inputs = [], {}, set(), []
    for batch in cohort["evaluation"]["pedigree"]["batches"]:
        rows = adapter.records(batch["samples_tsv"])
        selected = [row for row in rows if row["sample_id"] in controls]
        if not selected:
            continue
        name = batch["batch"]
        plan_path = plan_root / f"plans/plans/{name}/pedigrees.json"
        plan = json.loads(plan_path.read_text())
        inputs.extend([plan_path, Path(batch["samples_tsv"]), Path(batch["vcf"])])
        for row in selected:
            if (abs(float(row["realized_EUR"]) - .5) > 1e-10
                    or abs(float(row["realized_EAS"]) + float(row["realized_SEA"]) - .5) > 1e-10):
                raise ValueError("Selected control is no longer a known 50/50 EUR/Asian individual")
            ancestors = ancestor_donors(plan, row["sample_id"])
            donor_ids.update(ancestors)
            samples.append(dict(sample_id=row["sample_id"], batch=name,
                                family_id=row["family_id"], ancestors=sorted(ancestors)))
        batches[name] = dict(plan=str(plan_path), vcf=batch["vcf"])
    if sorted(row["sample_id"] for row in samples) != sorted(controls):
        raise ValueError("Matched controls are missing or duplicated")
    real = Path("/data/flare2_test/source_model/real_user_validation/lonni/lonni.autosomes.vcf.gz")
    observed_file = output / "real_observed_keys.tsv"
    with observed_file.open("w") as stream:
        subprocess.run(["bcftools", "query", "-i", "INFO/IMP!=1",
                        "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(real)], check=True, stdout=stream)
    # Exclude every development founder and pure-validation individual, not just the five controls.
    eligible = plan_root / "inputs/eligible_heldout_founders.tsv"
    eligible_rows = pd.read_csv(eligible, sep="\t")
    if "sample" not in eligible_rows:
        raise ValueError("Unrecognized original eligible-founder metadata")
    excluded = set(eligible_rows["sample"].astype(str)) | donor_ids
    for batch in cohort["evaluation"]["pure"]["batches"]:
        excluded.update(row["sample_id"] for row in adapter.records(batch["samples_tsv"]))
    exclude_path = output / "excluded_reference_ids.txt"
    exclude_path.write_text("".join(sample + "\n" for sample in sorted(excluded)))
    pd.DataFrame(samples).assign(
        ancestors=lambda frame: frame.ancestors.map(lambda values: ",".join(values))).to_csv(
            output / "controls.tsv", sep="\t", index=False)
    inputs.extend([Path(cohort["markers_tsv"]), real, eligible, beagle,
                   Path(adapter.__file__), Path(__file__), screen / "cohort.json",
                   PROJECT / "results_v3_real_transfer_audit_03/audit.json"])
    for chrom in range(1, 23):
        ref = reference / "beagle_references" / REFERENCE_PATTERN.format(chrom=chrom)
        inputs.extend([ref, Path(str(ref) + ".tbi")])
        inputs.append(reference / f"plink.GRCh38.map/plink.chr{chrom}.GRCh38.map")
        model = screen / f"candidates/{CANDIDATE}/models/chr{chrom}.pkl"
        inputs.extend([model, model.with_suffix(".json")])
        for name in batches:
            inputs.append(plan_root / f"chromosomes/pdb_{name}_chr{chrom}.bp")
    audit = dict(
        status="prepared", candidate=CANDIDATE, seed=adapter.SEED,
        screen=str(screen), reference=str(reference), beagle=str(beagle),
        samples=samples, donor_ids=sorted(donor_ids), batches=batches,
        plan_root=str(plan_root), excluded_reference_ids=sorted(excluded),
        original_cohort_sha256=adapter.input_hashes(cohort),
        input_sha256=fingerprint_files(inputs),
        singleton_imputation=True, no_new_families=True,
        new_ancestry_models_trained=False, production_changed=False,
        truth_extension="Original Ped-sim breakpoints applied to the same dense phased founder haplotypes; original frozen-model genotypes must match exactly.",
        unphased_array_input=True, reference_identity_disjoint=True,
        original_historical_beagle_binary_provenance_verified=False,
        runtime_note="Uses Beagle extracted from the currently installed imputer image; historical real-sample binary/seed not recoverable from QC VCF headers.",
        limitations=[
            "Identity-disjoint does not establish kinship-disjoint reference panels.",
            "Five controls share some donors and are not independent replicates.",
            "No real-sample phase or raw-array genotype truth is supplied.",
            "Controls are known 50/50 EUR/Asian, not necessarily EUR/Filipino.",
            "Array SNPs absent from the dense reference or outside saved breakpoint spans cannot be reconstructed.",
            "The observed-site mask is learned from one real sample; this is a transfer diagnostic, not a new holdout.",
        ])
    audit["prepared_output_sha256"] = fingerprint_files(
        [observed_file, exclude_path, output / "controls.tsv"])
    adapter.write_json(output / "pilot.json", audit)
    print(f"Prepared {len(samples)} controls from {len(donor_ids)} original donors", flush=True)


def chromosome(screen, output, chrom):
    pilot = json.loads((output / "pilot.json").read_text())
    cohort = json.loads((screen / "cohort.json").read_text())
    directory = output / f"chr{chrom}"
    checkpoint = directory / "audit.json"
    if checkpoint.exists():
        previous = json.loads(checkpoint.read_text())
        if previous["status"] == "complete":
            verify_hashes(previous["output_sha256"])
            verify_hashes(previous["input_sha256"])
            if previous["pilot_sha256"] != adapter.digest(output / "pilot.json"):
                raise ValueError("Pilot provenance changed since chromosome completion")
            print(f"chr{chrom}: verified completed checkpoint; no repeated imputation", flush=True)
            return
        raise ValueError(f"chr{chrom}: incomplete prior execution; select a new run or inspect explicitly")
    directory.mkdir(parents=True, exist_ok=True)
    adapter.load_upstream()
    reference = Path(pilot["reference"])
    ref = reference / "beagle_references" / REFERENCE_PATTERN.format(chrom=chrom)
    map_path = reference / f"plink.GRCh38.map/plink.chr{chrom}.GRCh38.map"
    plan_root = Path(pilot["plan_root"])
    inputs = [ref, Path(str(ref) + ".tbi"), map_path, Path(pilot["beagle"]),
              Path(__file__), Path(adapter.__file__), screen / "cohort.json"]
    model_path = screen / f"candidates/{CANDIDATE}/models/chr{chrom}.pkl"
    inputs += [model_path, model_path.with_suffix(".json")]
    inputs += [plan_root / f"chromosomes/pdb_{name}_chr{chrom}.bp" for name in pilot["batches"]]
    inputs += [Path(batch[field]) for batch in pilot["batches"].values() for field in ("plan", "vcf")]
    for path in inputs:
        if adapter.digest(path) != pilot["input_sha256"][str(path)]:
            raise ValueError(f"Pilot input changed: {path}")
    verify_hashes(pilot["prepared_output_sha256"])
    bundle, model_metadata = adapter.saved_chromosome(
        cohort, adapter.config_for(CANDIDATE), chrom, screen / "candidates" / CANDIDATE)
    keys, layout = bundle["keys"], bundle["layout"]
    with (output / "real_observed_keys.tsv").open() as stream:
        observed = [(fields[0], int(fields[1]), fields[2], fields[3]) for line in stream
                    if (fields := line.rstrip().split("\t"))[0] == str(chrom)]
    observed_set = set(observed)
    span = cohort["chromosome_spans"][str(chrom)]
    wanted = sorted(set(keys) | {key for key in observed if span[0] <= key[1] < span[1]},
                    key=lambda key: (key[1], key[2], key[3]))
    regions = directory / "truth.positions.tsv"
    regions.write_text("".join(f"{chrom}\t{position}\n" for position in sorted({k[1] for k in wanted})))
    wanted_set, found, seen, donor_values = set(wanted), [], set(), []
    command = ["bcftools", "query", "-R", str(regions), "-s", ",".join(pilot["donor_ids"]),
               "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n", str(ref)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    for line in result.stdout.splitlines():
        row = line.split("\t")
        key = (row[0], int(row[1]), row[2], row[3])
        if key not in wanted_set:
            continue
        values = dict(field.split("=") for field in row[4:])
        if key in seen:
            raise ValueError("Duplicate dense-reference marker identity")
        seen.add(key)
        found.append(key)
        donor_values.append([allele for sample in pilot["donor_ids"]
                             for allele in adapter.phased_gt(values[sample])])
    if not set(keys).issubset(found):
        raise ValueError("Full imputation reference does not cover every frozen-model SNP")
    # Gnomix and Beagle both require unique positions; do not silently choose alternate records.
    if len({key[1] for key in found}) != len(found):
        raise ValueError("Ambiguous multiple allele records at a required array position")
    donors = np.asarray(donor_values, dtype=np.int8).T
    if np.any(donors > 1):
        raise ValueError("Founder truth contains missing genotypes")
    index = {key: i for i, key in enumerate(found)}
    model_indices = np.asarray([index[key] for key in keys])
    positions = np.asarray([key[1] for key in found])
    truth = []
    for sample in pilot["samples"]:
        batch = pilot["batches"][sample["batch"]]
        plan = json.loads(Path(batch["plan"]).read_text())
        bp = load_bp(plan_root / f"chromosomes/pdb_{sample['batch']}_chr{chrom}.bp", chrom)
        x = pedigree_haplotypes(plan, bp, sample["sample_id"], pilot["donor_ids"], donors, positions)
        original = adapter.load_haplotypes(batch["vcf"], [sample["sample_id"]], chrom, keys)
        if not np.array_equal(x[:, model_indices], original):
            raise ValueError("Extending the original pedigree changed its frozen-model haplotypes")
        truth.append(x)
    truth = np.asarray(truth)
    np.save(directory / "truth_haplotypes.npy", truth)
    pd.DataFrame(found, columns=["chrom", "pos", "ref", "alt"]).to_csv(
        directory / "truth_markers.tsv", sep="\t", index=False)
    truth_vcf = directory / "truth.vcf"
    with truth_vcf.open("w") as stream:
        names = "\t".join(sample["sample_id"] for sample in pilot["samples"])
        stream.write(f'##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{names}\n')
        for index, key in enumerate(found):
            genotypes = "\t".join(f"{x[0, index]}|{x[1, index]}" for x in truth)
            stream.write(f"{chrom}\t{key[1]}\t.\t{key[2]}\t{key[3]}"
                         f"\t.\tPASS\t.\tGT\t{genotypes}\n")
    checked(["bcftools", "view", "-Oz", "-o", directory / "truth.vcf.gz", truth_vcf],
            directory, "compress_truth")
    checked(["bcftools", "index", "-t", directory / "truth.vcf.gz"], directory, "index_truth")
    truth_vcf.unlink()
    original_ids = subprocess.check_output(["bcftools", "query", "-l", str(ref)], text=True).splitlines()
    excluded = set(pilot["excluded_reference_ids"])
    keep = [sample for sample in original_ids if sample not in excluded]
    if len(set(keep)) != len(keep) or set(keep) & set(pilot["donor_ids"]):
        raise ValueError("Imputation reference overlaps held-out founder identities")
    (directory / "reference.keep.txt").write_text("".join(sample + "\n" for sample in keep))
    disjoint = directory / "reference.disjoint.vcf.gz"
    checked(["bcftools", "view", "-S", directory / "reference.keep.txt", "-Oz",
             "--threads", "2", "-o", disjoint, ref], directory, "prepare_reference")
    actual_ids = subprocess.check_output(["bcftools", "query", "-l", str(disjoint)], text=True).splitlines()
    if actual_ids != keep or set(actual_ids) & excluded:
        raise ValueError("Held-out identities survived imputation reference extraction")
    array = np.asarray([key in observed_set for key in found])
    model_observed = np.asarray([key in observed_set for key in keys])
    audit = dict(status="running", chromosome=chrom, n_truth_sites=len(found),
                 n_array_sites=int(array.sum()), n_model_sites=len(keys),
                 n_model_observed_sites=int(model_observed.sum()),
                 n_real_observed_sites=len(observed),
                 n_array_sites_unavailable=len(observed_set - set(found)),
                 n_reference_before=len(original_ids), n_reference_after=len(actual_ids),
                 excluded_reference_present=sorted(set(original_ids) & excluded),
                 reference_disjoint=True, original_haplotypes_reproduced=True,
                 pilot_sha256=adapter.digest(output / "pilot.json"),
                 input_sha256={str(path): pilot["input_sha256"][str(path)] for path in inputs},
                 model_sha256=model_metadata["model_sha256"], samples=[])
    adapter.write_json(checkpoint, audit)
    predictions, scores = [], []
    model = bundle["model"]
    try:
        with threadpool_limits(limits=1):
            for sample_index, sample in enumerate(pilot["samples"]):
                sample_id = sample["sample_id"]
                work = directory / sample_id
                work.mkdir()
                target = work / "array.unphased.vcf"
                with target.open("w") as stream:
                    stream.write(f'##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n'
                                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                                 f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n')
                    for i in np.flatnonzero(array):
                        key = found[i]
                        genotype = sorted(truth[sample_index, :, i])
                        stream.write(f"{chrom}\t{key[1]}\t.\t{key[2]}\t{key[3]}"
                                     f"\t.\tPASS\t.\tGT\t{genotype[0]}/{genotype[1]}\n")
                prefix = work / "beagle"
                seconds = checked([
                    "java", "-Xmx16g", "-jar", pilot["beagle"], f"gt={target}",
                    f"ref={disjoint}", f"map={map_path}", f"out={prefix}",
                    "nthreads=6", f"seed={pilot['seed']}", "gp=true",
                ], work, "impute")
                imputed = Path(str(prefix) + ".vcf.gz")
                checked(["bcftools", "index", "-t", imputed], work, "index_imputed")
                frozen = work / "truth_panel.vcf.gz"
                checked(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Oz", "-o", frozen,
                         imputed, directory / "truth.vcf.gz"], work, "retain_truth_panel")
                checked(["bcftools", "index", "-t", frozen], work, "index_model_panel")
                recovered = adapter.load_haplotypes(frozen, [sample_id], chrom, keys)
                if np.any(recovered > 1):
                    raise ValueError("Beagle left missing frozen-model genotypes; no padding allowed")
                original = truth[sample_index][:, model_indices]
                stats = genotype_strata(original, recovered, model_observed)
                if stats["observed"]["n_genotype_errors"]:
                    raise ValueError("Imputation changed a supplied observed diploid genotype")
                inferred_array = adapter.load_haplotypes(frozen, [sample_id], chrom,
                                                         [found[i] for i in np.flatnonzero(array)])
                if not np.array_equal(inferred_array.sum(axis=0), truth[sample_index][:, array].sum(axis=0)):
                    raise ValueError("Imputation altered observed genotypes outside the model panel")
                # Repair only identifiable phase: correctly imputed true heterozygotes.
                repaired = recovered.copy()
                identifiable = (original.sum(axis=0) == 1) & (recovered.sum(axis=0) == 1)
                repaired[:, identifiable] = original[:, identifiable]
                if not np.array_equal(repaired.sum(axis=0), recovered.sum(axis=0)):
                    raise ValueError("Oracle phase control changed diploid genotype calls")
                for stage, haplotypes in (("original", original), ("array_imputed", recovered),
                                          ("oracle_phase_at_correct_hets", repaired)):
                    base = model.base.predict_proba(haplotypes)
                    probabilities = adapter.ordered_probabilities(model, base)
                    values, _ = adapter.aggregate_probabilities(probabilities, layout)
                    predictions.append(dict(sample_id=sample_id, chromosome=chrom, stage=stage,
                                            n_markers=len(keys),
                                            **dict(zip(adapter.SOURCES, values[0] / len(keys)))))
                for stratum, values in stats.items():
                    scores.append(dict(sample_id=sample_id, chromosome=chrom, stratum=stratum, **values))
                audit["samples"].append(dict(sample_id=sample_id, fit_seconds=seconds,
                                             genotype_errors=stats["all"]["n_genotype_errors"],
                                             switches=stats["all"]["switches"],
                                             switch_pairs=stats["all"]["switch_pairs"]))
                adapter.write_json(checkpoint, audit)
                # Keep exact truth/model/array genotypes; dense output is reproducible from logged inputs.
                # Delete only this sample's verified, successfully evaluated dense temporary output.
                imputed.unlink()
                Path(str(imputed) + ".tbi").unlink()
                print(f"chr{chrom} {sample_id}: imputed in {seconds:.1f}s, "
                      f"genotype_errors={stats['all']['n_genotype_errors']}, "
                      f"switches={stats['all']['switches']}/{stats['all']['switch_pairs']}", flush=True)
        pd.DataFrame(predictions).to_csv(directory / "predictions.tsv", sep="\t", index=False)
        pd.DataFrame(scores).to_csv(directory / "genotype_phase_metrics.tsv", sep="\t", index=False)
        verify_hashes(audit["input_sha256"])
        if adapter.digest(screen / f"candidates/{CANDIDATE}/models/chr{chrom}.pkl") != audit["model_sha256"]:
            raise ValueError("Frozen ancestry model changed")
        if adapter.digest(output / "pilot.json") != audit["pilot_sha256"]:
            raise ValueError("Pilot configuration changed during execution")
        # The disjoint reference can be regenerated from the pinned original plus audited keep list.
        disjoint.unlink()
        audit.update(status="complete", output_sha256=fingerprint_files(
            [p for p in directory.rglob("*") if p.is_file() and p != checkpoint]))
        adapter.write_json(checkpoint, audit)
    except Exception as error:
        audit.update(status="failed", error=str(error))
        adapter.write_json(checkpoint, audit)
        raise


def summarize(screen, output):
    rows, qc, audits = [], [], []
    for chrom in range(1, 23):
        directory = output / f"chr{chrom}"
        audit = json.loads((directory / "audit.json").read_text())
        if audit["status"] != "complete":
            raise ValueError(f"Chromosome {chrom} is not complete")
        verify_hashes(audit["output_sha256"])
        audits.append(audit)
        rows.append(pd.read_csv(directory / "predictions.tsv", sep="\t"))
        qc.append(pd.read_csv(directory / "genotype_phase_metrics.tsv", sep="\t"))
    predictions = pd.concat(rows, ignore_index=True)
    if predictions.duplicated(["sample_id", "chromosome", "stage"]).any():
        raise ValueError("Duplicate chromosome predictions")
    for _, group in predictions.groupby(["sample_id", "stage"]):
        if set(group.chromosome) != set(range(1, 23)):
            raise ValueError("Incomplete chromosome coverage")
    for source in adapter.SOURCES:
        predictions[source] *= predictions.n_markers
    total = predictions.groupby(["sample_id", "stage"])[[*adapter.SOURCES, "n_markers"]].sum()
    total[list(adapter.SOURCES)] = total[list(adapter.SOURCES)].div(total.n_markers, axis=0)
    total = total.reset_index()
    truth = np.array([.5, .5, 0, 0, 0, 0, 0])
    total["TV_distance"] = np.abs(total[list(adapter.SOURCES)] - truth).sum(axis=1) / 2
    saved = pd.read_csv(screen / f"candidates/{CANDIDATE}/predictions.tsv", sep="\t")
    before = total[total.stage.eq("original")].merge(saved, on="sample_id", validate="one_to_one")
    error = np.max(np.abs(before[list(adapter.SOURCES)].to_numpy()
                          - before[["pred_" + s for s in adapter.SOURCES]].to_numpy()))
    if len(before) != 5 or error > 1e-10:
        raise ValueError("Original full-genome control predictions did not reproduce the frozen screen")
    total.to_csv(output / "predictions.tsv", sep="\t", index=False)
    qc = pd.concat(qc, ignore_index=True)
    qc.to_csv(output / "genotype_phase_metrics.tsv", sep="\t", index=False)
    sums = qc.groupby(["sample_id", "stratum"]).sum(numeric_only=True).drop(columns="chromosome")
    sums["genotype_error_rate"] = sums.n_genotype_errors / sums.n_sites
    sums["switch_error_rate"] = sums.switches / sums.switch_pairs.where(sums.switch_pairs.gt(0))
    sums["phase_hamming_error_rate"] = sums.phase_hamming_errors / sums.n_correct_heterozygotes.where(
        sums.n_correct_heterozygotes.gt(0))
    sums.to_csv(output / "aggregate_genotype_phase_metrics.tsv", sep="\t")
    summary = dict(
        status="complete", n_controls=5, n_chromosomes=22, candidate=CANDIDATE,
        reference_identity_disjoint=True, original_genotypes_reproduced=True,
        maximum_original_prediction_error=float(error),
        stage_mean_tv=total.groupby("stage").TV_distance.mean().to_dict(),
        stage_mean_predictions=total.groupby("stage")[list(adapter.SOURCES)].mean().to_dict("index"),
        n_observed_input_sites=sum(a["n_array_sites"] for a in audits),
        n_real_observed_sites=sum(a["n_real_observed_sites"] for a in audits),
        oracle_phase_limitation="Only correctly imputed true heterozygotes have their original phase restored; genotype-error sites are not a complete oracle decomposition.",
        switch_metric_definition="Adjacent correctly imputed truth heterozygotes within each chromosome; excludes genotype-error sites, allows global haplotype swap.",
        no_production_changes=True, no_ancestry_training=True)
    summary["output_sha256"] = fingerprint_files(
        [output / name for name in ("predictions.tsv", "genotype_phase_metrics.tsv",
                                    "aggregate_genotype_phase_metrics.tsv")])
    adapter.write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "chromosome", "summarize"), required=True)
    parser.add_argument("--screen", type=Path, default=Path("/results"))
    parser.add_argument("--outdir", type=Path, default=Path("/pilot"))
    parser.add_argument("--reference", type=Path, default=Path("/imputer"))
    parser.add_argument("--beagle", type=Path, default=Path("/pilot/beagle.jar"))
    parser.add_argument("--chromosome", type=int, choices=range(1, 23))
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.screen, args.outdir, args.reference, args.beagle)
    elif args.stage == "chromosome":
        if args.chromosome is None:
            parser.error("--chromosome is required for chromosome execution")
        chromosome(args.screen, args.outdir, args.chromosome)
    else:
        summarize(args.screen, args.outdir)


if __name__ == "__main__":
    main()
