#!/usr/bin/env python3
"""Intersect FLARE held-out additive-source individuals with actual high-density VCF IDs.

Never silently substitute training/reference-panel donors for missing validation donors.
Output is compatible with make_plans.py (sample, fine, source_group).
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import subprocess
from collections import Counter
from pathlib import Path

SOURCES = ('EUR','EAS','SAS','SEA','GME','NAM','OCE','SSA')


def sample_ids(vcf: str, bcftools='bcftools') -> set[str]:
    res=subprocess.run([bcftools,'query','-l',str(vcf)],check=True,capture_output=True,text=True)
    arr=[x.strip() for x in res.stdout.splitlines() if x.strip()]
    if len(arr)!=len(set(arr)) or not arr:
        raise ValueError('High-density VCF has empty or duplicate sample IDs')
    return set(arr)


def parse_tsv(path):
    with Path(path).open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'sample','source_group'}.issubset(reader.fieldnames or []):
            raise ValueError(f'Expected sample and source_group columns in {path}')
        result={}
        for r in reader:
            sample=r['sample'].strip()
            if not sample: continue
            group=r['source_group'].strip().upper()
            if sample in result:
                raise ValueError(f'Duplicate held-out sample {sample} in {path}')
            result[sample]=(group,str(r.get('fine','')).strip())
    return result


def parse_ids(path):
    with Path(path).open(newline='',encoding='utf-8') as f:
        first=f.readline().strip().split('\t');f.seek(0)
        if 'sample' in first:
            return {x['sample'].strip() for x in csv.DictReader(f,delimiter='\t') if x.get('sample','').strip()}
        return {l.split()[0] for l in f if l.strip() and not l.startswith('#')}


def load_label_metadata(path):
    labels = {}
    populations = {}
    with Path(path).open(newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not {'sample', 'fine', 'source_group'}.issubset(reader.fieldnames or []):
            raise ValueError(f'Label metadata needs sample, fine and source_group: {path}')
        for row in reader:
            sample = row['sample'].strip()
            label = (row['source_group'].strip().upper(), row['fine'].strip())
            if not sample or not label[1]:
                raise ValueError(f'Empty sample or fine-population label in {path}')
            if sample in labels and labels[sample] != label:
                raise ValueError(f'Conflicting authoritative labels for {sample}')
            labels[sample] = label
            if label[0]:
                if label[1] in populations and populations[label[1]] != label[0]:
                    raise ValueError(f'Conflicting source labels for population {label[1]}')
                populations[label[1]] = label[0]
    if not labels:
        raise ValueError('Empty authoritative label metadata')
    return labels


def validate_labels(rows, labels, description):
    for sample, label in rows.items():
        if sample not in labels:
            raise ValueError(f'{description}: no authoritative label for {sample}')
        if labels[sample] != label:
            raise ValueError(f'{description}: label mismatch for {sample}: '
                             f'provided={label}, authoritative={labels[sample]}')


def mosaic_donor_ids(path):
    donors = set()
    with Path(path).open(newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'donor_weights_json' not in (reader.fieldnames or []):
            raise ValueError(f'Mosaic metadata needs donor_weights_json: {path}')
        for line, row in enumerate(reader, 2):
            weights = json.loads(row['donor_weights_json'])
            if not isinstance(weights, dict) or not weights:
                raise ValueError(f'{path}:{line}: empty or invalid mosaic donors')
            values = [float(v) for v in weights.values()]
            if (any(not isinstance(s, str) or not s.strip() for s in weights)
                    or any(not math.isfinite(v) or v < 0 for v in values)
                    or abs(sum(values) - 1) > 1e-5):
                raise ValueError(f'{path}:{line}: invalid mosaic donor weights')
            donors.update(weights)
    if not donors:
        raise ValueError('Empty mosaic donor exclusion set')
    return donors


def read_reference_panel(path):
    labels = {}
    with Path(path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) != 2 or fields[0] in labels:
                raise ValueError(f'Invalid or duplicate reference-panel row: {line.strip()}')
            labels[fields[0]] = fields[1]
    return labels


def check_existing_output(outdir, design):
    outdir = Path(outdir)
    if not (outdir / 'pedigree').exists():
        return
    existing = outdir / 'pedigree/inputs/eligible_heldout_founders.tsv.summary.json'
    if not existing.exists():
        artifacts = ('evaluations', 'truth', 'plans', 'panel_genomes', 'chromosomes')
        if any((outdir / 'pedigree' / name).is_dir()
               and any((outdir / 'pedigree' / name).iterdir()) for name in artifacts):
            raise ValueError('Existing pedigree output has no founder audit; use a new output directory')
        return
    audit = json.loads(existing.read_text())
    if (audit.get('strict_reference_disjoint') is not True or audit.get('labels_validated') is not True
            or audit.get('mosaic_training_exclusion_checked') is not True
            or audit.get('flare_reference_overlap_count') != 0):
        raise ValueError('Historical pedigree output is not validated/reference-disjoint; '
                         'use a new output directory rather than overwriting it')
    manifest = outdir / 'pedigree/plans/plans/manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()).get('design') != design:
        raise ValueError('Existing pedigree design differs; use a separately labeled output directory')


def prepare(validation,training,vcf_samples,out,min_per_source=2,required_sources=SOURCES,
            reference_samples=None, strict_reference_disjoint=True, additional_metadata=None,
            label_metadata=None, reference_labels=None, mosaic_metadata=None):
    if reference_samples is None:
        raise ValueError('Reference sample IDs are required for disjoint founder selection')
    if not strict_reference_disjoint:
        raise ValueError('Reference-disjoint founder selection cannot be disabled')
    if label_metadata is None:
        raise ValueError('Authoritative label metadata is required')
    labels = load_label_metadata(label_metadata)
    held=parse_tsv(validation)
    validate_labels(held, labels, 'Validation metadata')
    validate_labels(parse_tsv(training), labels, 'Training metadata')
    if reference_labels is not None:
        if set(reference_labels) != set(reference_samples):
            raise ValueError('Reference VCF and reference-panel sample IDs differ')
        for sample, fine in reference_labels.items():
            if sample not in labels or labels[sample][1] != fine:
                raise ValueError(f'Reference-panel fine label mismatch for {sample}: {fine}')
    excluded=parse_ids(training)
    if not excluded:
        raise ValueError('Training sample exclusion set is empty; refusing potential data leakage')
    mosaic_donors = mosaic_donor_ids(mosaic_metadata) if mosaic_metadata else set()
    excluded |= mosaic_donors
    additional = parse_tsv(additional_metadata) if additional_metadata else {}
    validate_labels(additional, labels, 'Additional founder metadata')
    for sample, label in additional.items():
        if sample in excluded:
            raise ValueError(f'Additional founder {sample} is a training donor')
        if sample not in vcf_samples:
            raise ValueError(f'Additional founder {sample} is absent from the high-density VCF')
        if label[0] not in SOURCES:
            raise ValueError(f'Additional founder {sample} has unsupported source {label[0]}')
        if sample in held and held[sample] != label:
            raise ValueError(f'Conflicting metadata for additional founder {sample}')
        held[sample] = label
    matches=[]
    overlaps = []
    for s, (source, fine) in sorted(held.items()):
        if s in vcf_samples and s not in excluded and source in SOURCES:
            if s in reference_samples:
                overlaps.append(s)
                continue
            matches.append((s,fine,source))
    counts=Counter(x[2] for x in matches)
    missing=[s for s in required_sources if counts[s]<min_per_source]
    if missing:
        raise ValueError('Insufficient *held-out high-density* founders in '+str(missing)+'. Counts: '+
                         str(dict(sorted(counts.items())))+'. Add external donor-disjoint reference genomes; '
                         'do not silently use FLARE training donors.')
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    # Preserve mtime for deterministic -resume cache hits.
    import io
    buf=io.StringIO()
    w=csv.writer(buf,delimiter='\t',lineterminator='\n')
    w.writerow(['sample','fine','source_group']);w.writerows(matches)
    new_text=buf.getvalue()
    if not Path(out).exists() or Path(out).read_text()!=new_text:
        Path(out).write_text(new_text)
    summary={'flare_reference_overlap_count':0,'flare_reference_overlap_examples':[],
             'excluded_reference_overlap_ids':overlaps,
             'excluded_reference_overlap_count':len(overlaps),
             'additional_founder_ids':sorted(additional),
             'additional_reference_overlap_ids':sorted(set(additional) & set(reference_samples or [])),
             'strict_reference_disjoint':True, 'labels_validated':True,
             'mosaic_training_exclusion_checked':mosaic_metadata is not None,
             'excluded_mosaic_donor_count':len(mosaic_donors),
             'available_highdensity_ids':len(vcf_samples),'heldout_ids':len(held),
             'excluded_training_ids':len(excluded),'selected_donors':len(matches),
             'selected_per_source':{s:counts[s] for s in SOURCES},'min_per_source':min_per_source}
    summ_path=Path(str(out)+'.summary.json')
    summary_content=json.dumps(summary,indent=2)+'\n'
    if not summ_path.exists() or summ_path.read_text()!=summary_content:
        summ_path.write_text(summary_content)
    return summary


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--validation',required=True);p.add_argument('--training',required=True)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--vcf');group.add_argument('--sample-list')
    p.add_argument('--output',required=True)
    p.add_argument('--min-per-source',type=int,default=2)
    p.add_argument('--reference-vcf',required=True)
    p.add_argument('--reference-panel',required=True)
    p.add_argument('--label-metadata',required=True)
    p.add_argument('--mosaic-metadata',required=True)
    p.add_argument('--benchmark-outdir',required=True)
    p.add_argument('--design',choices=('full','reduced'),required=True)
    p.add_argument('--strict-reference-disjoint',action='store_true')
    p.add_argument('--additional-metadata', help='Explicitly approved supplementary donor TSV; reference overlap is audited')
    a=p.parse_args()
    if a.min_per_source<1:p.error('min-per-source must be positive')
    check_existing_output(a.benchmark_outdir,a.design)
    samples=sample_ids(a.vcf) if a.vcf else set(Path(a.sample_list).read_text().split())
    refs=sample_ids(a.reference_vcf)
    print(json.dumps(prepare(a.validation,a.training,samples,a.output,a.min_per_source,
                             reference_samples=refs,strict_reference_disjoint=True,
                             additional_metadata=a.additional_metadata,
                             label_metadata=a.label_metadata,
                             reference_labels=read_reference_panel(a.reference_panel),
                             mosaic_metadata=a.mosaic_metadata),indent=2))
if __name__=='__main__':main()
