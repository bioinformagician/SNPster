#!/usr/bin/env python3
"""Intersect FLARE held-out additive-source individuals with actual high-density VCF IDs.

Never silently substitute training/reference-panel donors for missing validation donors.
Output is compatible with make_plans.py (sample, fine, source_group).
"""
from __future__ import annotations
import argparse
import csv
import json
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


def prepare(validation,training,vcf_samples,out,min_per_source=2,required_sources=SOURCES, reference_samples=None, strict_reference_disjoint=False):
    held=parse_tsv(validation)
    excluded=parse_ids(training)
    if not excluded:
        raise ValueError('Training sample exclusion set is empty; refusing potential data leakage')
    matches=[]
    for s, (source, fine) in sorted(held.items()):
        if s in vcf_samples and s not in excluded and source in SOURCES:
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
    overlaps = sorted({r[0] for r in matches} & set(reference_samples or []))
    if overlaps and strict_reference_disjoint:
        raise ValueError(f'{len(overlaps)} candidate founders also appear in the FLARE reference; strict reference exclusion enabled')
    summary={'flare_reference_overlap_count':len(overlaps),'flare_reference_overlap_examples':overlaps[:10],
             'strict_reference_disjoint':strict_reference_disjoint,
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
    p.add_argument('--reference-vcf')
    p.add_argument('--strict-reference-disjoint',action='store_true')
    a=p.parse_args()
    if a.min_per_source<1:p.error('min-per-source must be positive')
    samples=sample_ids(a.vcf) if a.vcf else set(Path(a.sample_list).read_text().split())
    refs=sample_ids(a.reference_vcf) if a.reference_vcf else None
    print(json.dumps(prepare(a.validation,a.training,samples,a.output,a.min_per_source,
                             reference_samples=refs,strict_reference_disjoint=a.strict_reference_disjoint),indent=2))
if __name__=='__main__':main()
