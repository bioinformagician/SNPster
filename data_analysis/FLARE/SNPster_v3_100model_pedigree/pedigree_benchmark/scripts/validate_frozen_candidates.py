#!/usr/bin/env python3
"""Validate saved candidate pairs and checksum-bound static pedigree inputs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import joblib

SOURCES = ('EUR', 'EAS', 'SAS', 'SEA', 'GME', 'NAM', 'OCE', 'SSA')
FIELDS = ('candidate_id', 'representation', 'candidate_dir', 'model_path', 'model_config')


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify_reusable_inputs(directory, genetic_map, reference_vcf):
    directory = Path(directory)
    manifest = json.loads((directory / 'manifest.json').read_text())
    for key, source in (('genetic_map_sha256', genetic_map), ('reference_vcf_sha256', reference_vcf)):
        if manifest.get(key) != sha256_file(source):
            raise ValueError(f'Reusable pedigree input source changed: {source}')
    expected = {'panel_sites.tsv'} | {f'maps/chr{chrom}.simmap' for chrom in range(1, 23)}
    if set(manifest.get('files', {})) != expected:
        raise ValueError('Reusable pedigree inputs must contain panel sites and all 22 chromosome maps')
    for name, expected_digest in manifest['files'].items():
        if sha256_file(directory / name) != expected_digest:
            raise ValueError(f'Reusable pedigree artifact changed: {name}')
    return manifest


def validate(manifest_path, output, report_path, reuse_dir=None, genetic_map=None, reference_vcf=None):
    with Path(manifest_path).open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not set(FIELDS).issubset(reader.fieldnames or []):
            raise ValueError(f'Frozen manifest needs columns {FIELDS}')
        rows = list(reader)
    if not rows or len({r['candidate_id'] for r in rows}) != len(rows):
        raise ValueError('Frozen manifest is empty or contains duplicate candidate IDs')
    normalized = []
    hashes = {}
    for row in rows:
        cid, rep = row['candidate_id'], row['representation']
        directory = Path(row['candidate_dir']).resolve(strict=True)
        model_path = Path(row['model_path']).resolve(strict=True)
        config_path = Path(row['model_config']).resolve(strict=True)
        config = json.loads(config_path.read_text())
        candidate = json.loads((directory / 'candidate.json').read_text())
        metrics = json.loads((directory / 'metrics.json').read_text())
        if (candidate.get('candidate_id') != cid or candidate.get('representation') != rep
                or metrics.get('candidate_id') != cid or metrics.get('representation') != rep
                or config.get('rep_id') != rep or directory.name != f'candidate_{cid}'):
            raise ValueError(f'Frozen artifact identity mismatch: {cid}')
        for key in ('min_maf', 'min_mac', 'flare_seed'):
            if key not in config:
                raise ValueError(f'{config_path}: missing inference setting {key}')
        saved = joblib.load(directory / 'model.joblib')
        if saved.get('candidate') != candidate or tuple(saved.get('sources', ())) != SOURCES:
            raise ValueError(f'Saved decoder identity/source ordering mismatch: {cid}')
        columns = saved.get('feature_columns')
        if not columns or len(columns) != len(set(columns)) or columns != metrics.get('feature_columns'):
            raise ValueError(f'Invalid saved feature ordering: {cid}')
        with model_path.open() as fh:
            ancestries = next(line.split() for line in fh if line.strip() and not line.startswith('#'))
        if ancestries != columns:
            raise ValueError(f'FLARE ancestry columns differ from decoder features: {cid}')
        with (directory / 'predictions.tsv').open(newline='') as fh:
            predictions = csv.DictReader(fh, delimiter='\t')
            required = {'sample_id', 'role'} | set(SOURCES) | {f'truth_{s}' for s in SOURCES}
            if not required.issubset(predictions.fieldnames or []):
                raise ValueError(f'Missing reference prediction columns: {cid}')
            if not any(r['role'] == 'validation' for r in predictions):
                raise ValueError(f'Missing saved reference validation predictions: {cid}')
        normalized.append({
            **row, 'candidate_dir': str(directory), 'model_path': str(model_path),
            'model_config': str(config_path), **{k: config[k] for k in ('min_maf', 'min_mac', 'flare_seed')},
        })
        hashes[cid] = {name: sha256_file(directory / name) for name in
                       ('candidate.json', 'metrics.json', 'model.joblib', 'predictions.tsv')}
        hashes[cid]['flare_model'] = sha256_file(model_path)
        hashes[cid]['model_config'] = sha256_file(config_path)
    reuse = verify_reusable_inputs(reuse_dir, genetic_map, reference_vcf) if reuse_dir else None
    with Path(output).open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(FIELDS) + ['min_maf', 'min_mac', 'flare_seed'], delimiter='\t')
        writer.writeheader()
        writer.writerows({k: row[k] for k in writer.fieldnames} for row in normalized)
    report = {'candidate_count': len(rows), 'representation_count': len({r['representation'] for r in rows}),
              'frozen_artifact_sha256': hashes, 'reused_static_inputs': reuse,
              'offspring_reused': False, 'models_retrained': False}
    Path(report_path).write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--report', required=True)
    parser.add_argument('--reusable-inputs')
    parser.add_argument('--genetic-map', required=True)
    parser.add_argument('--reference-vcf', required=True)
    args = parser.parse_args()
    report = validate(args.manifest, args.output, args.report, args.reusable_inputs,
                      args.genetic_map, args.reference_vcf)
    print(f"Validated {report['candidate_count']} frozen pairs across "
          f"{report['representation_count']} representations")


if __name__ == '__main__':
    main()
