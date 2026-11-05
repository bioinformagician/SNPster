#!/usr/bin/env python3
import itertools
import json
from pathlib import Path


def n_product(grid):
    n = 1
    for values in grid.values():
        n *= len(values)
    return n

cfg_path = Path(__import__('os').environ.get('SEARCH_SPACE', 'config/search_space.json'))
search = json.loads(cfg_path.read_text())

panel_n = n_product(search['panel_search'])
print(f"Search space: {cfg_path}")
print(f"Panel-probability configurations : {panel_n}")
print()

total_models = 0
for strategy in search['flare_model_strategies']:
    n = panel_n * n_product(strategy.get('grid', {}))
    total_models += n
    print(f"{strategy['name']:28s} {n:8d}  ({strategy['family']})")
print('-' * 72)
print(f"Trainable FLARE models          : {total_models}")
print()

st2 = search['stage2']
decoders = 0
for d in st2['decoders']:
    decoders += n_product(d.get('grid', {}))
calibrators = sum(n_product(c.get('grid', {})) for c in st2['calibrators'])
weights = len(st2['training_weights'])
per_rep = decoders * calibrators * weights
keep = int(search['stage1']['retain_top_representations'])
forced = len(st2.get('fixed_candidates', []))

print(f"Decoder variants / representation: {decoders}")
print(f"Calibrator variants              : {calibrators}")
print(f"Training-weight variants         : {weights}")
print(f"Stage2 candidates / representation: {per_rep}")
print(f"Stage1 retained representations  : {keep} (+ up to {forced} forced fixed)")
print(f"Max default stage2 candidates    : {per_rep * (keep + forced)}")
print(f"Stage2 retained candidates       : {st2['retain_top_candidates']}")
print(f"Finalists                        : {search['stage3']['retain_finalists']}")
