#!/usr/bin/env python3
"""Build unrelated, ancestry-labelled founder trees in Ped-sim .def syntax.

All branches have exactly one simulated individual (i1); every parental
relationship is explicit. New pure-source founders are born in the generation
just before their child, so deep backcrosses do not mate across generations.
"""
import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

SOURCES = ('EUR','EAS','SAS','SEA','GME','NAM','OCE','SSA')

@dataclass
class Node:
    source: str | None = None
    left: 'Node | None' = None
    right: 'Node | None' = None
    generation: int = 0
    branch: int = 0
    sex: str = 'M'
    person_id: str = ''
    donor: str = ''

    @property
    def founder(self): return self.source is not None

    def leaves(self):
        if self.founder: return [self]
        return self.left.leaves() + self.right.leaves()

    def expected(self):
        if self.founder: return {self.source: 1.0}
        a, b = self.left.expected(), self.right.expected()
        return {s: (a.get(s,0.0)+b.get(s,0.0))/2 for s in SOURCES if a.get(s,0.0)+b.get(s,0.0)}


def leaf(s): return Node(source=s)
def cross(a,b): return Node(left=a, right=b)


def candidates(rng):
    """Cycle through benchmarks, including all distinct source pairs and uneven trees."""
    pair_list = [(a,b) for i,a in enumerate(SOURCES) for b in SOURCES[i+1:]]
    pair_list.remove(('EUR','SEA'))
    rng.shuffle(pair_list)
    pair_list.insert(0, ('EUR','SEA'))
    counter = 0
    while True:
        kind = counter % 9
        a,b = pair_list[(counter//9) % len(pair_list)]
        if kind == 0:
            # F1 (ensure direct EUR/Filipino example at beginning)
            if counter == 0: a,b = 'EUR','SEA'
            tree = cross(leaf(a),leaf(b)); cls='f1'
        elif kind == 1:
            # Five independent mating events => 1/32 from original founder.
            if counter == 1: a,b = 'SEA','EUR'
            tree = cross(leaf(a), leaf(b))
            for _ in range(4): tree=cross(tree,leaf(b))
            cls='dilution_3pct'
        elif kind == 2:
            order = ['EUR','SSA','EAS','SEA'] if counter == 2 else rng.sample(list(SOURCES),4)
            tree = cross(cross(leaf(order[0]),leaf(order[1])),cross(leaf(order[2]),leaf(order[3])))
            cls='four_way_balanced'
        elif kind == 3:
            dominant, minority = ('SSA','EUR') if counter==3 else (a,b)
            tree = cross(cross(leaf(minority),leaf(dominant)),cross(leaf(dominant),leaf(dominant)))
            cls='four_way_repeated_founders'
        elif kind == 4:
            order = rng.sample(list(SOURCES),4)
            tree=cross(cross(leaf(order[0]),leaf(order[1])),cross(leaf(order[2]),leaf(order[3])))
            cls='admixed_x_admixed'
        elif kind == 5:
            order = ['EAS','SEA','EUR','SSA'] if counter==5 else rng.sample(list(SOURCES),4)
            # 75% last group, 12.5% third, 6.25% first two
            tree=cross(cross(cross(leaf(order[0]),leaf(order[1])),leaf(order[2])),leaf(order[3]))
            tree=cross(tree,leaf(order[3])); cls='four_way_uneven'
        elif kind == 6:
            order=rng.sample(list(SOURCES),8)
            ab=[cross(leaf(order[i]),leaf(order[i+1])) for i in range(0,8,2)]
            tree=cross(cross(ab[0],ab[1]),cross(ab[2],ab[3]))
            cls='eight_way_admixed_x_admixed'
        elif kind == 7:
            order = rng.sample(list(SOURCES),4)
            tree=cross(cross(leaf(order[0]),leaf(order[1])),cross(leaf(order[2]),leaf(order[3])))
            tree=cross(tree,leaf(order[0]))
            tree=cross(tree,leaf(order[0]))
            tree=cross(tree,leaf(order[0]))
            cls='rare_3pct_multisource'
        else:
            # Serial mating entirely with new F1 admixed partners (not pure donors).
            # Four-way descendant starts at 25% per original source. Three
            # admixed-parent matings dilute original minority lines to 3.125%.
            o=rng.sample(list(SOURCES),4)
            tree=cross(cross(leaf(o[0]),leaf(o[1])),cross(leaf(o[2]),leaf(o[3])))
            protected=o[0]  # never reintroduce this original source during dilution
            for _ in range(3):
                p,q=rng.sample([s for s in SOURCES if s!=protected],2)
                tree=cross(tree,cross(leaf(p),leaf(q)))
            cls='serial_admixed_mating_3pct'
        counter+=1
        yield cls,tree


def load_ids(path):
    """Accept source_train_metadata TSV or plain sample-id files."""
    if not path: return set()
    with open(path, encoding='utf-8') as fh:
        first = fh.readline().strip().split('\t')
        fh.seek(0)
        if 'sample' in first:
            return {r['sample'].strip() for r in csv.DictReader(fh, delimiter='\t') if r.get('sample','').strip()}
        return {l.split()[0] for l in fh if l.strip() and not l.startswith('#')}


def load_pool(path, exclusions):
    pools=defaultdict(list)
    observed_samples=set()
    with open(path, newline='', encoding='utf-8') as fh:
        reader=csv.DictReader(fh, delimiter='\t')
        if not {'sample','source_group'}.issubset(reader.fieldnames or []):
            raise ValueError('Donor TSV needs columns sample and source_group; input must be source-labelled')
        for r in reader:
            sample=r['sample'].strip()
            group=r['source_group'].strip().upper()
            if group not in SOURCES or sample in exclusions or not sample: continue
            if sample in observed_samples: raise ValueError(f'Duplicate eligible sample across groups: {sample}')
            observed_samples.add(sample)
            pools[group].append({'sample':sample,'fine':str(r.get('fine') or '').strip()})
    for group, donors in pools.items():
        pools[group]=sorted(donors,key=lambda r:r['sample'])
    missing=[s for s in SOURCES if not pools[s]]
    if missing: raise ValueError(f'No eligible founder donors for {missing}; supply reference samples from every additive group')
    return pools


def max_depth(n):
    return 1 if n.founder else 1+max(max_depth(n.left),max_depth(n.right))


def assign_tree(root,name):
    """Assign all edges to adjacent generations and every node a branch."""
    root.generation=max_depth(root)
    generations=defaultdict(list)
    def visit(n):
        generations[n.generation].append(n)
        if not n.founder:
            n.left.generation=n.right.generation=n.generation-1
            n.left.sex='M'; n.right.sex='F'
            visit(n.left); visit(n.right)
    visit(root)
    for gen in sorted(generations):
        for i,n in enumerate(generations[gen],1):
            n.branch=i
            n.person_id=f'{name}1_g{gen}-b{i}-i1'
    return generations


def write_family(name, root, cls, rng, pools, allow_reuse_global, used_global):
    generations=assign_tree(root,name)
    used_here=set()
    fine_used_here=defaultdict(set)
    founders={}
    for n in root.leaves():
        options=[x for x in pools[n.source] if x['sample'] not in used_here and
                 (allow_reuse_global or x['sample'] not in used_global)]
        if not options: raise ValueError(f'Insufficient {n.source} held-out donors for pedigree {name}. '
                                         'Use smaller batches or --allow-reuse-across-families for development only.')
        newfine=[x for x in options if x['fine'] and x['fine'] not in fine_used_here[n.source]]
        chosen=rng.choice(newfine or options)
        donor=chosen['sample']
        fine_used_here[n.source].add(chosen['fine'])
        used_here.add(donor); used_global.add(donor)
        n.donor=donor
        founders[n.person_id]={'vcf_sample': donor,'source':n.source,'fine':chosen['fine']}
    lines=[f'def {name} 1 {root.generation}']
    metadata=[]
    for gen in sorted(generations):
        nodes=generations[gen]
        specs=[]
        for n in nodes:
            specs.append(f'{n.branch}s{n.sex}')
            if gen>1:
                if n.founder: specs.append(f'{n.branch}:')
                else: specs.append(f'{n.branch}:{n.left.branch}_{n.right.branch}')
            metadata.append({
                'sample_id': n.person_id, 'family_id': name, 'class': cls,
                'generation':gen, 'father_id':n.left.person_id if not n.founder else '',
                'mother_id':n.right.person_id if not n.founder else '',
                'is_founder':n.founder, 'founder_donor':n.donor,
                'founder_source':n.source or '', 'expected': n.expected(),
                'number_of_sources':len(n.expected()),'number_of_unique_founders':len(n.leaves()),
            })
        lines.append(f"{gen} 1 {len(nodes)} " + ' '.join(specs))
    return '\n'.join(lines)+'\n', founders, metadata


def generate(metadata, excludes, outdir, families, batch_size, seed, reuse):
    rng=random.Random(seed)
    exclusion=set()
    for file in excludes: exclusion.update(load_ids(file))
    pool=load_pool(metadata,exclusion)
    outdir.mkdir(parents=True,exist_ok=True)
    used=set()
    planned=[]; by_batch=defaultdict(list)
    for i, (cls,root) in zip(range(families),candidates(rng)):
        name=f'fam{i+1:06d}'
        definition, founders, rows=write_family(name,root,cls,rng,pool,reuse,used)
        batch_id=f'batch_{i//batch_size+1:05d}'
        by_batch[batch_id].append((definition,founders,rows))
        planned.extend(rows)
    for batch, members in by_batch.items():
        directory=outdir/batch;directory.mkdir()
        founders_all={};metadata_all=[]
        with (directory/'pedigrees.def').open('w') as f:
            for definition,founders,rows in members:
                f.write(definition+'\n');founders_all.update(founders);metadata_all.extend(rows)
        with (directory/'set_founders.tsv').open('w') as f:
            for fid,v in sorted(founders_all.items()): f.write(f"{fid}\t{v['vcf_sample']}\n")
        (directory/'founder_samples.txt').write_text(''.join(f'{x}\n' for x in sorted({v['vcf_sample'] for v in founders_all.values()})))
        (directory/'descendant_samples.txt').write_text(''.join(f"{r['sample_id']}\n" for r in metadata_all if not r['is_founder']))
        (directory/'pedigrees.json').write_text(json.dumps({'batch':batch,'founders':founders_all,'nodes':metadata_all},indent=2))
    (outdir/'manifest.json').write_text(json.dumps({'seed':seed,'n_families':families,
        'n_batches':len(by_batch),'batch_size':batch_size,'eligible_donors':{x:len(pool[x]) for x in SOURCES},
        'excluded_ids':len(exclusion),'reuse_across_families':reuse,
        'classes':dict(Counter(r['class'] for r in planned if r['generation']==max(t['generation'] for t in planned if t['family_id']==r['family_id'])))},indent=2))
    with (outdir/'expected_ancestry.tsv').open('w') as f:
        w=csv.DictWriter(f, delimiter='\t',extrasaction='ignore',fieldnames=['sample_id','family_id','class','generation','father_id','mother_id','is_founder','founder_donor','founder_source']+[f'expected_{s}' for s in SOURCES])
        w.writeheader()
        for r in planned:
            w.writerow({**r,**{f'expected_{s}':r['expected'].get(s,0.) for s in SOURCES}})


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--exclude',type=Path,action='append',default=[])
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--families',type=int,default=20)
    p.add_argument('--batch-size',type=int,default=5)
    p.add_argument('--seed',type=int,default=20261008)
    p.add_argument('--allow-reuse-across-families',action='store_true')
    a=p.parse_args()
    if a.families<1 or a.batch_size<1: p.error('families and batch-size must be positive')
    generate(a.metadata,a.exclude,a.outdir,a.families,a.batch_size,a.seed,a.allow_reuse_across_families)
if __name__=='__main__':main()
