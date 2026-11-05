#!/usr/bin/env python3
"""Convert an existing GRCh38 4-column PLINK map into 22 Ped-sim maps.
Only provides SEX-AVERAGED map positions; original build must already be GRCh38.
"""
import argparse
from pathlib import Path
from collections import defaultdict

def convert(infile, outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    rows=defaultdict(list)
    for ln in Path(infile).open():
        s=ln.split()
        if not s or s[0].startswith('#'): continue
        if len(s)!=4: raise ValueError('Expected PLINK: CHROM MARKER cM POS (four fields)')
        chrom=s[0].removeprefix('chr')
        if not chrom.isdigit() or not 1<=int(chrom)<=22: continue
        cM=float(s[2]);pos=int(s[3]);rows[chrom].append((pos,cM))
    for chrom,data in rows.items():
        data=sorted(set(data))
        if len(data)<2: raise ValueError(f'Fewer than two map points on chromosome {chrom}')
        last_pos=0;last_cm=-float('inf')
        lines=['#chrom\tposition\tsex_averaged_cM\n']
        for pos,cm in data:
            if pos<=last_pos or cm<last_cm: raise ValueError(f'Non-monotonic map chr{chrom}:{pos}')
            last_pos=pos;last_cm=cm
            lines.append(f'{chrom}\t{pos}\t{cm:.8f}\n')
        target=outdir/f'chr{chrom}.simmap'
        content=''.join(lines)
        if not target.exists() or target.read_text()!=content:
            target.write_text(content)
    print(f'Generated {len(rows)} sex-averaged Ped-sim maps; these are not sex-specific maps')

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--plink-map',required=True)
    a.add_argument('--outdir',type=Path,required=True)
    p=a.parse_args();convert(p.plink_map,p.outdir)
if __name__=='__main__':main()
