#!/usr/bin/env python3
"""Ensure Ped-sim actually used precisely the specified labelled reference donors."""
import argparse

def read_map(path):
    out={}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'): continue
            parts=line.split()
            if len(parts)!=2: raise ValueError(f'Expected two columns in {path}: {line}')
            if parts[0] in out: raise ValueError(f'Duplicate founder {parts[0]} in {path}')
            out[parts[0]]=parts[1]
    return out

def run(expected,actual):
    a,b=read_map(expected),read_map(actual)
    if a!=b:
        missing={k:v for k,v in a.items() if k not in b}
        extra={k:v for k,v in b.items() if k not in a}
        wrong={k:(v,b[k]) for k,v in a.items() if k in b and b[k]!=v}
        raise ValueError(f'Ped-sim founder mismatch: missing={missing}; extra={extra}; wrong={wrong}')
    print(f'Founder assignment verified: {len(a)} unique, fully specified Ped-sim founders')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--expected',required=True);p.add_argument('--actual',required=True)
    a=p.parse_args();run(a.expected,a.actual)
