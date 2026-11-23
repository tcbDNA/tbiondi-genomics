#!/usr/bin/env python3
"""Check gene termini: start codon at alignment start, distance to first in-frame stop after
alignment end (yeast mito code: TAA, TAG), and whether an upstream ATG exists.
Usage: termini.py genome_circ.fna refs.faa miniprot.paf"""
import re, sys
def readfa(p):
    d = {}; n = None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l)
    return {k: ''.join(v).upper() for k, v in d.items()}
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
G = readfa(sys.argv[1]); R = readfa(sys.argv[2])
for line in open(sys.argv[3]):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[4] not in '+-' or '_v2' in f[0] or '_v3' in f[0]: continue
    q, qlen, qs, qe, strand, t, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])
    g = G[t]
    if strand == '+': cod = g[ts:te]; down = g[te:te+300]; up = g[max(0,ts-60):ts]
    else: cod = rc(g[ts:te]); down = rc(g[max(0,te-0-0):te][:0]) or rc(g[max(0,ts-300):ts]); up = rc(g[te:te+60])
    stop_d = None
    for i in range(0, len(down)-2, 3):
        if down[i:i+3] in ('TAA', 'TAG'): stop_d = i // 3; break
    print('%s\t%s\tq=%d-%d/%d\tstart_codon=%s\tlast_codon=%s\tcodons_to_stop_after_end=%s\tstop=%s' % (
        t, q, qs, qe, qlen, cod[:3], cod[-3:], stop_d, down[3*stop_d:3*stop_d+3] if stop_d is not None else '-'))
