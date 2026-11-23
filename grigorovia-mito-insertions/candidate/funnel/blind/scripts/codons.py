#!/usr/bin/env python3
"""Collect in-frame codons of each gene from the best miniprot alignment per gene
(M and D ops only; frameshift nt, introns and split codons skipped).
Usage: codons.py genome_circ.fna genome_len refs.faa miniprot.paf
Prints codon counts for selected codons and lists every in-frame occurrence of
motif codons with their coordinates."""
import re, sys
from collections import Counter, defaultdict
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
gfa, glen, rfa, paf = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
genome = readfa(gfa)
best = {}
for line in open(paf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[4] not in '+-': continue
    if '_v2' in f[0] or '_v3' in f[0]: continue
    gene = f[0].split('_', 1)[1].replace('_v1', '')
    ms = int(next(x[5:] for x in f[12:] if x.startswith('ms:i:')))
    if gene not in best or ms > best[gene][0]: best[gene] = (ms, f)
allc = Counter(); occ = defaultdict(list); fs_next = []
for gene, (ms, f) in sorted(best.items()):
    q, qs, strand, t, ts, te = f[0], int(f[2]), f[4], f[5], int(f[7]), int(f[8])
    cg = next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
    cod = genome[t][ts:te] if strand == '+' else rc(genome[t][ts:te])
    def fwd(c):
        x = ts + c + 1 if strand == '+' else te - c
        return (x - 1) % glen + 1
    c = 0; n_cod = 0
    for n, op in re.findall(r'(\d+)([MIDFGNUV])', cg):
        n = int(n)
        if op in 'MD':
            for k in range(n):
                cd = cod[c+3*k:c+3*k+3]; allc[cd] += 1; n_cod += 1
                occ[cd].append('%s:%d' % (gene, fwd(c+3*k)))
            c += 3*n
        elif op in 'FGNUV': c += n
    print('#', gene, q, strand, '%d-%d' % (fwd(0) if strand=='+' else fwd(te-ts-1), fwd(te-ts-1) if strand=='+' else fwd(0)), 'codons=%d' % n_cod, 'cg=' + cg)
tot = sum(allc.values())
print('total codons', tot)
for cd in ['AGG', 'AGA', 'TGG', 'TGA', 'GGC', 'GGA', 'GGT', 'GGG', 'CTA', 'CTT', 'CTC', 'CTG', 'CGA', 'CGC', 'CGG', 'CGT', 'TGC', 'TGT', 'CCA', 'CCC', 'CCG', 'CCT', 'ATA', 'TAA', 'TAG']:
    print(cd, allc[cd], ' '.join(occ[cd]) if cd in ('AGG', 'CTA', 'CTT', 'CTC', 'CTG', 'CGA', 'CGC', 'CGG', 'CGT', 'TGG') else '')
