#!/usr/bin/env python3
"""Build frameshift-free coding sequences (as read by the ribosome under the canonical placements)
for the 7 genes of one genome from miniprot PAF, then report (a) base composition, (b) per-position
frequency of a motif in the genomic CDS span, (c) in-frame codon counts for listed codons.
Usage: motif_null.py genome_circ.fna glen refs.faa paf MOTIF[,MOTIF...] CODON[,CODON...]"""
import re, sys
from collections import Counter
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
G = readfa(sys.argv[1]); glen = int(sys.argv[2]); paf = sys.argv[4]
motifs = sys.argv[5].split(','); codons = sys.argv[6].split(',')
best = {}
for line in open(paf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[4] not in '+-' or '_v2' in f[0] or '_v3' in f[0]: continue
    gene = f[0].split('_', 1)[1].replace('_v1', '')
    ms = int(next(x[5:] for x in f[12:] if x.startswith('ms:i:')))
    if gene not in best or ms > best[gene][0]: best[gene] = (ms, f)
span_seq = ''   # genomic coding-strand sequence of exons incl. extra bases
comp = Counter(); mot = Counter(); npos = 0
for gene, (ms, f) in best.items():
    strand, t, ts, te = f[4], f[5], int(f[7]), int(f[8])
    cg = next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
    cod = G[t][ts:te] if strand == '+' else rc(G[t][ts:te])
    c = 0; exon = ''
    for n, op in re.findall(r'(\d+)([MIDFGNUV])', cg):
        n = int(n)
        if op in 'MD': exon += cod[c:c+3*n]; c += 3*n
        elif op in 'FG': exon += cod[c:c+n]; c += n
        elif op == 'N': c += n
        elif op in 'UV': c += n   # phase-1/2 intron: skipped (loses one split codon)
    comp.update(exon)
    for i in range(len(exon)):
        npos += 1
        for m in motifs:
            if exon[i:i+len(m)] == m: mot[m] += 1
tot = sum(comp.values())
print('coding-strand composition:', {b: round(comp[b]/tot, 3) for b in 'ACGT'}, 'n=%d' % tot)
for m in motifs: print('motif %s: %d occurrences in %d positions (%.4f per position)' % (m, mot[m], npos, mot[m]/npos))
