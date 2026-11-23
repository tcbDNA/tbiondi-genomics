#!/usr/bin/env python3
"""List every occurrence of a motif inside the 7 standard genes (best miniprot alignment per gene),
with its position relative to the reading frame actually supported by the alignment.
codon_pos 0/1/2 = motif starts at 1st/2nd/3rd base of a codon in the aligned frame;
'x' = motif starts on a base miniprot treats as an extra (frameshift) nucleotide.
Usage: motif_scan.py genome_circ.fna glen refs.faa paf MOTIF"""
import re, sys
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
s = ''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
glen = int(sys.argv[2]); paf = sys.argv[4]; motif = sys.argv[5]
best = {}
for line in open(paf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[4] not in '+-' or '_v2' in f[0] or '_v3' in f[0]: continue
    gene = f[0].split('_', 1)[1].replace('_v1', '')
    ms = int(next(x[5:] for x in f[12:] if x.startswith('ms:i:')))
    if gene not in best or ms > best[gene][0]: best[gene] = (ms, f)
summary = {}
for gene, (ms, f) in sorted(best.items()):
    strand, ts, te = f[4], int(f[7]), int(f[8])
    cg = next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
    cod = s[ts:te] if strand == '+' else rc(s[ts:te])
    def fwd(c):
        x = ts + c + 1 if strand == '+' else te - c
        return (x - 1) % glen + 1
    c = 0; pos = {}
    for n, op in re.findall(r'(\d+)([MIDFGNUV])', cg):
        n = int(n)
        if op in 'MD':
            for k in range(3*n): pos[c+k] = k % 3
            c += 3*n
        elif op == 'F':
            for k in range(n): pos[c+k] = 'x'
            c += n
        elif op in 'GNUV': c += n
        # 'I' consumes no genome nt
    for m in re.finditer('(?=%s)' % motif, cod):
        p = m.start()
        if p not in pos: continue
        summary[pos[p]] = summary.get(pos[p], 0) + 1
        print('%s\t%s\t%s\tfwd=%d\tcodon_pos=%s\tctx=%s%s%s' % (f[5], gene, motif, fwd(p), pos[p],
              cod[max(0, p-6):p].lower(), cod[p:p+len(motif)], cod[p+len(motif):p+len(motif)+6].lower()))
print('# summary codon_pos counts:', summary)
