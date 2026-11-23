#!/usr/bin/env python3
"""Show, for a reference protein and residue window, the codons each genome aligns to
those residues (from miniprot cg strings), so the three genomes can be compared at each site.
Usage: compare_sites.py refs.faa REF START END paf1:circ1:len1 [paf2:circ2:len2 ...]
START/END are 1-based reference residue positions. '>' marks extra genome nt (F op)."""
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
refs = readfa(sys.argv[1]); R = sys.argv[2]; a, b = int(sys.argv[3]), int(sys.argv[4])
print('%-10s' % 'ref', ' '.join('%4s' % ('%d%s' % (i, refs[R][i-1]))[-4:] for i in range(a, b+1)))
for spec in sys.argv[5:]:
    paf, circ, glen = spec.split(':'); glen = int(glen)
    G = readfa(circ)
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12 or f[0] != R or f[4] not in '+-': continue
        qs, strand, t, ts, te = int(f[2]), f[4], f[5], int(f[7]), int(f[8])
        cg = next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
        cod = G[t][ts:te] if strand == '+' else rc(G[t][ts:te])
        c = 0; qi = qs; m = {}; extra = {}
        for n, op in re.findall(r'(\d+)([MIDFGNUV])', cg):
            n = int(n)
            if op == 'M':
                for k in range(n): m[qi+k] = cod[c+3*k:c+3*k+3]
                c += 3*n; qi += n
            elif op == 'I':
                for k in range(n): m[qi+k] = '---'
                qi += n
            elif op == 'D': c += 3*n
            elif op == 'F': extra[qi] = cod[c:c+n].lower(); c += n
            elif op == 'G': m[qi] = cod[c:c+n] + '?'; c += n; qi += 1
            elif op in 'NUV': c += n; qi += (0 if op == 'N' else 1)
        row = []
        for i in range(a-1, b):
            s = m.get(i, '   ')
            if i in extra: s = '>' + extra[i] + s
            row.append('%4s' % s)
        print('%-10s' % t, ' '.join(row))
