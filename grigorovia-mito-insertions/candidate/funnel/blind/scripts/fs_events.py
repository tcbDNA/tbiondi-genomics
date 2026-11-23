#!/usr/bin/env python3
"""Extract frameshift events (cg ops F and G) from miniprot PAF output.
Usage: fs_events.py genome.fna refs.faa miniprot.paf [genome_len]
Walks the cg CIGAR in protein order, tracks genome position on the coding strand,
and reports each frameshift with forward-strand 1-based coordinates, the coding-strand
context in the original reading frame and the reference residue aligned right after."""
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

# yeast mito code (table 3)
B = 'TCAG'
AA3 = 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a+b+c: AA3[16*i+4*j+k] for i,a in enumerate(B) for j,b in enumerate(B) for k,c in enumerate(B)}
def tr(s): return ''.join(CODE.get(s[i:i+3], 'X') for i in range(0, len(s)-2, 3))

genome = readfa(sys.argv[1]); refs = readfa(sys.argv[2])
GLEN = int(sys.argv[4]) if len(sys.argv) > 4 else None
for line in open(sys.argv[3]):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[4] not in ('+', '-'): continue   # skip unmapped (-u) lines
    q, qlen, qs, qe, strand, t, tlen, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8])
    cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), '')
    g = genome[t]
    # coding-strand coordinate system: c = offset along coding strand from alignment start
    def fwd(c):  # coding offset -> forward 1-based coordinate (wrapped if genome_len given)
        x = ts + c + 1 if strand == '+' else te - c
        return (x - 1) % GLEN + 1 if GLEN else x
    cod = g[ts:te] if strand == '+' else rc(g[ts:te])
    c = 0; qi = qs
    for n, op in re.findall(r'(\d+)([MIDFGNUV])', cg):
        n = int(n)
        if op == 'M': c += 3*n; qi += n
        elif op == 'I': qi += n
        elif op == 'D': c += 3*n
        elif op in 'NUV':
            c += n; qi += (0 if op == 'N' else 1)
        elif op in 'FG':
            if op == 'F':   # n extra nt in genome, no residue consumed
                kind = '+%d (extra %d nt)' % (n, n); shift = n; nxt = qi
            else:           # one residue aligned to n (<3) nt
                kind = '-%d (codon of %d nt)' % (3-n, n); shift = n - 3; nxt = qi
            start_c = c
            up = cod[max(0, start_c-12):start_c]
            # original frame codons across the site (frame set by preceding codons)
            after = cod[start_c:start_c+9]
            orig = ' '.join(up[i:i+3] for i in range(0, len(up), 3)) + ' | ' + ' '.join(after[i:i+3] for i in range(0, len(after), 3))
            ext = cod[start_c:start_c+n]
            if strand == '+': coord = '%d-%d' % (fwd(start_c), fwd(start_c+n-1)) if n > 1 else str(fwd(start_c))
            else: coord = '%d-%d' % (fwd(start_c+n-1), fwd(start_c)) if n > 1 else str(fwd(start_c))
            r = refs[q]
            nxt_aa = r[nxt] if nxt < len(r) else '-'
            print('\t'.join(map(str, [t, q, strand, op, n, kind, coord, ext, 'aa_after=%d%s' % (nxt+1, nxt_aa), orig])))
            c += n
            if op == 'G': qi += 1
