#!/usr/bin/env python3
"""For each frameshift (cg op F/G) in a miniprot PAF, test every placement of the
extra/missing base(s) within +-W nt of miniprot's call, and report which placements
give the best identity to the reference protein in a fixed ungapped window.
Usage: fs_window.py genome_circ.fna genome_len refs.faa miniprot.paf [W]
Coordinates are reported 1-based on the forward strand of the original genome
(positions beyond genome_len are wrapped)."""
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
B = 'TCAG'
AA3 = 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a+b+c: AA3[16*i+4*j+k] for i,a in enumerate(B) for j,b in enumerate(B) for k,c in enumerate(B)}
def tr(s): return ''.join(CODE.get(s[i:i+3], 'X') for i in range(0, len(s)-2, 3))

gfa, glen, rfa, paf = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
W = int(sys.argv[5]) if len(sys.argv) > 5 else 18
genome = readfa(gfa); refs = readfa(rfa)
for line in open(paf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 12 or f[4] not in '+-': continue
    q, qs, strand, t, ts, te = f[0], int(f[2]), f[4], f[5], int(f[7]), int(f[8])
    cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), '')
    if not re.search(r'\d+[FG]', cg): continue
    g = genome[t]; r = refs[q]
    cod = g[ts:te] if strand == '+' else rc(g[ts:te])
    def fwd(c):
        x = ts + c + 1 if strand == '+' else te - c
        return (x - 1) % glen + 1
    # build residue -> codon-start map along coding strand, ignoring frameshift ops
    c = 0; qi = qs; res2c = {}; events = []
    for n, op in re.findall(r'(\d+)([MIDFGNUV])', cg):
        n = int(n)
        if op == 'M':
            for k in range(n): res2c[qi + k] = c + 3*k
            c += 3*n; qi += n
        elif op == 'I': qi += n
        elif op == 'D': c += 3*n
        elif op in 'NUV': c += n; qi += (0 if op == 'N' else 1)
        elif op == 'F': events.append((c, n, qi, 'F')); c += n
        elif op == 'G': events.append((c, n, qi, 'G')); c += n; qi += 1
    for ei, (c0, n, qi0, op) in enumerate(events):
        prev_q = events[ei-1][2] if ei > 0 else -1
        next_q = events[ei+1][2] if ei + 1 < len(events) else 10**9
        # number of nt to delete to restore frame: F n -> delete n ; G n -> insert 3-n (we test deletion of n for F only)
        if op == 'G':
            print('\t'.join([t, q, strand, 'G%d' % n, str(fwd(c0)), 'aa_after=%d%s' % (qi0+1, r[qi0]), 'G-op (missing base) not window-tested']))
            continue
        # reference residues upstream (aligned in frame0) and downstream (frame shifted)
        up_res = [k for k in range(qi0-1, qi0-9, -1) if k in res2c and k >= prev_q]
        dn_res = [k for k in range(qi0, qi0+8) if k in res2c and k < next_q]
        if not up_res or not dn_res: continue
        lo_res, hi_res = min(up_res), max(dn_res)
        start = res2c[lo_res]                   # coding-strand start of first window codon
        end = res2c[hi_res] + 3 + n             # (dn codon coords already include the shift)
        # fix: res2c for downstream residues was recorded after the shift, so they already include +n
        end = res2c[hi_res] + 3
        seg = cod[start:end]
        refseg = r[lo_res:hi_res+1]
        best = []; scores = {}
        for p in range(start, end - n):
            s2 = cod[start:p] + cod[p+n:end]
            pep = tr(s2)
            m = sum(1 for a, b in zip(pep, refseg) if a == b)
            stops = pep.count('*')
            scores[p] = (m, -stops, pep)
        mx = max(scores.values())[:2]
        okp = sorted(p for p, v in scores.items() if v[:2] == mx)
        coords = sorted(set(fwd(p) for p in okp))
        print('\t'.join(map(str, [t, q, strand, 'F%d' % n, 'miniprot_pos=%d' % fwd(c0), 'aa_after=%d%s' % (qi0+1, r[qi0]),
              'ref=%s' % refseg, 'best_id=%d/%d' % (mx[0], len(refseg)),
              'best_placements(fwd)=%s' % ','.join(map(str, coords)),
              'peps=%s' % ','.join(sorted(set(scores[p][2] for p in okp)))])))
