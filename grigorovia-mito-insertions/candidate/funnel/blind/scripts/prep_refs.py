#!/usr/bin/env python3
"""Give duplicate reference headers unique names (refA_cob x3 -> refA_cob_v1/_v2/_v3).
Usage: prep_refs.py reference_proteins.faa > refs_renamed.faa"""
import sys
seqs = []
for l in open(sys.argv[1]):
    l = l.strip()
    if not l: continue
    if l.startswith('>'): seqs.append([l[1:], ''])
    else: seqs[-1][1] += l
cnt = {}
for n, s in seqs:
    if sum(1 for x, _ in seqs if x == n) > 1:
        cnt[n] = cnt.get(n, 0) + 1; n = '%s_v%d' % (n, cnt[n])
    print('>%s\n%s' % (n, s))
