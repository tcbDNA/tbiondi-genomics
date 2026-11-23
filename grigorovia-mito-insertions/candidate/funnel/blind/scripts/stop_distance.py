#!/usr/bin/env python3
"""How many codons the ORIGINAL (unshifted) frame runs from the shift codon to the first TAA/TAG.
Usage: stop_distance.py genome.fna POS STRAND [POS STRAND ...] (POS = canonical extra base)"""
import sys
COMP = str.maketrans('ACGT', 'TGCA')
s = ''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
L = len(s); s3 = s + s + s; a = sys.argv[2:]
for i in range(0, len(a), 2):
    pos, st = int(a[i]), a[i+1]
    d = s3[L+pos-1:L+pos-1+900] if st == '+' else s3[L+pos-900:L+pos].translate(COMP)[::-1]
    n = next((j//3 for j in range(0, len(d)-2, 3) if d[j:j+3] in ('TAA', 'TAG')), None)
    print('%d%s first in-frame stop of the original frame at codon %s (shift codon = 1)' % (pos, st, None if n is None else n + 1))
