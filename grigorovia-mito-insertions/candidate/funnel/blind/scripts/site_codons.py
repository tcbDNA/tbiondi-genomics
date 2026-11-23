#!/usr/bin/env python3
"""Print original-frame and shifted-frame codons around a +1 frameshift site.
Usage: site_codons.py genome.fna POS STRAND [POS STRAND ...]
POS = 1-based forward coordinate of the extra base, taken as the first base of the
original-frame codon at the shift (canonical placement)."""
import sys
COMP = str.maketrans('ACGT', 'TGCA')
B = 'TCAG'; AA3 = 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a+b+c: AA3[16*i+4*j+k] for i,a in enumerate(B) for j,b in enumerate(B) for k,c in enumerate(B)}
s = ''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
L = len(s); s3 = s + s + s
args = sys.argv[2:]
for i in range(0, len(args), 2):
    pos, st = int(args[i]), args[i+1]
    if st == '+':
        up = s3[L+pos-1-9:L+pos-1]; dn = s3[L+pos-1:L+pos-1+12]
    else:
        # coding strand runs toward lower coordinates
        up = s3[L+pos:L+pos+9].translate(COMP)[::-1]
        dn = s3[L+pos-12:L+pos].translate(COMP)[::-1]
    orig = [up[j:j+3] for j in range(0, 9, 3)] + [dn[j:j+3] for j in range(0, 12, 3)]
    shifted = [dn[1+j:4+j] for j in range(0, 9, 3)]
    fmt = lambda cs: ' '.join('%s(%s)' % (c, CODE[c]) for c in cs)
    print('%s%s  original frame: %s [%s] %s   | after +1 shift (extra %s skipped): %s' % (
        pos, st, fmt(orig[:3]), fmt(orig[3:4]), fmt(orig[4:]), dn[0], fmt(shifted)))
