#!/usr/bin/env python3
"""Coding-strand 81-nt context (site in upper case) and its GC fraction.
Usage: local_context.py genome.fna POS STRAND [POS STRAND ...]"""
import sys
COMP = str.maketrans('ACGT', 'TGCA')
s = ''.join(l.strip() for l in open(sys.argv[1]) if not l.startswith('>')).upper()
L = len(s); s3 = s + s + s; a = sys.argv[2:]
for i in range(0, len(a), 2):
    pos, st = int(a[i]), a[i+1]
    w = s3[L+pos-1-40:L+pos+40]
    if st == '-': w = w.translate(COMP)[::-1]
    gc = sum(1 for x in w if x in 'GC') / len(w)
    print('%d%s GC=%.2f %s%s%s' % (pos, st, gc, w[:40].lower(), w[40], w[41:].lower()))
