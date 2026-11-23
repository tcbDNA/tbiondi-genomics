#!/usr/bin/env python3
"""Length and base composition of each genome. Usage: composition.py genome.fna [...]"""
import sys
from collections import Counter
for p in sys.argv[1:]:
    s = ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()
    c = Counter(s)
    print(p, len(s), {b: c[b] for b in 'ACGT'}, 'GC=%.3f' % ((c['G'] + c['C']) / len(s)))
