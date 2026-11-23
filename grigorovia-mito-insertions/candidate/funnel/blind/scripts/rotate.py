#!/usr/bin/env python3
"""Rotate a circular genome so it starts at its midpoint. Usage: rotate.py genome.fna > rot.fna"""
import sys
name = None; s = []
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]
    else: s.append(l.upper())
s = ''.join(s); h = len(s) // 2
print('>%s_rot_offset%d\n%s' % (name, h, s[h:] + s[:h]))
