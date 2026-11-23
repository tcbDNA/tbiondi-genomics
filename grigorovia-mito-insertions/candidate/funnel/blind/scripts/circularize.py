#!/usr/bin/env python3
"""Append the first N bases of a (circular) genome to its end so genes crossing the
sequence end are contiguous. Usage: circularize.py genome.fna N > genome_circ.fna"""
import sys
name = None; s = []
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]
    else: s.append(l.upper())
s = ''.join(s); n = int(sys.argv[2])
print('>%s\n%s' % (name, s + s[:n]))
