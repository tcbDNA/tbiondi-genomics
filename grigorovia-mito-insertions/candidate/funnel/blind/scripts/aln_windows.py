#!/usr/bin/env python3
"""Print windows of miniprot --aln blocks around frameshift columns.
Usage: aln_windows.py file.aln [flank_columns]"""
import sys
flank = int(sys.argv[2]) if len(sys.argv) > 2 else 45
blk = None
def show(b):
    if not b or 'ATN' not in b: return
    atn, ata, aas, aqa = b['ATN'], b['ATA'], b['AAS'], b['AQA']
    paf = b['PAF'].split('\t')
    # frameshift columns: in ATN a lowercase/odd column; miniprot marks frameshift nt where AQA has no residue triple.
    # Find runs in ATN where ATA shows '$' or '!' (miniprot uses '$' for frameshift, '!' for stop)
    cols = [i for i, ch in enumerate(ata) if ch in '$!']
    cols = sorted(set(cols))
    done = set()
    for c in cols:
        if any(abs(c - d) < 6 for d in done): continue
        done.add(c)
        lo, hi = max(0, c - flank), c + flank
        print('#', paf[0], paf[5], paf[4], 'qspan', paf[2], paf[3], 'tspan', paf[7], paf[8], 'col', c)
        for tag, s in (('ATN', atn), ('ATA', ata), ('AAS', aas), ('AQA', aqa)):
            print('%s %s' % (tag, s[lo:hi]))
        print()
for line in open(sys.argv[1]):
    line = line.rstrip('\n')
    if line and not line.startswith('#'):
        show(blk); blk = {'PAF': line}
    elif line.startswith('##') and blk is not None:
        tag, s = line[2:].split('\t', 1)
        blk[tag] = s
show(blk)
