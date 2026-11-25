#!/bin/sh
# after calibration: run each leave-one-order-out CM on the whole panel, glocal
cd "$(dirname "$0")/.."; D=results/gate2/trny_diag
until grep -q "calibration done" $D/log.txt; do sleep 5; done
for c in nem_minus_Mermithida nem_minus_Dorylaimida nem_minus_Trichinellida nem_minus_Rhabditida nem_all; do
  cols=$(env/prefix/bin/cmemit -c $D/$c.cm | grep -v '>' | tr -d '\n' | env/prefix/bin/python -I -c "
import sys,re
s=sys.stdin.read().upper().replace('T','U'); m=re.search(r'UGUA[AG][AG]',s) or re.search(r'GUA',s); g=m.start()+(2 if m.group(0)[0]=='U' else 1)   # 1-based col of G34
print(f'({g-1},{g},range({g-7},{g-2}),range({g+5},{g+10}),range({g-2},{g+5}))')")
  echo "$(date '+%T') C $c cols $cols" >> $D/log.txt
  env/prefix/bin/python -I src/trny_variants.py C $D/$c.cm $D/C_$c.jsonl 6 results/gate2/trny_diag/ids_all.txt "$cols" >> $D/log.txt 2>&1
done; echo "$(date '+%T') C end" >> $D/log.txt
