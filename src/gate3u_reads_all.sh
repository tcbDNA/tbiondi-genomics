#!/bin/sh
cd "$(dirname "$0")/.."; O=results/gate3_unverified/reads_hits.txt; : > $O
r(){ echo "== $1" >> $O; nice -n 5 env/prefix/bin/python -I src/gate3u_reads.py "$@" >> $O 2>&1; }
r PV962435.1 SRR10079992
r PV962563.1 SRR10080023
r PV962571.1 SRR11014335 SRR10079898
r KM893424.1 SRR10251503
r MW776437.1 SRR19044081
r MW255583.1 SRR14069289
r MW810062.1 SRR8998268
echo "ALL DONE" >> $O
