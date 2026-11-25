#!/bin/sh
# Post-hoc read checks, smallest first; each run is deleted after use. Output: results/gate3/reads_b_<acc>.txt
cd "$(dirname "$0")/.."; P="env/prefix/bin/python -I src/gate3b_reads.py"; O=results/gate3
r(){ a=$1; shift; echo "$(date '+%F %T') start $a $*" >> $O/reads_b.log; nice -n 10 $P $a "$@" > $O/reads_b_$a.txt 2>&1; echo "$(date '+%F %T') done $a" >> $O/reads_b.log; }
r NC_013436.1 SRR22396609
r KU726971.1 SRR15202050 SRR15202047
r PV357903.1 SRR25465297
r PX804232.1 SRR36182985
r NC_013627.1 SRR36227421
r NC_013628.2 SRR38094758
r CM163399.1 SRR37037826
r PZ495457.1 SRR38368693
r NC_047225.1 SRR9586127 SRR9586188
r OZ318138.1 SRR33975184 ERR14792883
r QB209682.1 SRR41519539
r OZ407450.1 ERR14857142 ERR14842678
r OZ318138.1 ERR14835550
echo "$(date '+%F %T') ALL DONE" >> $O/reads_b.log
