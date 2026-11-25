#!/bin/sh
# calibrate the cross-validated nematode CMs, then run each fold's model on the WHOLE panel (variant C)
cd "$(dirname "$0")/.."; D=results/gate2/trny_diag
for c in nem_all nem_minus_Mermithida nem_minus_Dorylaimida nem_minus_Trichinellida nem_minus_Rhabditida; do echo "$(date '+%T') calibrate $c" >> $D/log.txt; env/prefix/bin/cmcalibrate --cpu 4 $D/$c.cm >> $D/calib.log 2>&1 || echo "FAIL $c" >> $D/log.txt; done; echo "$(date '+%T') calibration done" >> $D/log.txt
