#!/bin/sh
# EXPLORATORY/POST-HOC: identical gate3_run.py pipeline on the UNVERIFIED set, separate output directory. Restarts up to 20 times (checkpointed).
cd /Users/tcb_dna/projects/random_biology_stuff/metazoan-uaa-tyr-check; export GATE3_UNLOCK=1; D=results/gate3_unverified/run
n=0; while [ $n -lt 20 ]; do
  env/prefix/bin/python -I src/gate3_run.py results/gate3_unverified/manifest.tsv $D 8 >> $D/stderr.log 2>&1; rc=$?
  tail -1 $D/progress.log | grep -q FINISHED && exit 0
  [ $rc = 2 ] && { echo "$(date '+%F %T') stopped: disk floor" >> $D/progress.log; exit 2; }
  n=$((n+1)); echo "$(date '+%F %T') runner exited rc=$rc, restart $n" >> $D/progress.log; sleep 10
done; echo "$(date '+%F %T') GAVE UP after 20 restarts" >> $D/progress.log; exit 1
