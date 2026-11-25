#!/bin/sh
# Gate 3 sweep wrapper, reordered manifest (non-arthropods first, 2026-11-24): restarts after a crash (checkpointed), at most 20 times; stops when FINISHED or on a disk stop.
cd "$(dirname "$0")/.."; export GATE3_UNLOCK=1
n=0; while [ $n -lt 20 ]; do
  env/prefix/bin/python -I src/gate3_run.py results/gate3/manifest_nonarthropoda_first.tsv results/gate3/run 8 >> results/gate3/run/stderr.log 2>&1; rc=$?
  tail -1 results/gate3/run/progress.log | grep -q FINISHED && exit 0
  [ $rc = 2 ] && { echo "$(date '+%F %T') stopped: disk floor" >> results/gate3/run/progress.log; exit 2; }
  n=$((n+1)); echo "$(date '+%F %T') runner exited rc=$rc, restart $n" >> results/gate3/run/progress.log; sleep 10
done; echo "$(date '+%F %T') GAVE UP after 20 restarts" >> results/gate3/run/progress.log; exit 1
