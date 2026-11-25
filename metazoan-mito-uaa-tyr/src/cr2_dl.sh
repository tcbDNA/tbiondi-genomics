#!/bin/sh
# Cold review 2 (post-hoc, 2026-11-25): download reads listed in $1 (one run per line) with src/dl.sh (30 GB floor, ledger). Progress in results/work/r2/dl.log
cd "$(dirname "$0")/.."
while read r; do [ -s results/work/$r.sra ] || sh src/dl.sh https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r results/work/$r.sra 6x; echo "$(date '+%F %T') done $r free $(df -g . | tail -1 | awk '{print $4}') GB"; done < "$1"
echo ALL_DL_DONE "$1"
