#!/bin/sh
# dl.sh URL OUTFILE GATE : resumable download, refuses if <30 GB would stay free, appends to the ledger
cd /Users/tcb_dna/projects/random_biology_stuff/metazoan-uaa-tyr-check
sz=$(curl -sSIL "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{s=$2} END{print s+0}'); free=$(df -k . | tail -1 | awk '{print $4*1024}')
[ $(echo "$free - $sz" | awk '{print ($1-$3 < 30000000000)}') = 1 ] && { echo "REFUSED: <30 GB would remain ($2)"; exit 2; }
curl -sSL --fail --retry 5 -C - -o "$2" "$1" || { echo "FAIL $2"; exit 1; }
got=$(stat -f %z "$2"); [ "$sz" -gt 0 ] && [ "$got" != "$sz" ] && { echo "SIZE_MISMATCH $2 $got $sz"; exit 3; }
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$3" "$1" "$2" "$got" "$(shasum -a 256 "$2" | cut -d' ' -f1)" "on_disk" >> data/manifests/download_ledger.tsv; echo "OK $2 $got"
