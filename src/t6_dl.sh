#!/bin/sh
# cold-review T6 reads: sequential downloads (dl.sh enforces the 30 GB floor and ledgers)
cd "$(dirname "$0")/.."
for r in SRR38145622 SRR13036352 SRR33975184; do sh src/dl.sh https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r results/work/$r.sra 5; done; echo ALL_DL_DONE
