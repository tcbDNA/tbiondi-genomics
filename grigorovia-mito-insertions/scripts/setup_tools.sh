#!/usr/bin/env bash
# External tools. Versions are those used for the results in data/.
#
# - DIAMOND 2.2.9 and miniprot 0.18 (r281): Homebrew (`brew install diamond miniprot`).
# - ARAGORN 1.2.41, Infernal 1.1.5 and tRNAscan-SE 2.0.12 (tRNA calls in the Grigorovia follow-up only), built from
#   the authors' source releases into tools/local; checksums below are those of the files used.
# - Grigorovia follow-up tests (candidate/followup/): minimap2 2.31, samtools 1.25, HMMER 3.4, MAFFT 7.526,
#   PAML 4.10.10 and sra-tools 3.4.1 from Homebrew (`brew install minimap2 samtools hmmer mafft paml sratoolkit`);
#   SPAdes 4.3.0 (official macOS arm64 binary) and compleasm 0.2.10 (GitHub tag), fetched below.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p tools/src tools/bin tools/local
cd tools/src

curl -fLO http://www.ansikte.se/ARAGORN/Downloads/aragorn1.2.41.c
echo "92a31cc5c0b0ad16d4d7b01991989b775f07d2815df135fe6e3eab88f5e97f4a  aragorn1.2.41.c" | shasum -a 256 -c
cc -O3 -w -o ../bin/aragorn aragorn1.2.41.c

curl -fLO http://eddylab.org/infernal/infernal-1.1.5.tar.gz
echo "ad4ddae02f924ca7c85bc8c4a79c9f875af8df96aeb726702fa985cbe752497f  infernal-1.1.5.tar.gz" | shasum -a 256 -c
tar -xzf infernal-1.1.5.tar.gz
(cd infernal-1.1.5 && ./configure --prefix="$PWD/../../local" && make -j8 && make install)

curl -fL -o tRNAscan-SE-2.0.12.tar.gz http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz
echo "96fa4af507cd918c1c623763d9260bd6ed055d091662b44314426f6bbf447251  tRNAscan-SE-2.0.12.tar.gz" | shasum -a 256 -c
tar -xzf tRNAscan-SE-2.0.12.tar.gz
# the archive's timestamps make make(1) try to rerun autotools; skip that
(cd tRNAscan-SE-2.0 && ./configure --prefix="$PWD/../../local" \
    && make AUTOHEADER=true ACLOCAL=true AUTOMAKE=true AUTOCONF=true \
    && make install AUTOHEADER=true ACLOCAL=true AUTOMAKE=true AUTOCONF=true)
cd ..
# SPAdes 4.3.0 (official binary; the Homebrew formula was not used)
mkdir -p spades_dl && (cd spades_dl \
    && curl -fLO https://github.com/ablab/spades/releases/download/v4.3.0/SPAdes-4.3.0-Darwin-arm64.tar.gz \
    && echo "2cb4a84dac08ec7c9751ad791b6c9b15f618cc4766da69d05f39014e137fb2f0  SPAdes-4.3.0-Darwin-arm64.tar.gz" | shasum -a 256 -c \
    && tar -xzf SPAdes-4.3.0-Darwin-arm64.tar.gz)
# compleasm 0.2.10 (two files are used: compleasm.py and _version.py)
mkdir -p compleasm_dl compleasm && (cd compleasm_dl \
    && curl -fL -o compleasm-0.2.10.tar.gz https://github.com/huangnengCSU/compleasm/archive/refs/tags/v0.2.10.tar.gz \
    && echo "76d7bdd3259b76b085a57f179df789397ff862b3107fedf6064b167c1442405f  compleasm-0.2.10.tar.gz" | shasum -a 256 -c \
    && tar -xzf compleasm-0.2.10.tar.gz && cp compleasm-0.2.10/compleasm.py compleasm-0.2.10/_version.py ../compleasm/)
echo "tools in tools/bin, tools/local/bin, tools/spades_dl and tools/compleasm"
