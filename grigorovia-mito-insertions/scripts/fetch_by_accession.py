"""Fetch nucleotide records by accession.version from NCBI E-utilities into a gzipped GenBank file and a gzipped FASTA.

    uv run python scripts/fetch_by_accession.py --accessions list.txt --gbff out.gbff.gz --fasta out.fna.gz

No API key: 100 records per request, pauses between requests, retries on errors. Resumable: records already in the
output files are skipped.
"""

from __future__ import annotations

import argparse
import gzip
import time
import urllib.request
from pathlib import Path

from Bio import SeqIO

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = {"User-Agent": "grigorovia-mito-insertions/0.1"}


def fetch(ids: list[str], rettype: str = "gbwithparts") -> str:
    url = f"{EFETCH}?db=nuccore&rettype={rettype}&retmode=text&id={','.join(ids)}"
    for attempt in range(6):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                return r.read().decode()
        except Exception:  # noqa: BLE001
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"failed: {ids[0]}..{ids[-1]}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--fasta", required=True)
    args = ap.parse_args()
    accs = [l.strip() for l in open(args.accessions) if l.strip()]
    done = set()
    if Path(args.gbff).exists():
        with gzip.open(args.gbff, "rt") as fh:
            done = {rec.id for rec in SeqIO.parse(fh, "genbank")}
    todo = [a for a in accs if a not in done]
    with gzip.open(args.gbff, "at") as gb, gzip.open(args.fasta, "at") as fa:
        for i in range(0, len(todo), 100):
            batch = todo[i:i + 100]
            gb.write(fetch(batch, "gbwithparts"))
            time.sleep(0.5)
            fa.write(fetch(batch, "fasta"))   # FASTA separately: CON records carry no sequence in GenBank format
            time.sleep(0.5)
            print(f"{min(i + 100, len(todo))}/{len(todo)}", flush=True)


if __name__ == "__main__":
    main()
