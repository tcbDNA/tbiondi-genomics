"""Fetch the 2,695 fungal mitogenome TPA records (Sarhan et al. 2025; BK072095-BK074789) as GenBank flat files.

    uv run python code/audit/fetch_fungal_tpa.py
Writes data/raw/fungal_tpa/BK072095-BK074789.gb.gz (NCBI efetch, 100 records per request, resumable).
"""

from __future__ import annotations

import gzip
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/raw/fungal_tpa"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id="
UA = {"User-Agent": "grigorovia-mito-insertions/0.1"}
FIRST, LAST = 72095, 74789


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ids = [f"BK{n:06d}" for n in range(FIRST, LAST + 1)]
    parts = []
    for i in range(0, len(ids), 100):
        part = OUT / f"part_{i // 100:03d}.gb"
        parts.append(part)
        if part.exists() and part.stat().st_size > 0:
            continue
        for attempt in range(6):
            try:
                req = urllib.request.Request(EFETCH + ",".join(ids[i:i + 100]), headers=UA)
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                n = data.count(b"\n//")
                if n != len(ids[i:i + 100]):
                    print(f"part {i // 100}: {n} records for {len(ids[i:i + 100])} ids", flush=True)
                part.write_bytes(data)
                break
            except Exception as e:  # noqa: BLE001
                print(f"part {i // 100} attempt {attempt}: {e}", flush=True)
                time.sleep(5 * (attempt + 1))
        time.sleep(0.4)
    with gzip.open(OUT / f"BK{FIRST:06d}-BK{LAST:06d}.gb.gz", "wb") as out:
        for p in parts:
            out.write(p.read_bytes())
    total = sum(p.read_bytes().count(b"\n//") for p in parts)
    for p in parts:
        p.unlink()
    print(f"records: {total} of {len(ids)}")


if __name__ == "__main__":
    main()
