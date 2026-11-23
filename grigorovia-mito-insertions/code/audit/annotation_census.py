"""Prior-art audit: how do GenBank/RefSeq mitochondrial records already write frameshifts?

Reads annotation only (feature tables), never scans sequence. For each record, counts CDS features with:
  ribosomal_slippage qualifier; exception qualifier (any value); a note mentioning a frameshift;
  artificial_location (NCBI's fix for low-quality sequence); pseudo/pseudogene; and CDS whose join() has
  adjacent parts separated by -2..+2 nt (a frameshift written as a join). Also counts intron features
  and "tiny" introns (<= 10 nt), and CDS joins whose gaps are 3..10 nt.

    uv run python code/audit/annotation_census.py <gbff[.gz]> <label>
Writes results/audit/annotation_<label>.tsv (one row per record) and prints a summary.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

ROOT = Path(__file__).resolve().parents[2]
FS_NOTE = re.compile(r"frame.?shift|slippage|frameshifted", re.I)


def group(lineage: list[str]) -> str:
    s = set(lineage)
    if "Embryophyta" in s:
        return "land_plants"
    for g in ("Metazoa", "Fungi", "Viridiplantae", "Rhodophyta", "Sar", "Discoba", "Amoebozoa"):
        if g in s:
            return g
    return "other_eukaryotes"


def join_gaps(loc) -> list[int]:
    """Gaps between consecutive parts in transcript order (negative = overlap)."""
    if not isinstance(loc, CompoundLocation):
        return []
    parts = list(loc.parts)
    gaps = []
    for a, b in zip(parts, parts[1:]):
        if a.strand != b.strand or a.ref or b.ref:
            gaps.append(10_000)  # trans-spliced or on another record
        elif a.strand == -1:
            gaps.append(int(a.start) - int(b.end))  # b lies upstream on the forward strand
        else:
            gaps.append(int(b.start) - int(a.end))
    return gaps


def main() -> None:
    path, label = Path(sys.argv[1]), sys.argv[2]
    out = ROOT / "results/audit"
    out.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    cols = ["accession", "group", "organism", "length", "n_cds", "slippage", "exception", "note_fs",
            "artificial_location", "pseudo", "join_shift", "join_gap3_10", "introns", "tiny_introns",
            "exception_values", "join_shift_genes"]
    rows = []
    with opener(path, "rt") as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            lineage = rec.annotations.get("taxonomy", [])
            r = {c: 0 for c in cols}
            r.update(accession=rec.id, group=group(lineage), organism=rec.annotations.get("organism", ""),
                     length=len(rec.seq), exception_values=set(), join_shift_genes=[])
            for f in rec.features:
                q = f.qualifiers
                if f.type == "intron":
                    r["introns"] += 1
                    if len(f.location) <= 10:
                        r["tiny_introns"] += 1
                if f.type not in ("CDS", "gene"):
                    continue
                if f.type == "CDS":
                    r["n_cds"] += 1
                    if "ribosomal_slippage" in q:
                        r["slippage"] += 1
                    if "exception" in q:
                        r["exception"] += 1
                        r["exception_values"].update(q["exception"])
                    if "artificial_location" in q:
                        r["artificial_location"] += 1
                    gaps = join_gaps(f.location)
                    if any(-2 <= g <= 2 for g in gaps):
                        r["join_shift"] += 1
                        r["join_shift_genes"].append(q.get("gene", q.get("product", ["?"]))[0])
                    if any(3 <= g <= 10 for g in gaps):
                        r["join_gap3_10"] += 1
                    if any(FS_NOTE.search(n) for n in q.get("note", [])):
                        r["note_fs"] += 1
                if "pseudo" in q or "pseudogene" in q:
                    r["pseudo"] += 1
            r["exception_values"] = "|".join(sorted(r["exception_values"]))
            r["join_shift_genes"] = ",".join(r["join_shift_genes"])
            rows.append(r)
    with open(out / f"annotation_{label}.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    def summarise(sub: list[dict], name: str) -> None:
        n = len(sub)
        any_fs = sum(1 for r in sub if r["slippage"] or r["note_fs"] or r["join_shift"] or r["tiny_introns"])
        print(f"{name:28s} records {n:6d} | any frameshift-type annotation {any_fs:5d} | "
              f"slippage {sum(r['slippage'] > 0 for r in sub):4d} | note 'frameshift' {sum(r['note_fs'] > 0 for r in sub):4d} | "
              f"join -2..+2 {sum(r['join_shift'] > 0 for r in sub):4d} | tiny intron {sum(r['tiny_introns'] > 0 for r in sub):4d} | "
              f"exception {sum(r['exception'] > 0 for r in sub):4d} | artificial_loc {sum(r['artificial_location'] > 0 for r in sub):4d} | "
              f"pseudo {sum(r['pseudo'] > 0 for r in sub):4d} | join gap 3-10 {sum(r['join_gap3_10'] > 0 for r in sub):4d}")

    print(f"== {label}")
    summarise(rows, "all")
    for g in sorted({r["group"] for r in rows}):
        summarise([r for r in rows if r["group"] == g], g)
    ex = {}
    for r in rows:
        for v in filter(None, r["exception_values"].split("|")):
            ex[v] = ex.get(v, 0) + 1
    print("exception values (records):", ex)


if __name__ == "__main__":
    main()
