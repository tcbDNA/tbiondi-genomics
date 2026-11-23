"""Final reference protein set and DIAMOND database, derived from prep.py's ref_proteins.faa/.tsv.

    uv run python code/refs.py [--length-filter]

Rules (DESIGN.md amendment 4, PREREG.md section 1):
  - EXCLUDE: reference genomes removed from the references (NC_079576, "Pholiota spumosa", is a
    mislabeled ascomycete; Strophariaceae has no other RefSeq genome with proteins, so it gets no
    representative).
  - --length-filter: drop core-gene reference proteins longer than LENGTH_FACTOR x the median length of
    that gene among references of the same kingdom group (Metazoa, Fungi, Viridiplantae, other). Such
    proteins usually carry intron-encoded sequence fused in frame.
Writes data/work/ref_final.faa, ref_final.dmnd and ref_final.json (what was dropped).
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "data/work"
EXCLUDE = {"NC_079576.1"}
LENGTH_FACTOR = 1.5
CORE = {"cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6"}


def group(meta: dict) -> str:
    k = meta["kingdom_name"]
    return k if k in ("Metazoa", "Fungi", "Viridiplantae") else "other"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--length-filter", action="store_true")
    args = ap.parse_args()
    G = {r["accession"]: r for r in csv.DictReader(open(WORK / "genomes.tsv"), delimiter="\t")}
    R = [r for r in csv.DictReader(open(WORK / "ref_proteins.tsv"), delimiter="\t")]
    drop = {r["ref_id"]: "excluded_genome" for r in R if r["accession"] in EXCLUDE}
    if args.length_filter:
        lens = defaultdict(list)
        for r in R:
            if r["gene"] in CORE and r["ref_id"] not in drop:
                lens[(r["gene"], group(G[r["accession"]]))].append(int(r["length"]))
        med = {k: statistics.median(v) for k, v in lens.items()}
        for r in R:
            k = (r["gene"], group(G[r["accession"]]))
            if r["gene"] in CORE and r["ref_id"] not in drop and int(r["length"]) > LENGTH_FACTOR * med[k]:
                drop[r["ref_id"]] = f"long:{r['length']}>{LENGTH_FACTOR}x{med[k]:.0f}"
    n = 0
    with open(WORK / "ref_proteins.faa") as fi, open(WORK / "ref_final.faa", "w") as fo:
        keep = False
        for line in fi:
            if line.startswith(">"):
                keep = line[1:].split()[0] not in drop
                n += keep
            if keep:
                fo.write(line)
    subprocess.run(["diamond", "makedb", "--in", str(WORK / "ref_final.faa"), "-d", str(WORK / "ref_final"), "--quiet"],
                   check=True)
    by_kind = defaultdict(int)
    for v in drop.values():
        by_kind[v.split(":")[0]] += 1
    info = {"length_filter": args.length_filter, "length_factor": LENGTH_FACTOR if args.length_filter else None,
            "exclude": sorted(EXCLUDE), "n_kept": n, "n_dropped": len(drop), "dropped_by_kind": by_kind,
            "dropped": drop}
    json.dump(info, open(WORK / "ref_final.json", "w"), indent=1)
    print(f"ref_final: {n} proteins kept, {len(drop)} dropped {dict(by_kind)}; length_filter={args.length_filter}")


if __name__ == "__main__":
    main()
