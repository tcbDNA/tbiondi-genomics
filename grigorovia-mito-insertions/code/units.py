"""Assembly/strain units (PREREG log: "unit for the strong rule"). Records sharing a BioSample or an Assembly
accession (RefSeq DBLINK; fungal TPA: BioSample from the Sarhan accession table), or with an identical organism
name that carries a strain (a token after the binomial containing a digit), form one unit. Bare binomials and
subspecies trinomials are not merged by name (PREREG log). Cached in data/work/units.json.

    from units import unit_of; u = unit_of()   # accession -> unit id
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data/work/units.json"


def _dblinks() -> dict[str, set[str]]:
    keys = defaultdict(set)
    cur, in_db = None, False
    with gzip.open(ROOT / "data/raw/refseq_mito/mitochondrion.1.genomic.gbff.gz", "rt") as fh:
        for line in fh:
            if line.startswith("VERSION"):
                cur = line.split()[1]
            elif line.startswith("DBLINK"):
                in_db = True
            elif in_db and not line.startswith("            "):
                in_db = False
            if in_db and cur:
                for tag in ("BioSample:", "Assembly:"):
                    if tag in line:
                        keys[cur].add(tag + line.split(tag, 1)[1].strip())
    for r in csv.DictReader(open(ROOT / "data/raw/fungal_tpa/All_fungal_mitochondrial_genomes_accessions.tsv"), delimiter="\t"):
        if r.get("BioSample") and r["BioSample"] != "NA":
            keys[r["NCBI_accession_number"] + ".1"].add("BioSample:" + r["BioSample"])
    return keys


def unit_of() -> dict[str, str]:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    G = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    keys = _dblinks()
    parent = {a: a for a in G}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    by_key = defaultdict(list)
    for a, r in G.items():
        toks = r["organism"].split()
        # merge by name only when the name carries a strain: a token after the binomial containing a digit
        if len(toks) > 2 and any(re.search(r"\d", t) for t in toks[2:]):
            by_key["org:" + r["organism"]].append(a)
        for k in keys.get(a, ()):
            by_key[k].append(a)
    for members in by_key.values():
        for b in members[1:]:
            union(members[0], b)
    u = {a: find(a) for a in G}
    CACHE.write_text(json.dumps(u))
    return u


if __name__ == "__main__":
    u = unit_of()
    sizes = defaultdict(int)
    for v in u.values():
        sizes[v] += 1
    multi = {k: v for k, v in sizes.items() if v > 1}
    print(f"{len(u)} records, {len(sizes)} units, {len(multi)} units with > 1 record")
    print([a for a, v in u.items() if v == u.get("NW_009276970.1")])
