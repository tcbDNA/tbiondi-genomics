"""Cut the NCBI taxonomy rows the analysis uses from the dump it was run with (2026-11-21, taxdump.tar.gz md5
6d18d6c9d60bbc20759869f6b3592eff), so a rerun reproduces the published values without that 80 MB file.

Rows kept, in the dump's own format:
- nodes.dmp: every node on the lineage (to the root) of each scanned genome's taxid (data/work/genomes.tsv, after
  merged-id replacement) and of each species in the BUSCO saccharomycetes_odb12 reference set (nuclear_genes.py);
- merged.dmp: rows whose old id is one of those query taxids;
- names.dmp: scientific-name rows for every kept node.

    uv run python code/release/make_taxonomy_extract.py OUT_DIR      (e.g. <release>/protocol/taxonomy_extract)
"""
import csv
import gzip
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TAX = ROOT / "data/raw/taxdump"


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    assert hashlib.md5((TAX / "taxdump.tar.gz").read_bytes()).hexdigest() == "6d18d6c9d60bbc20759869f6b3592eff"
    query = {int(r["taxid"]) for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    with gzip.open(ROOT / "data/busco_lineages/saccharomycetes_odb12/refseq_db.faa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                query.add(int(line[1:].split("_")[1]))
    merged, merged_rows = {}, {}
    for line in open(TAX / "merged.dmp"):
        a, b = line.rstrip("\t|\n").split("\t|\t")
        merged[int(a)] = int(b)
        merged_rows[int(a)] = line
    parent, rows = {}, {}
    for line in open(TAX / "nodes.dmp"):
        f = line.split("\t|\t")
        parent[int(f[0])] = int(f[1])
        rows[int(f[0])] = line
    keep = set()
    for t in query:
        t = merged.get(t, t)
        while t in parent and t not in keep:
            keep.add(t)
            if t == 1:
                break
            t = parent[t]
    with open(out / "nodes.dmp", "w") as o:
        for t in sorted(keep):
            o.write(rows[t])
    with open(out / "merged.dmp", "w") as o:
        for t in sorted(query):
            if t in merged_rows:
                o.write(merged_rows[t])
    n = 0
    with open(out / "names.dmp", "w") as o:
        for line in open(TAX / "names.dmp"):
            if "scientific name" in line and int(line.split("\t|\t")[0]) in keep:
                o.write(line)
                n += 1
    (out / "README.md").write_text(
        "# NCBI taxonomy rows used by the analysis\n\n"
        f"Cut by `code/release/make_taxonomy_extract.py` from the NCBI taxonomy dump downloaded on 2026-11-21 "
        f"(`taxdump.tar.gz`, md5 6d18d6c9d60bbc20759869f6b3592eff): {len(keep)} nodes (the lineages of {len(query)} "
        f"query taxa: the scanned genomes and the BUSCO saccharomycetes_odb12 reference species), "
        f"{sum(1 for t in query if t in merged_rows)} merged-id rows and {n} scientific names, in the dump's own "
        "format. `prep.py`, `flag_context.py`, `explore/nuclear_genes.py` and `explore/grigorovia_h1h2.py` read these "
        "files when they are present, so a rerun reproduces the published taxonomy assignments.\n")
    print(f"query taxa {len(query)}, nodes {len(keep)}, names {n}")


if __name__ == "__main__":
    main(sys.argv[1])
