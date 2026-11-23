"""Prepare the scan inputs from RefSeq mitochondrion 238 and the 2,695 fungal TPA records.

Outputs (data/work/, not committed):
  genomes.tsv        one row per record: set, taxonomy pools, genetic code, length, CDS count
  cds.tsv            annotated CDS (RefSeq only): normalised gene, strand, parts, protein id
  introns.tsv        annotated intron edges (intron features and CDS-join gaps > 30 nt), RefSeq only
  queries/code_<n>.fna   genome sequences grouped by mitochondrial genetic code (main + land plants)
  ref_proteins.faa   reference proteins: one representative RefSeq genome per family pool
  ref_proteins.tsv   reference protein -> gene, family pool, accession

    uv run python code/prep.py
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
WORK = ROOT / "data/work"

GENE_PATTERNS = [  # (normalised name, regex on gene/product), first match wins
    ("nad4l", r"\bn(a)?d(h)?4l\b|nadh dehydrogenase subunit 4l|subunit 4L"),
    ("cox1", r"\bco(x)?1\b|\bcoi\b|cytochrome c oxidase subunit (1|i)\b"),
    ("cox2", r"\bco(x)?2\b|\bcoii\b|cytochrome c oxidase subunit (2|ii)\b"),
    ("cox3", r"\bco(x)?3\b|\bcoiii\b|cytochrome c oxidase subunit (3|iii)\b"),
    ("cob", r"\bcob\b|\bcyt ?b\b|cytochrome b\b|apocytochrome b"),
    ("atp6", r"\batp(ase)? ?6\b|atp synthase (f0|fo)? ?subunit (6|a)\b"),
    ("atp8", r"\batp(ase)? ?8\b|atp synthase (f0|fo)? ?subunit 8\b"),
    ("atp9", r"\batp(ase)? ?9\b|atp synthase (f0|fo)? ?subunit (9|c)\b"),
    ("nad1", r"\bn(a)?d(h)?1\b|nadh dehydrogenase subunit 1\b"),
    ("nad2", r"\bn(a)?d(h)?2\b|nadh dehydrogenase subunit 2\b"),
    ("nad3", r"\bn(a)?d(h)?3\b|nadh dehydrogenase subunit 3\b"),
    ("nad4", r"\bn(a)?d(h)?4\b|nadh dehydrogenase subunit 4\b"),
    ("nad5", r"\bn(a)?d(h)?5\b|nadh dehydrogenase subunit 5\b"),
    ("nad6", r"\bn(a)?d(h)?6\b|nadh dehydrogenase subunit 6\b"),
]
GENE_RE = [(g, re.compile(p, re.I)) for g, p in GENE_PATTERNS]


def norm_gene(q: dict) -> str:
    text = " | ".join(q.get("gene", []) + q.get("product", []))
    for g, rx in GENE_RE:
        if rx.search(text):
            return g
    return "other:" + (q.get("gene", q.get("product", ["?"]))[0][:30].replace("\t", " "))


class Taxonomy:
    def __init__(self, d: Path):
        self.parent, self.rank, self.mgc = {}, {}, {}
        with open(d / "nodes.dmp") as fh:
            for line in fh:
                f = line.split("\t|\t")
                t = int(f[0])
                self.parent[t], self.rank[t], self.mgc[t] = int(f[1]), f[2], int(f[8])
        self.merged = {}
        with open(d / "merged.dmp") as fh:
            for line in fh:
                a, b = line.rstrip("\t|\n").split("\t|\t")
                self.merged[int(a)] = int(b)
        self.names: dict[int, str] = {}
        self._names_file = d / "names.dmp"

    def fix(self, t: int) -> int:
        return self.merged.get(t, t)

    def lineage(self, t: int) -> list[int]:
        t, out = self.fix(t), []
        while t in self.parent and t != 1:
            out.append(t)
            t = self.parent[t]
        return out  # leaf first

    def load_names(self, wanted: set[int]) -> None:
        with open(self._names_file) as fh:
            for line in fh:
                if "scientific name" not in line:
                    continue
                f = line.split("\t|\t")
                t = int(f[0])
                if t in wanted:
                    self.names[t] = f[1]

    def pools(self, t: int) -> dict:
        lin = self.lineage(t)
        by_rank = {self.rank[x]: x for x in lin}
        # nearest named clade above genus where a rank is missing
        above_genus = lin
        if "genus" in by_rank:
            above_genus = lin[lin.index(by_rank["genus"]) + 1:]
        elif "species" in by_rank:
            above_genus = lin[lin.index(by_rank["species"]) + 1:]
        fam = by_rank.get("family") or (above_genus[0] if above_genus else t)
        fam_i = lin.index(fam) if fam in lin else 0
        order = by_rank.get("order") or (lin[fam_i + 1] if fam_i + 1 < len(lin) else fam)
        return {"family_pool": fam, "family_is_rank": int("family" in by_rank), "order_pool": order,
                "order_is_rank": int("order" in by_rank), "class": by_rank.get("class", 0),
                "phylum": by_rank.get("phylum", 0), "kingdom": by_rank.get("kingdom", 0),
                "superkingdom": by_rank.get("superkingdom", by_rank.get("domain", 0))}


def parts_of(loc) -> list[tuple[int, int, int]]:
    ps = loc.parts if isinstance(loc, CompoundLocation) else [loc]
    return [(int(p.start), int(p.end), p.strand or 1) for p in ps]


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "queries").mkdir(exist_ok=True)
    tax = Taxonomy(RAW / "taxdump")
    print("taxonomy loaded", flush=True)

    # sequences of all RefSeq records come from the .fna (205 gbff records lack sequence)
    seqs = {}
    with gzip.open(RAW / "refseq_mito/mitochondrion.1.1.genomic.fna.gz", "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            seqs[rec.id] = str(rec.seq).upper()

    genomes, cds_rows, intron_rows, prot = [], [], [], {}
    with gzip.open(RAW / "refseq_mito/mitochondrion.1.genomic.gbff.gz", "rt") as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            src = next((f for f in rec.features if f.type == "source"), None)
            taxid = next((int(x.split(":")[1]) for x in (src.qualifiers.get("db_xref", []) if src else [])
                          if x.startswith("taxon:")), 0)
            lineage_names = rec.annotations.get("taxonomy", [])
            codes = Counter()
            n_cds = 0
            for f in rec.features:
                q = f.qualifiers
                if f.type == "CDS":
                    n_cds += 1
                    if "transl_table" in q:
                        codes[int(q["transl_table"][0])] += 1
                    ps = parts_of(f.location)
                    if any(p for p in f.location.parts if getattr(p, "ref", None)):
                        continue
                    gene = norm_gene(q)
                    cds_rows.append([rec.id, gene, ps[0][2], ";".join(f"{a}-{b}" for a, b, _ in ps),
                                     q.get("protein_id", [""])[0], int("ribosomal_slippage" in q)])
                    if "translation" in q:
                        prot[(rec.id, q.get("protein_id", [f"cds{n_cds}"])[0])] = (gene, q["translation"][0])
                    # CDS-join gaps > 30 nt are introns (or bypass elements) for calibration
                    for (a1, b1, s), (a2, b2, _) in zip(ps, ps[1:]):
                        gap = (a2 - b1) if s != -1 else (a1 - b2)
                        if gap > 30:
                            lo, hi = (b1, a2) if s != -1 else (b2, a1)
                            intron_rows.append([rec.id, gene, lo, hi, "cds_join"])
                elif f.type == "intron":
                    intron_rows.append([rec.id, norm_gene(q), int(f.location.start), int(f.location.end), "intron_feature"])
            code = codes.most_common(1)[0][0] if codes else tax.mgc.get(tax.fix(taxid), 0)
            is_plant = "Embryophyta" in lineage_names
            genomes.append({"accession": rec.id, "set": "refseq_landplant" if is_plant else "refseq_main",
                            "taxid": taxid, "organism": rec.annotations.get("organism", ""),
                            "code": code, "code_source": "cds" if codes else "taxonomy",
                            "length": len(seqs.get(rec.id, "")), "n_cds": n_cds,
                            "has_seq": int(rec.id in seqs)})
    print(f"RefSeq parsed: {len(genomes)} records", flush=True)

    with gzip.open(RAW / "fungal_tpa/BK072095-BK074789.gb.gz", "rt") as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            src = next((f for f in rec.features if f.type == "source"), None)
            taxid = next((int(x.split(":")[1]) for x in (src.qualifiers.get("db_xref", []) if src else [])
                          if x.startswith("taxon:")), 0)
            seqs[rec.id] = str(rec.seq).upper()
            genomes.append({"accession": rec.id, "set": "fungal_tpa", "taxid": taxid,
                            "organism": rec.annotations.get("organism", ""),
                            "code": tax.mgc.get(tax.fix(taxid), 0), "code_source": "taxonomy",
                            "length": len(rec.seq), "n_cds": 0, "has_seq": 1})
    print(f"with fungal TPA: {len(genomes)} records", flush=True)

    wanted = set()
    for g in genomes:
        wanted.update(tax.lineage(g["taxid"]))
    tax.load_names(wanted)
    for g in genomes:
        p = tax.pools(g["taxid"])
        g.update(p)
        for k in ("family_pool", "order_pool", "class", "phylum", "kingdom", "superkingdom"):
            g[k + "_name"] = tax.names.get(p[k], "")
    cols = list(genomes[0].keys())
    with open(WORK / "genomes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(genomes)
    with open(WORK / "cds.tsv", "w") as fh:
        fh.write("accession\tgene\tstrand\tparts\tprotein_id\tslippage\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in cds_rows)
    with open(WORK / "introns.tsv", "w") as fh:
        fh.write("accession\tgene\tstart\tend\tsource\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in intron_rows)

    # queries grouped by genetic code
    by_code = defaultdict(list)
    for g in genomes:
        if g["accession"] in seqs and g["code"]:
            by_code[g["code"]].append(g["accession"])
    for code, accs in by_code.items():
        with open(WORK / f"queries/code_{code}.fna", "w") as fh:
            for a in accs:
                fh.write(f">{a}\n{seqs[a]}\n")
    print("genomes per code:", {c: len(a) for c, a in sorted(by_code.items())})

    # reference proteins: one representative RefSeq genome per family pool (most proteins, then longest)
    gi = {g["accession"]: g for g in genomes}
    per_acc = defaultdict(list)
    for (acc, pid), (gene, aa) in prot.items():
        per_acc[acc].append((pid, gene, aa))
    best = {}
    for acc, ps in per_acc.items():
        fam = gi[acc]["family_pool"]
        key = (len(ps), gi[acc]["length"])
        if fam not in best or key > best[fam][0]:
            best[fam] = (key, acc)
    n = 0
    with open(WORK / "ref_proteins.faa", "w") as fa, open(WORK / "ref_proteins.tsv", "w") as ft:
        ft.write("ref_id\taccession\tgene\tfamily_pool\torder_pool\tlength\n")
        for fam, (_, acc) in sorted(best.items()):
            for pid, gene, aa in per_acc[acc]:
                rid = f"{acc}|{pid}"
                fa.write(f">{rid}\n{aa}\n")
                ft.write(f"{rid}\t{acc}\t{gene}\t{fam}\t{gi[acc]['order_pool']}\t{len(aa)}\n")
                n += 1
    print(f"reference proteins: {n} from {len(best)} family representatives")


if __name__ == "__main__":
    main()
