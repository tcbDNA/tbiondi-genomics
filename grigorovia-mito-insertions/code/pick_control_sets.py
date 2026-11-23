"""Fix the control genome sets in advance (written to protocol/controls/, committed before any control run).

- null_genomes.txt: >= 100 RefSeq genomes with no known shift, up to 8 per phylum, seed 7. Excluded:
  any record with a frameshift-type annotation (results/audit/annotation_summary.txt definition), and
  every lineage in results/audit/controls_known_cases.tsv (named below), plus all birds and turtles.
- plant_genomes.txt: genomes for the planted null, stratified by genetic code (seed 11), disjoint from the
  null set and from the known lineages, with annotated CDS (positions are planted into annotated CDS).
- intron_genomes.txt: RefSeq fungal genomes with annotated CDS introns (all of them, except NC_079576).

    uv run python code/pick_control_sets.py
"""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "protocol/controls"

KNOWN_TAXA = {  # family_pool or order_pool names excluded from null/planted sets (known shift lineages)
    "Aphrocallistidae", "Euplectellidae", "Rossellidae", "Caulophacidae", "Hyalonematidae", "Sceptrulophora",
    "Lyssacinosida", "Amphidiscosida", "Perkinsidae", "Perkinsida", "Saccharomycetaceae", "Dipodascaceae",
    "Testudines", "Veneridae", "Arcticidae", "Ostreidae", "Anomiidae", "Placunidae", "Plicatulidae",
    "Formicidae", "Phaeodactylaceae", "Rhizophydiaceae", "Spizellomycetaceae", "Physaraceae",
}
EXCLUDE_CLASS_NAMES = {"Aves", "Hexactinellida"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    g = list(csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t"))
    ann = {r["accession"]: r for r in csv.DictReader(open(ROOT / "results/audit/annotation_refseq238.tsv"), delimiter="\t")}
    junction_accs = {r["accession"] for r in csv.DictReader(open(ROOT / "results/audit/known_junctions_refseq238.tsv"), delimiter="\t")}

    def known(r: dict) -> bool:
        a = ann.get(r["accession"], {})
        annotated = (r["accession"] in junction_accs or int(a.get("slippage", 0) or 0) > 0
                     or int(a.get("note_fs", 0) or 0) > 0 or int(a.get("tiny_introns", 0) or 0) > 0)
        taxa = {r["family_pool_name"], r["order_pool_name"]}
        lineage_hit = bool(taxa & KNOWN_TAXA) or r["class_name"] in EXCLUDE_CLASS_NAMES or r["order_pool_name"] == "Testudines"
        return annotated or lineage_hit

    main_ok = [r for r in g if r["set"] == "refseq_main" and int(r["length"]) > 5000 and int(r["n_cds"]) >= 8
               and int(r["code"]) > 0 and not known(r)]
    rng = random.Random(7)
    by_phylum = defaultdict(list)
    for r in main_ok:
        by_phylum[r["phylum_name"] or r["kingdom_name"] or "unranked"].append(r)
    null = []
    for ph in sorted(by_phylum):
        rs = by_phylum[ph]
        rng.shuffle(rs)
        null += rs[:8]
    rng.shuffle(null)
    null = sorted(null[:max(120, min(len(null), 120))], key=lambda r: r["accession"])
    null_ids = {r["accession"] for r in null}

    rng2 = random.Random(11)
    quota = {"2": 120, "5": 120, "4": 160, "3": 60, "9": 60, "1": 60}
    by_code = defaultdict(list)
    for r in main_ok:
        if r["accession"] not in null_ids:
            by_code[r["code"]].append(r)
    plant = []
    for code, rs in sorted(by_code.items()):
        rng2.shuffle(rs)
        plant += rs[:quota.get(code, len(rs))]

    intron_accs = {r["accession"] for r in csv.DictReader(open(ROOT / "data/work/introns.tsv"), delimiter="\t")}
    # NC_079576 ("Pholiota spumosa") removed 2026-11-22: mislabeled ascomycete (DESIGN.md amendment 4)
    intron = [r for r in g if r["set"] == "refseq_main" and r["kingdom_name"] == "Fungi" and r["accession"] in intron_accs
              and r["accession"] not in {"NC_079576.1"}]

    for name, rows in (("null_genomes", null), ("plant_genomes", plant), ("intron_genomes", intron)):
        with open(OUT / f"{name}.txt", "w") as fh:
            fh.write("\n".join(r["accession"] for r in rows) + "\n")
        with open(OUT / f"{name}.tsv", "w") as fh:
            fh.write("accession\tcode\tphylum\tfamily\torganism\n")
            fh.writelines(f"{r['accession']}\t{r['code']}\t{r['phylum_name']}\t{r['family_pool_name']}\t{r['organism']}\n" for r in rows)
        print(f"{name}: {len(rows)} genomes; codes {dict(sorted(defaultdict(int, {c: sum(1 for x in rows if x['code']==c) for c in {x['code'] for x in rows}}).items()))}; phyla {len({r['phylum_name'] for r in rows})}")


if __name__ == "__main__":
    main()
