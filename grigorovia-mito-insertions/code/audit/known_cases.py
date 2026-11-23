"""Prior-art audit: every frameshift already written into RefSeq mitochondrial CDS joins, with codon context.

A CDS written as join(a..b, c..d) on one strand with c - b - 1 in {-2, -1, +1, +2} is read as a
frameshift (+1 = one base skipped, -1 = one base re-read, +2 = two skipped, -2 = two re-read).
Gap 0 (parts touching) is reported separately and not counted as a shift. Annotation only; no scan.

For each junction it reports the genome coordinate, the skipped/re-read bases, and the old-frame
context: the codon holding the first skipped base, plus one codon either side, in the upstream frame.

    uv run python code/audit/known_cases.py data/raw/refseq_mito/mitochondrion.1.genomic.gbff.gz refseq238
Writes results/audit/known_junctions_<label>.tsv
"""

from __future__ import annotations

import gzip
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    path, label = Path(sys.argv[1]), sys.argv[2]
    rows, gap0, no_seq = [], 0, []
    with gzip.open(path, "rt") as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            tax = rec.annotations.get("taxonomy", [])
            group = "land_plants" if "Embryophyta" in tax else next(
                (g for g in ("Metazoa", "Fungi", "Viridiplantae", "Rhodophyta", "Sar", "Discoba", "Amoebozoa") if g in tax), "other")
            fam = next((t for t in reversed(tax) if t.endswith("idae") or t.endswith("aceae")), "")
            n = len(rec.seq)
            try:
                bytes(rec.seq[:1])
            except Exception:  # noqa: BLE001 - record without sequence (e.g. a CON record)
                no_seq.append(rec.id)
                continue
            for f in rec.features:
                if f.type != "CDS" or not isinstance(f.location, CompoundLocation):
                    continue
                parts = list(f.location.parts)
                if any(p.ref for p in parts) or len({p.strand for p in parts}) != 1:
                    continue
                strand = parts[0].strand
                seq = rec.seq if strand != -1 else rec.seq.reverse_complement()
                # transcript-order coordinates on `seq`
                tparts = [(int(p.start), int(p.end)) if strand != -1 else (n - int(p.end), n - int(p.start)) for p in parts]
                cs = int(f.qualifiers.get("codon_start", ["1"])[0]) - 1
                q = f.qualifiers
                gene = q.get("gene", q.get("product", ["?"]))[0]
                coding = -cs
                for (s1, e1), (s2, e2) in zip(tparts, tparts[1:]):
                    coding += e1 - s1
                    gap = s2 - e1
                    if gap == 0:
                        gap0 += 1
                        continue
                    if gap not in (-2, -1, 1, 2):
                        continue
                    phase = coding % 3  # bases of the current old-frame codon already used before e1
                    cstart = e1 - phase  # old-frame codon that holds the first skipped / re-read base
                    ctx = seq[cstart - 3: cstart + 6]
                    gpos = (e1 + 1) if strand != -1 else (n - e1)  # 1-based genome coordinate of first skipped base
                    rows.append({
                        "accession": rec.id, "organism": rec.annotations.get("organism", ""), "group": group,
                        "family": fam, "gene": gene, "strand": strand, "genome_pos": gpos,
                        "cds_pos": coding + cs + 1, "shift": f"+{gap}" if gap > 0 else str(gap),
                        "bases": str(seq[e1:e1 + gap]) if gap > 0 else str(seq[s2:e1]),
                        "old_frame_codons": f"{ctx[0:3]} [{ctx[3:6]}] {ctx[6:9]}",
                        "slippage": int("ribosomal_slippage" in q),
                        "note": "; ".join(q.get("note", []))[:120],
                    })
    out = ROOT / f"results/audit/known_junctions_{label}.tsv"
    cols = list(rows[0])
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(rows)} shift junctions in {len({r['accession'] for r in rows})} records "
          f"(+ {gap0} touching joins with gap 0, not counted; {len(no_seq)} records without sequence skipped: {no_seq[:5]}) -> {out.relative_to(ROOT)}")
    by = Counter((r["group"], r["gene"].upper().replace("NADH DEHYDROGENASE SUBUNIT 3", "ND3").replace("NAD3", "ND3"), r["shift"]) for r in rows)
    for k, v in by.most_common(40):
        print(f"  {v:5d}  {k}")
    print("records by group:", Counter(r["group"] for r in {x["accession"]: x for x in rows}.values()))
    cod = Counter(r["old_frame_codons"].split("[")[1][:3] for r in rows if r["group"] == "Metazoa" and "3" in r["gene"])
    print("Metazoa ND3 junctions, old-frame codon holding the skipped base:", cod.most_common(8))


if __name__ == "__main__":
    main()
