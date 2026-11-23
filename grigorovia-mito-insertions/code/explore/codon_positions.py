"""Codon counts per gene in G. jiainica's seven core mitochondrial genes (extra bases removed), from the codon alignments
used for dN/dS: codons, GGC codons, glycine codons by type, and the base at the third codon position.

    uv run python code/explore/codon_positions.py   -> results/explore/grigorovia/h1h2/genus/codon_positions.md
"""
import collections as C
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ALN = ROOT / "results/explore/grigorovia/h1h2/genus/dnds"   # published as candidate/followup/genus/codon_alignments/
GENES = [("cox1", "cox1.phy"), ("cob", "cob.phy"), ("cox2, cox3, atp6, atp8, atp9", "cox2_cox3_atp6_atp8_atp9.phy")]


def read_phy(p):
    seqs = {}
    for line in [l.rstrip("\n") for l in open(p) if l.strip()][1:]:
        name, s = line.split(None, 1)
        seqs[name] = seqs.get(name, "") + s.replace(" ", "")
    return seqs


def main():
    L = ["# Codons in G. jiainica's seven core mitochondrial genes (extra bases removed; post hoc)", "",
         "From the codon alignments used for dN/dS (`codon_alignments/`; G. jiainica rows; gap codons skipped).", "",
         "| genes | codons | GGC | GGT | GGA | GGG | third position A | C | G | T |", "|---|---|---|---|---|---|---|---|---|---|"]
    tot_c, tot_t, tot_g = 0, C.Counter(), C.Counter()
    for label, f in GENES:
        s = read_phy(ALN / f)
        j = next(v for k, v in s.items() if "jiainica" in k.lower())
        cod = [c for c in (j[i:i + 3].upper() for i in range(0, len(j) - 2, 3)) if "-" not in c]
        t, g = C.Counter(c[2] for c in cod), C.Counter(c for c in cod if c.startswith("GG"))
        tot_c += len(cod); tot_t += t; tot_g += g
        L.append(f"| {label} | {len(cod)} | {g['GGC']} | {g['GGT']} | {g['GGA']} | {g['GGG']} | "
                 + " | ".join(f"{t[b]}" for b in "ACGT") + " |")
    L.append(f"| **all seven** | {tot_c} | {tot_g['GGC']} | {tot_g['GGT']} | {tot_g['GGA']} | {tot_g['GGG']} | "
             + " | ".join(f"{tot_t[b]} ({tot_t[b] / tot_c:.1%})" for b in "ACGT") + " |")
    L += ["", f"Glycine codons other than the {tot_g['GGC']} GGC: {tot_g['GGT'] + tot_g['GGA'] + tot_g['GGG']}, "
          f"none of them ending in C."]
    out = ROOT / "results/explore/grigorovia/h1h2/genus/codon_positions.md"
    out.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
