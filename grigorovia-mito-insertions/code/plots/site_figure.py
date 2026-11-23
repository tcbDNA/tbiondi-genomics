"""Static figure: the Grigorovia jiainica cob double site (plots/grigorovia_cob_site.png).

    uv run python code/plots/site_figure.py [--data candidate/funnel/figure_site.json] [--out plots]

The site data (codons, their reading in the yeast mitochondrial code, the S. cerevisiae COB residue aligned to each
codon, read counts) are written by `code/explore/grigorovia_followup.py figure` to figure_site.json. Fonts are
matplotlib's bundled DejaVu, so the script runs anywhere.
"""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
       "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
       "Y": "Tyr", "V": "Val", "*": "stop"}
PEACH, ORANGE, INK, GREY, LIGHT = "#f8dcc6", "#e8650c", "#1d1d1d", "#8a8a85", "#b9b9b4"
MONO, SANS = "DejaVu Sans Mono", "DejaVu Sans"


def window(d: dict, first: int = 3, drop_last: int = 0) -> list[dict]:
    items = d["items"][first:]
    return items[:len(items) - drop_last] if drop_last else items


def site_figure(d: dict, out: Path) -> None:
    items = window(d, 3, 2)
    W, H = 16, 9
    fig = plt.figure(figsize=(W, H), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    ax.text(0.55, 8.15, "Two extra bases in cob (DNA only)", fontsize=30, weight="bold", family=SANS, color=INK)
    ax.text(0.55, 7.5, f"Mitochondrial cytochrome b gene of the yeast {d['species']}: removing the two bases restores the "
            f"residues of other yeasts", fontsize=16, family=SANS, color=GREY)
    x0, box_w, skip_w, gap = 3.3, 0.86, 0.36, 0.09
    y_dna, y_pro, y_ref = 5.55, 4.35, 2.85
    for lbl, y in (("DNA", y_dna), ("Protein", y_pro)):
        ax.text(0.55, y, lbl, fontsize=19, weight="bold", family=SANS, va="center", color=INK)
    ax.text(0.55, y_ref, "Same protein in\nS. cerevisiae", fontsize=17, weight="bold", family=SANS, va="center",
            color=INK, linespacing=1.3)
    x = x0
    n_same = n_ref = 0
    after_skip = False
    for it in items:
        if it["kind"] == "skip":
            ax.add_patch(Rectangle((x, y_dna - 0.36), skip_w, 0.72, fill=False, ls=(0, (3, 2)), ec=LIGHT, lw=1.6))
            ax.text(x + skip_w / 2, y_dna, it["base"], fontsize=20, family=MONO, ha="center", va="center", color=LIGHT)
            ax.text(x + skip_w / 2, y_dna + 1.0, "extra\nbase", fontsize=12, family=SANS, ha="center", va="center", color=INK)
            ax.annotate("", xy=(x + skip_w + 0.5, y_dna + 0.6), xytext=(x - 0.05, y_dna + 0.6),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
            x += skip_w + gap
            after_skip = True
            continue
        hi = after_skip and it["aa"] == "G"
        after_skip = False
        ax.add_patch(Rectangle((x, y_dna - 0.36), box_w, 0.72, color=ORANGE if hi else PEACH, ec="none"))
        ax.text(x + box_w / 2, y_dna, it["codon"], fontsize=20, family=MONO, ha="center", va="center",
                color="white" if hi else INK, weight="bold" if hi else "normal")
        col = ORANGE if hi else INK
        ax.plot([x + box_w / 2] * 2, [y_dna - 0.5, y_pro + 0.3], color=ORANGE if hi else LIGHT, lw=1.2)
        ax.text(x + box_w / 2, y_pro, AA3[it["aa"]], fontsize=19, family=SANS, ha="center", va="center", color=col,
                weight="bold" if hi else "normal")
        if it.get("ref"):
            n_ref += 1
            same = it["ref"] == it["aa"]
            n_same += same
            if same:
                ax.plot([x + box_w / 2] * 2, [y_pro - 0.3, y_ref + 0.3], color=ORANGE if hi else INK, lw=1.6)
            ax.text(x + box_w / 2, y_ref, AA3[it["ref"]], fontsize=19, family=SANS, ha="center", va="center",
                    color=(ORANGE if hi else INK) if same else LIGHT, weight="bold" if hi else "normal")
        x += box_w + gap
    wo = d["without_skip"]
    ax.text(0.55, 1.6, "Read straight on:", fontsize=17, weight="bold", family=SANS, color=INK, va="center")
    ax.text(3.3, 1.6, "  ".join(c for c, _ in wo), fontsize=16, family=MONO, color=INK, va="center")
    ax.text(3.3, 1.12, "= " + " ".join(AA3[a] for _, a in wo) + ": the protein stops short", fontsize=15, family=SANS,
            color=INK, va="center")
    r = d["reads_with_without"]
    ax.text(0.55, 0.55,
            f"Genome {d['genome']} (strain {d['strain']}; Sarhan et al. 2025), yeast mitochondrial code. Bars: same amino acid in "
            f"S. cerevisiae COB ({n_same} of {n_ref}).\nIllumina reads behind the assembly carry both extra bases "
            f"({r['26496'][0]:,} and {r['26477'][0]:,} reads; {r['26496'][1] + r['26477'][1]} without). Cause unknown: a translational "
            f"skip, RNA editing, or disabling mutations.",
            fontsize=12.5, family=SANS, color=GREY, va="center", linespacing=1.5)
    fig.savefig(out, dpi=100)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="candidate/funnel/figure_site.json")
    ap.add_argument("--out", default="plots")
    args = ap.parse_args()
    d = json.load(open(args.data))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    site_figure(d, out / "grigorovia_cob_site.png")
    print("wrote", out / "grigorovia_cob_site.png")


if __name__ == "__main__":
    main()
