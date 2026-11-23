"""Parse miniprot GFF (with ##PAF lines) for the nuclear gene panel and call present / partial / broken / absent.
Calls (logged in PREREG before running): present = best alignment covers >= 70% of the S. cerevisiae protein with no
frameshift or in-frame stop; broken = a frameshift or in-frame stop within the aligned region; absent = no alignment
covering >= 30%; otherwise partial.

    uv run python code/explore/nuclear_calls.py results/explore/grigorovia/h1h2/nuclear/*.gff
"""
import sys
from pathlib import Path


def parse(path):
    best = {}
    for line in open(path):
        if not line.startswith("##PAF"):
            continue
        f = line.rstrip("\n").split("\t")[1:]
        q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        tags = {t.split(":")[0]: t.split(":")[-1] for t in f[12:] if t.count(":") >= 2}
        rec = {"qlen": qlen, "cov": (qe - qs) / qlen, "fs": int(tags.get("fs", 0)), "st": int(tags.get("st", 0)),
               "contig": f[5], "start": int(f[7]), "end": int(f[8]), "strand": f[4], "score": int(tags.get("ms", 0)),
               "ident": int(f[9]) / max(1, int(f[10])) if f[9].isdigit() else 0.0, "rank": tags.get("cg", "")}
        best.setdefault(q, []).append(rec)
    return best


def call(recs):
    if not recs:
        return "absent", None
    r = max(recs, key=lambda x: x["score"])
    if r["cov"] < 0.30:
        return "absent", r
    if r["fs"] or r["st"]:
        return "broken", r
    if r["cov"] >= 0.70:
        return "present", r
    return "partial", r


if __name__ == "__main__":
    for p in sys.argv[1:]:
        b = parse(p)
        print(Path(p).stem, len(b), "queries aligned")


# ------------------------------------------------------------------ congener-protein queries
def longest_orf_over(seq, a, b, strand):
    """Longest stop-to-stop open frame (standard code) overlapping [a, b) on the given strand; returns protein."""
    from Bio.Seq import Seq
    s = seq if strand == "+" else str(Seq(seq).reverse_complement())
    if strand == "-":
        a, b = len(seq) - b, len(seq) - a
    best = ""
    for f in range(3):
        prot = str(Seq(s[f:len(s) - (len(s) - f) % 3]).translate())
        pos = 0
        for piece in prot.split("*"):
            st, en = f + 3 * pos, f + 3 * (pos + len(piece))
            if en > a and st < b:
                m = piece.find("M")
                cand = piece[m:] if m >= 0 and (f + 3 * (pos + m)) < b else piece
                if len(cand) > len(best):
                    best = cand
            pos += len(piece) + 1
    return best
