"""Funnel step "a second assembly carries it", pre-approved by the user only for V. dahliae VdLs.17 and only if it
comes out strong (PREREG log, conditional second-assembly step). Local RefSeq records only; no reads.

    uv run python code/second_assembly_check.py --scan results/scan/main --flagged NW_009276968.1 NW_009276969.1 NW_009276970.1 \
        --others NC_008248.1 NC_029238.1 --out results/scan/analysis/second_assembly_VdLs17.tsv

For every clean shift on the flagged records: align a 121-nt window (+-60) locally to each independent assembly.
Comparable if >= 100 query bases aligned at >= 85% identity. -1 shift (query lacks a base): the other assembly
"carries the base" if the alignment has a gap in the query within 3 nt of the shift. +1 shift: it "carries the
correct form" if the alignment has a gap in the other assembly there. A position counts as carried when the first
other assembly carries it, or, where that one is not comparable, when the next does. Killed if carried at > 50%
of comparable positions.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import Align

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scan as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
COMP = str.maketrans("ACGTN", "TGCAN")


def aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -5, -2
    return a


def check(window: str, shift_at: int, shift: int, other: str, al: Align.PairwiseAligner) -> dict:
    """window: query bases; shift_at: index in window of the shift; other: full other genome (both strands tried)."""
    best = None
    for strand, target in ((1, other + other[:200]), (-1, (other + other[:200]).translate(COMP)[::-1])):  # circular
        aln = al.align(window, target)[0]
        if best is None or aln.score > best[0].score:
            best = (aln, strand)
    aln, strand = best
    # Biopython: align(seqA, seqB) -> aligned[0] = seqA (here the window), aligned[1] = seqB (the other genome)
    qblocks, tblocks = aln.aligned
    other_seq = aln.query
    q_cov = sum(e - s for s, e in qblocks)
    ident = sum(1 for (qs, qe), (ts, te) in zip(qblocks, tblocks) for i in range(qe - qs) if window[qs + i] == other_seq[ts + i])
    ident /= max(1, q_cov)
    if q_cov < 100 or ident < 0.85:
        return {"comparable": 0, "carried": "", "identity": round(ident, 3), "aligned": q_cov}
    gaps_in_query, gaps_in_other = [], []
    for (qs1, qe1), (ts1, te1), (qs2, qe2), (ts2, te2) in zip(qblocks, tblocks, qblocks[1:], tblocks[1:]):
        dq, dt = qs2 - qe1, ts2 - te1
        if dt > dq:      # other assembly has extra bases here: a gap in the query
            gaps_in_query.append((qe1, dt - dq))
        elif dq > dt:    # query has extra bases: a gap in the other assembly
            gaps_in_other.append((qe1, dq - dt))
    near = (lambda gl: any(abs(p - shift_at) <= 3 and n % 3 != 0 for p, n in gl))
    carried = near(gaps_in_query) if shift == -1 else near(gaps_in_other)
    return {"comparable": 1, "carried": int(carried), "identity": round(ident, 3), "aligned": q_cov}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", nargs="+", required=True)
    ap.add_argument("--flagged", nargs="+", required=True)
    ap.add_argument("--others", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    G = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    need = set(args.flagged) | set(args.others)
    seqs = {}
    for code in {G[a]["code"] for a in need}:
        seqs.update(S.read_fasta(ROOT / f"data/work/queries/code_{code}.fna", need))
    ev = []
    for d in args.scan:
        ev += [e for e in csv.DictReader(open(Path(d) / "events.tsv"), delimiter="\t") if e["accession"] in args.flagged]
    al = aligner()
    rows = []
    for e in ev:
        if e["class"] != "clean":
            continue
        s = seqs[e["accession"]]
        p = int(e["pos_fwd"])
        lo, hi = max(0, p - 60), min(len(s), p + 61)
        win = s[lo:hi]
        shift = int(e["shift"])
        # work on the coding strand so that +1/-1 keep their meaning; gaps are strand-symmetric for this test
        r = {"accession": e["accession"], "gene": e["gene"], "pos_fwd": p, "shift": shift}
        carried_final = ""
        for o in args.others:
            c = check(win, p - lo, shift, seqs[o], al)
            r.update({f"{o}_comparable": c["comparable"], f"{o}_carried": c["carried"], f"{o}_identity": c["identity"]})
            if carried_final == "" and c["comparable"]:
                carried_final = c["carried"]
        r["carried"] = carried_final
        rows.append(r)
    comp = [r for r in rows if r["carried"] != ""]
    n_car = sum(r["carried"] == 1 for r in comp)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["accession"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    verdict = ("strong, killed at the second-assembly step: assembly error" if comp and n_car > 0.5 * len(comp)
               else "not killed: funnel waits for the user")
    print(f"{len(rows)} clean shift positions; comparable {len(comp)}; carried by an independent assembly {n_car}; "
          f"verdict: {verdict}")


if __name__ == "__main__":
    main()
