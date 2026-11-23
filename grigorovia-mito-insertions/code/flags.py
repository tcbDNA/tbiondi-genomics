"""Confirmatory test (single genome): codon pile-up of clean shifts, with expectations from planted shifts.

One window rule (DESIGN.md amendment 4, item 1; PREREG.md section 1):
    A clean shift VOTES for the old-frame codon holding its rightmost protein-equivalent placement and the
    one holding its leftmost (one vote if they are the same codon), provided its placement window
    (rightmost - leftmost, nt) is <= the window cap. Shifts with wider windows count in n but do not vote.
For each genome, direction d (+1/-1) and codon c:
    n = clean shifts of direction d in core genes (voting or not); k = those voting for c; g = distinct genes
        (split-hit shifts are included only if they passed the probation rule, PREREG log 07:10)
    q = (votes for c + 0.5) / (N + 32) among N clean recovered planted random shifts of direction d in the
        genome's stratum (code x GC half, backoff code, then all), same pipeline and same rule
    flag if k >= 2, g >= 2 and P(Binomial(n, q) >= k) <= alpha,
    alpha = ALPHA_SENSE for sense codons, ALPHA_STOP for stop codons of the genome's code.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

CORE = {"cox1", "cox2", "cox3", "cob", "atp6", "atp8", "atp9", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6"}
CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
ALPHA_SENSE = 1e-3   # amendment 4, item 3: the pre-written rule's values, fixed
ALPHA_STOP = 3.2e-3
POOL_SITE_RESIDUES = 50   # pool "same site": same gene and direction, subject positions within 50 residues
THRESHOLDS = Path(__file__).resolve().parents[1] / os.environ.get("CTRL_OUT", "results/controls") / "thresholds.json"


def binom_sf(k: int, n: int, q: float) -> float:
    """P(X >= k) for X ~ Binomial(n, q)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    return min(1.0, sum(math.comb(n, i) * q ** i * (1 - q) ** (n - i) for i in range(k, n + 1)))


def votes(e: dict, cap: int) -> set[str]:
    """Codons a shift votes for under the window rule (empty if the window is wider than the cap)."""
    if int(e["window"]) > cap:
        return set()
    return {e["codon"], e["codon_left"]} - {"NNN"}


class Null:
    """Planted-null vote rates by stratum, with backoff and pseudo-counts.

    split_votes: whether split-hit shifts vote (PREREG log 07:10: decided by `bars`, stored in
    thresholds.json). The null is built the same way: with split_votes False, a planted shift counts only if
    the single-hit path recovers it clean."""

    def __init__(self, path: Path, split_votes: bool | None = None):
        d = json.loads(Path(path).read_text())
        if split_votes is None:
            split_votes = bool(json.loads(THRESHOLDS.read_text())["split_votes"])
        self.split_votes = split_votes
        v = d["variants"]["split" if split_votes else "nosplit"]
        self.counts = {k: Counter(x) for k, x in v["counts"].items()}  # key: stratum|dir -> votes per codon
        self.totals = {k: int(x) for k, x in v["totals"].items()}      # key: stratum|dir -> N clean plants
        self.cap = int(v["window_cap"])
        self.gc_cut = d["gc_cut"]  # code -> GC median of its planted genomes
        self.min_events = d.get("min_events", 500)

    def stratum(self, code: int, gc: float) -> str:
        cut = self.gc_cut.get(str(code))
        if cut is None:
            return f"{code}|all"
        return f"{code}|{'lo' if gc < cut else 'hi'}"

    def key(self, code: int, gc: float, d: int) -> str:
        for key in (f"{self.stratum(code, gc)}|{d}", f"{code}|all|{d}", f"ALL|{d}"):
            if self.totals.get(key, 0) >= self.min_events:
                return key
        return f"ALL|{d}"

    def q(self, code: int, gc: float, d: int, codon: str) -> tuple[float, str]:
        key = self.key(code, gc, d)
        return (self.counts.get(key, Counter()).get(codon, 0) + 0.5) / (self.totals.get(key, 0) + 32.0), key


def clean_core(events: list[dict], split_votes: bool = True) -> list[dict]:
    return [e for e in events if e["class"] == "clean" and e["gene"] in CORE and (split_votes or e["via"] != "split")]


def genome_tests(events: list[dict], code: int, gc: float, null: Null, stops: set[str],
                 alpha_sense: float = ALPHA_SENSE, alpha_stop: float = ALPHA_STOP) -> list[dict]:
    """Return one row per (dir, codon) with k >= 2; 'flag' marks those passing."""
    clean = clean_core(events, null.split_votes)
    out = []
    for d in (1, -1):
        sub = [e for e in clean if int(e["shift"]) == d]
        n = len(sub)
        by = defaultdict(list)
        for e in sub:
            for c in votes(e, null.cap):
                by[c].append(e)
        for codon, es in by.items():
            k = len(es)
            if k < 2:
                continue
            genes = {e["gene"] for e in es}
            q, src = null.q(code, gc, d, codon)
            p = binom_sf(k, n, q)
            is_stop = codon in stops
            alpha = alpha_stop if is_stop else alpha_sense
            out.append({"shift": d, "codon": codon, "k": k, "n": n, "genes": len(genes),
                        "q": q, "q_src": src, "p": p, "is_stop": int(is_stop),
                        "flag": int(k >= 2 and len(genes) >= 2 and p <= alpha)})
    return out


def kmin(n: int, q: float, alpha: float) -> int | None:
    """Smallest k >= 2 with P(Binomial(n, q) >= k) <= alpha, or None."""
    for k in range(2, n + 1):
        if binom_sf(k, n, q) <= alpha:
            return k
    return None


def expected_chance(n: int, qs: dict[str, float], alpha_of: dict[str, float]) -> float:
    """Expected number of (codon) results reaching alpha_of[codon] in one genome-direction with n clean
    shifts, if the shifts were random (planted-null rates qs). Ignores the >= 2 genes condition, so it is
    an upper bound."""
    if n < 2:
        return 0.0
    tot = 0.0
    for c, q in qs.items():
        k = kmin(n, q, alpha_of[c])
        if k is not None:
            tot += binom_sf(k, n, q)
    return tot


def genome_dir_inputs(events: list[dict], code: int, gc: float, null: Null, stops: set[str]):
    """(n, qs, alpha_of) per direction, for expected_chance()."""
    clean = clean_core(events, null.split_votes)
    out = []
    for d in (1, -1):
        n = sum(int(e["shift"]) == d for e in clean)
        qs = {c: null.q(code, gc, d, c)[0] for c in CODONS}
        alpha_of = {c: (ALPHA_STOP if c in stops else ALPHA_SENSE) for c in CODONS}
        out.append((n, qs, alpha_of))
    return out
