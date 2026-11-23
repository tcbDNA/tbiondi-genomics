"""Core of the frameshift scan, shared by the controls, the planted null and the real scan.

Everything that decides what counts as a shift and which codon it is assigned to lives here, so that
real, planted and null shifts go through identical code.

DIAMOND frameshift notation (checked on planted indels, 2026-11-21):
  '\\' in the query column of BTOP = an extra base in the DNA; the base at the cursor is skipped (+1).
  '/'  = a missing base; the cursor steps back one base, which is read twice (-1).
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from Bio.Data import CodonTable

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "data/work/ref_final"   # built by code/refs.py from prep.py's ref_proteins
FIELDS = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "slen", "qlen", "qframe", "bitscore",
          "evalue", "btop"]

# ---- parameters fixed in PREREG section 1 (values here are the calibrated ones) -------------------
P = {
    "F": 15,                  # DIAMOND frameshift penalty
    "k": 10,                  # targets per query range before own-family exclusion
    "evalue": 1e-5,
    "sensitivity": "--sensitive",
    "flank_cols": 15,         # aligned columns required on each side of a shift
    "flank_identity": 0.20,   # minimum identity in each flank
    "pair_window": 150,       # opposite-direction shift within this many nt of the same HSP = paired
    "homopolymer_run": 5,     # run length that kills a shift ...
    "homopolymer_reach": 6,   # ... if it reaches within this many nt of the placement window
    "gene_start_aa": 3,       # shifts within the first 3 subject residues are gene-start cases
    "gene_end_aa": 10,        # and within the last 10 are gene-end cases
    "insert_min": 30,         # unaligned insert > 30 nt between HSPs of one subject = intron/insert
    "mask_nt": 30,            # shifts within 30 nt of an insert edge are masked
    "bypass_max": 100,        # inserts of 20..100 nt are bypass-like (tally only)
    "bypass_min": 20,
    "norm_reach": 9,          # placement window searched +-9 nt
    # frame change between two consecutive hits of the same reference protein (amendment 4, item 4)
    "split_overlap_aa": 3,    # protein coordinates contiguous: overlap of at most 3 residues ...
    "split_missing_aa": 1,    # ... or at most 1 residue unaligned between the two hits
    "split_overlap_nt": 10,   # the hits abut on the genome or overlap by at most 10 nt
}
# ----------------------------------------------------------------------------------------------------


def codon_table(code: int) -> dict[str, str]:
    t = CodonTable.unambiguous_dna_by_id[code]
    out = dict(t.forward_table)
    for s in t.stop_codons:
        out[s] = "*"
    return out


def run_diamond(fasta: Path, code: int, out: Path, threads: int = 12, extra: list[str] | None = None) -> None:
    cmd = ["diamond", "blastx", "-d", str(DB), "-q", str(fasta), "--query-gencode", str(code),
           "-F", str(P["F"]), "--range-culling", "-k", str(P["k"]), "-e", str(P["evalue"]),
           "--threads", str(threads), "--quiet", "--outfmt", "6", *FIELDS]
    if P["sensitivity"]:
        cmd.insert(2, P["sensitivity"])
    if extra:
        cmd += extra
    with open(out, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


@dataclass
class Hit:
    q: str
    s: str
    qstart: int
    qend: int
    sstart: int
    send: int
    slen: int
    qlen: int
    qframe: int
    bits: float
    btop: str
    # oriented coordinates: 0-based start/end on the strand of the hit
    ostart: int = 0
    oend: int = 0

    def __post_init__(self):
        if self.qframe > 0:
            self.ostart, self.oend = self.qstart - 1, self.qend
        else:
            self.ostart, self.oend = self.qlen - self.qstart, self.qlen - self.qend + 1

    @property
    def strand(self) -> int:
        return 1 if self.qframe > 0 else -1

    def fwd_interval(self) -> tuple[int, int]:
        return (self.qstart - 1, self.qend) if self.qframe > 0 else (self.qend - 1, self.qstart)


def read_hits(path: Path) -> list[Hit]:
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append(Hit(f[0], f[1], int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]),
                           int(f[8]), float(f[9]), f[11]))
    return out


def cull(hits: list[Hit], own_family: str, ref_family: dict[str, str]) -> list[Hit]:
    """Drop hits to the query's own family, then keep the best hit per genomic locus
    (greedy by bitscore; a hit is kept if it overlaps kept hits on the same strand by < 50% of its length)."""
    keep: list[Hit] = []
    for h in sorted((h for h in hits if ref_family.get(h.s) != own_family), key=lambda h: -h.bits):
        a, b = h.fwd_interval()
        ok = True
        for k in keep:
            if k.strand != h.strand:
                continue
            c, d = k.fwd_interval()
            ov = min(b, d) - max(a, c)
            if ov > 0.5 * (b - a):
                ok = False
                break
        if ok:
            keep.append(h)
    return keep


TOKEN = re.compile(r"\d+|..", re.S)


@dataclass
class Event:
    q: str
    s: str
    strand: int
    pos: int            # oriented 0-based: +1 -> skipped base; -1 -> cursor before stepping back
    shift: int          # +1 or -1 (net, after merging adjacent ops; +2 -> -1, -2 -> +1)
    raw: int            # raw net offset before folding (+1, +2, -1, -2)
    spos: int           # subject residue index (1-based) at the shift
    slen: int
    left_cols: list = field(default_factory=list)   # identity flags of aligned columns before
    right_cols: list = field(default_factory=list)  # ... and after
    left_gap: bool = False
    right_gap: bool = False
    via: str = "hsp"    # "hsp": inside one DIAMOND frameshift alignment; "split": between two hits


def walk(h: Hit, seq_o: str) -> tuple[list[Event], list[tuple[int, str]], list[tuple[int, int]], list]:
    """Walk the BTOP of one hit. Returns shift events, in-frame stop columns (oriented pos, subject aa),
    query insertions inside the HSP (runs of subject-gap columns, merged across adjacent frameshift
    ops) as oriented (start, end), and the column list (type, identical, cursor_after, spos_after)."""
    cols: list[tuple[str, int, int, int]] = []  # (type, identical, cursor_after, spos_after)
    cur, sp = h.ostart, h.sstart
    events: list[Event] = []
    stops: list[tuple[int, str]] = []
    qins: list[list[int]] = []
    in_ins = False
    pending = 0
    pend_pos = pend_sp = None
    for tok in TOKEN.findall(h.btop):
        is_qins = (not tok.isdigit()) and tok[1] == "-" and tok[0] not in "-\\/"
        is_fs = (not tok.isdigit()) and tok[0] in "\\/"
        if is_qins:
            if in_ins:
                qins[-1][1] = cur + 3
            else:
                qins.append([cur, cur + 3])
                in_ins = True
        elif not is_fs:
            in_ins = False
        if tok.isdigit():
            n = int(tok)
            if pending:
                events.append(Event(h.q, h.s, h.strand, pend_pos, 0, pending, pend_sp, h.slen))
                cols.append(("fs", 0, cur, sp))
                pending = 0
            for _ in range(n):
                cur += 3
                sp += 1
                cols.append(("aa", 1, cur, sp))
            continue
        qa, sa = tok[0], tok[1]
        if qa in "\\/":
            if not pending:
                pend_pos, pend_sp = cur, sp
            pending += 1 if qa == "\\" else -1
            cur += 1 if qa == "\\" else -1
            continue
        if pending:
            events.append(Event(h.q, h.s, h.strand, pend_pos, 0, pending, pend_sp, h.slen))
            cols.append(("fs", 0, cur, sp))
            pending = 0
        if qa == "-":
            sp += 1
            cols.append(("gap", 0, cur, sp))
        elif sa == "-":
            cur += 3
            cols.append(("gap", 0, cur, sp))
        else:
            if qa == "*":
                stops.append((cur, sa))
            cur += 3
            sp += 1
            cols.append(("aa", int(qa == sa), cur, sp))
    if pending:
        events.append(Event(h.q, h.s, h.strand, pend_pos, 0, pending, pend_sp, h.slen))
        cols.append(("fs", 0, cur, sp))
    # fold raw offsets: +1/-2 -> +1 ; -1/+2 -> -1 ; multiples of 3 dropped
    out = []
    for e in events:
        r = e.raw % 3
        if r == 0:
            continue
        e.shift = 1 if r == 1 else -1
        # flanks: aligned columns before and after this event, skipping other fs columns
        idx = next(i for i, c in enumerate(cols) if c[0] == "fs" and c[3] == e.spos and abs(c[2] - e.pos) <= 3)
        left, right, lgap, rgap = [], [], False, False
        for c in reversed(cols[:idx]):
            if len(left) >= P["flank_cols"]:
                break
            if c[0] == "gap":
                lgap = True
                break
            if c[0] == "aa":
                left.append(c[1])
        for c in cols[idx + 1:]:
            if len(right) >= P["flank_cols"]:
                break
            if c[0] == "gap":
                rgap = True
                break
            if c[0] == "aa":
                right.append(c[1])
        e.left_cols, e.right_cols, e.left_gap, e.right_gap = left, right, lgap, rgap
        out.append(e)
    return out, stops, [(a, b) for a, b in qins], cols


def _edge_cols(cols: list, from_end: bool) -> tuple[list, bool]:
    """Identity flags of up to flank_cols aligned columns at one end of a hit, stopping at a gap column
    and skipping frameshift columns (as walk() does for in-HSP flanks)."""
    out, gap = [], False
    for c in (reversed(cols) if from_end else cols):
        if len(out) >= P["flank_cols"]:
            break
        if c[0] == "gap":
            gap = True
            break
        if c[0] == "aa":
            out.append(c[1])
    return out, gap


def split_events(kept: list[Hit], cols_of: dict[int, list]) -> list[tuple[int, Event]]:
    """Frame changes between two consecutive hits (genome order, same strand) of the SAME reference protein.
    With a = upstream hit, b = downstream hit (oriented coordinates):
      d_s = b.sstart - a.send - 1   residues between the hits on the protein; must be in
                                    [-split_overlap_aa, +split_missing_aa] (contiguous protein coordinates)
      d_q = b.ostart - a.oend       bases between the hits on the genome; must be >= -split_overlap_nt
      r   = d_q - 3 * d_s           frame residual; must be +-1 or +-2 (folded: +1/-2 -> +1, -1/+2 -> -1)
    so the genome gap is at most what the unaligned residues and one shift explain (the hits abut).
    Wider gaps are left to inserts() (introns, bypass-like inserts). The shift is placed at the middle of
    the junction; trigger() then finds the protein-equivalent placements around it, as for in-HSP shifts."""
    out = []
    by = {}
    for i, h in enumerate(kept):
        by.setdefault((h.strand, h.s), []).append(i)
    for (strand, subj), idx in by.items():
        idx.sort(key=lambda i: kept[i].ostart)
        for i, j in zip(idx, idx[1:]):
            a, b = kept[i], kept[j]
            d_s = b.sstart - a.send - 1
            d_q = b.ostart - a.oend
            r = d_q - 3 * d_s
            if not (-P["split_overlap_aa"] <= d_s <= P["split_missing_aa"]) or d_q < -P["split_overlap_nt"]:
                continue
            if r not in (-2, -1, 1, 2):
                continue
            shift = 1 if r % 3 == 1 else -1
            jn = (a.oend + b.ostart) // 2            # base at the junction
            pos = jn if shift == 1 else jn + 1       # same convention as walk(): +1 skipped base; -1 cursor
            e = Event(a.q, subj, strand, pos, shift, r, (a.send + 1 + b.sstart) // 2, a.slen, via="split")
            e.left_cols, e.left_gap = _edge_cols(cols_of[i], True)
            e.right_cols, e.right_gap = _edge_cols(cols_of[j], False)
            out.append((i, e))   # i = index of the upstream hit in kept
    return out


def translate_skip(seq: str, start: int, t: int, stop: int, table: dict, mode: int) -> str:
    """Translate seq[start:stop] in the frame starting at `start`, with base t skipped (mode +1)
    or re-read (mode -1)."""
    if mode == 1:
        s = seq[start:t] + seq[t + 1:stop]
    else:
        s = seq[start:t + 1] + seq[t:stop]
    return "".join(table.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def trigger(e: Event, seq_o: str, table: dict, frame0: int | None = None) -> dict:
    """Shift-codon convention (one rule for real, planted and null shifts):
    find every placement within +-norm_reach nt that yields the same protein over the local window as
    DIAMOND's placement, take the RIGHTMOST such placement, and report the old-frame codon containing it."""
    R = P["norm_reach"]
    if e.shift == 1:
        base = e.pos            # skipped base
    else:
        base = e.pos - 1        # base read twice
    if frame0 is None:
        frame0 = e.pos % 3      # old-frame codons start at positions == e.pos (mod 3) for DIAMOND events
    # (frame0 is given explicitly only when computing a site's own codons from its annotated/planted position)
    lo = max(frame0, ((base - 24 - frame0) // 3) * 3 + frame0)
    hi = min(len(seq_o), base + 27)
    if lo < 0 or hi - lo < 30:
        return {"codon": "NNN", "codon_left": "NNN", "window": 0, "tpos": base}
    ref = translate_skip(seq_o, lo, base, hi, table, e.shift)
    eq = [t for t in range(max(lo + 1, base - R), min(hi - 4, base + R) + 1)
          if translate_skip(seq_o, lo, t, hi, table, e.shift)[:-1] == ref[:-1]]
    if base not in eq:
        eq.append(base)
    right, left = max(eq), min(eq)

    def codon_at(t: int) -> str:
        cs = t - ((t - frame0) % 3)
        return seq_o[cs:cs + 3]

    return {"codon": codon_at(right), "codon_left": codon_at(left), "window": right - left, "tpos": right}


def old_frame0(parts_fwd: list[tuple[int, int]], strand: int, L: int, opos: int) -> int:
    """Old-frame phase (codon starts == frame0 mod 3, oriented coordinates) at oriented position opos of a CDS
    given as forward parts; the phase of the part holding (or last before) opos, continued from earlier parts."""
    parts = sorted(parts_fwd) if strand == 1 else sorted((L - y, L - x) for x, y in parts_fwd)
    prev = 0
    k = 0
    for i, (a, b) in enumerate(parts):
        if a <= opos:
            k = i
    for a, b in parts[:k]:
        prev += b - a
    return (parts[k][0] + (3 - prev % 3) % 3) % 3


def homopolymer(seq_o: str, a: int, b: int, shift: int = 1) -> tuple[bool, str]:
    """True if a run of identical bases reaches within homopolymer_reach of [a, b]. The run length that
    counts is homopolymer_run for +1 shifts and one less for -1 shifts, because a base deleted from a run
    of 5 leaves a run of 4 in the sequence we see."""
    run = P["homopolymer_run"] - (1 if shift == -1 else 0)
    lo, hi = max(0, a - P["homopolymer_reach"] - run), min(len(seq_o), b + P["homopolymer_reach"] + run)
    region = seq_o[lo:hi]
    for m in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % ((run,) * 4), region):
        s, t = lo + m.start(), lo + m.end()
        if s <= b + P["homopolymer_reach"] and t >= a - P["homopolymer_reach"]:
            return True, m.group(0)[0]
    return False, ""


def subject_to_query(h: Hit, s_target: int) -> int:
    """Oriented query position at which hit h reaches subject residue s_target (walking its BTOP)."""
    cur, sp = h.ostart, h.sstart
    if s_target <= sp:
        return cur
    for tok in TOKEN.findall(h.btop):
        if tok.isdigit():
            n = int(tok)
            if sp + n >= s_target:
                return cur + 3 * (s_target - sp)
            cur, sp = cur + 3 * n, sp + n
            continue
        qa, sa = tok[0], tok[1]
        if qa in "\\/":
            cur += 1 if qa == "\\" else -1
        elif qa == "-":
            sp += 1
        elif sa == "-":
            cur += 3
        else:
            cur, sp = cur + 3, sp + 1
        if sp >= s_target:
            return cur
    return cur


def inserts(kept: list[Hit], gene_of: dict[str, str]) -> list[dict]:
    """Unaligned inserts between consecutive HSPs of the same subject gene on one strand.
    Three cases (oriented coordinates):
      (1) query gap > 0 with the subject roughly continuous (-15..+30 residues): the gap;
      (2) the next HSP restarts at an earlier subject residue (overlap > 15 residues), as when the
          upstream HSP runs into an intron because the reference protein carries intron-encoded
          sequence: from where the upstream HSP reaches the downstream HSP's first residue, to the
          downstream HSP's start;
      (3) query gap > insert_min with subject residues missing (an exon not hit): the gap."""
    out = []
    by = {}
    for h in kept:
        by.setdefault((h.strand, gene_of.get(h.s, h.s)), []).append(h)
    for (strand, gene), hs in by.items():
        hs.sort(key=lambda h: h.ostart)
        for a, b in zip(hs, hs[1:]):
            gap = b.ostart - a.oend
            sgap = b.sstart - a.send
            left = None
            if gap > 0 and -15 <= sgap <= 30:
                left, case = a.oend, 1
            elif sgap < -15 and gap > -300:
                left, case = min(a.oend, subject_to_query(a, b.sstart)), 2
            elif gap > P["insert_min"] and sgap > 30:
                left, case = a.oend, 3
            if left is None:
                continue
            n = b.ostart - left
            if n <= 0:
                continue
            out.append({"strand": strand, "gene": gene, "left": left, "right": b.ostart, "len": n, "case": case,
                        "class": ("intron" if n > P["bypass_max"] else
                                  "bypass_like" if n >= P["bypass_min"] else "short")})
    return out
