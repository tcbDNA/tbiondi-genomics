"""Exploratory tests weighing H1 (translational skip) against H2 (disabling mutations) for Grigorovia jiainica.
Predictions were logged in PREREG.md before these ran. Outputs go to results/explore/grigorovia/h1h2/.

    uv run python code/explore/grigorovia_h1h2.py sister      # test 1: congener codons at the four sites
    uv run python code/explore/grigorovia_h1h2.py introns     # test 2a/2b: intron null and per-site junctions
"""

from __future__ import annotations

import csv
import subprocess
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Align, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fscore as F  # noqa: E402
import scan as S  # noqa: E402

ROOT = F.ROOT
OUT = ROOT / "results/explore/grigorovia/h1h2"
COMP = str.maketrans("ACGTN", "TGCAN")
T3 = F.codon_table(3)
SITES = {"cox1 24711": 24710, "cob 26496": 26495, "cob 26477": 26476, "cob 26266": 26265}  # 0-based fwd, minus strand
REF_GLY = {"cox1 24711": "G 8, F 1, A 1 (10 references)", "cob 26496": "G 10 (10)", "cob 26477": "G 10 (10)",
           "cob 26266": "G 6, S 2, N 1 (9) under the AGG placement; T 9, N 1 at the scan's placement"}


def rc(s):
    return s.translate(COMP)[::-1]


def aligner(mode="local"):
    a = Align.PairwiseAligner()
    a.mode = mode
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -5, -2
    return a


def load():
    gj = S.read_fasta(ROOT / "data/work/queries/code_3.fna", {"BK073075.1"})["BK073075.1"]
    yak = str(SeqIO.read(ROOT / "data/raw/OR703946/OR703946.1.gb", "genbank").seq).upper()
    return gj, yak


def best_local(q, target):
    al = aligner()
    best = None
    for strand, t in ((1, target), (-1, rc(target))):
        a = al.align(q, t)[0]
        if best is None or a.score > best[0].score:
            best = (a, strand, t)
    return best


def mapping(a):
    """query index -> target index for aligned (match/mismatch) columns."""
    m = {}
    for (qs, qe), (ts, te) in zip(*a.aligned):
        for i in range(qe - qs):
            m[qs + i] = ts + i
    return m


def _step_b_residues():
    """Residues of distant references at the codon read after the skip (results/explore/grigorovia/step_b_residue.md)."""
    out = {}
    txt = (ROOT / "results/explore/grigorovia/step_b_residue.md").read_text()
    for pos, counts in re.findall(r"## \w+ @ (\d+):.*?\n- reference residues at that position \(\d+ hits\): (\{[^}]*\})", txt):
        out[int(pos)] = counts
    return out


def a_covers(r):
    """The scan-placement alignment spans the window apart from <= 6 nt at either end (otherwise, e.g. at an intron in
    the congener, it covers only one side and the congener codons are not homologous across the site)."""
    return r["qs"] <= 24 and r.get("qe", 0) >= 36


def sister():
    """Test 1. For each site, remove each candidate base within 3 nt of the extra base, align the G. jiainica coding
    window (+-30 nt) to the congener OR703946 (both strands, local), and rank placements by score then mismatches.
    The congener codons are then read in G. jiainica's frame (best alignments are gapless)."""
    OUT.mkdir(parents=True, exist_ok=True)
    gj, yak = load()
    L = len(gj)
    q = rc(gj)
    refres = _step_b_residues()
    rows, detail = [], []
    for name, pf in SITES.items():
        o = L - 1 - pf
        win = q[o - 30:o + 31]                            # extra base at index 30
        res = []
        for k in range(27, 34):
            cand = win[:k] + win[k + 1:]
            a, strand, t = best_local(cand, yak)
            qb, tb = a.aligned
            (qs, qe), (ts, te) = qb[0], tb[0]
            mm = sum(1 for i in range(qe - qs) if cand[qs + i] != t[ts + i])
            res.append({"offset": k - 30, "base": win[k], "score": a.score, "mismatches": mm, "gapless": len(qb) == 1,
                        "cand": cand, "qs": qs, "qe": qb[-1][1], "ts": ts, "t": t, "strand": strand})
        res.sort(key=lambda r: (-r["score"], r["mismatches"]))
        top = [r for r in res if r["score"] == res[0]["score"] and r["mismatches"] == res[0]["mismatches"]]
        # codons in G. jiainica's old frame: the extra base is at win index 30; frame = codon starts at o - (o - f) % 3
        # where f is the gene frame on the coding strand (cox1/cob read in frame with the base removed)
        r0 = next(r for r in res if r["offset"] == 0)     # the scan's placement (the extra base itself removed)
        shift = r0["ts"] - r0["qs"]
        gseg = r0["cand"][30 - 12:30 + 15]
        yseg = r0["t"][30 - 12 + shift:30 + 15 + shift]
        comparable = r0["gapless"] and len(yseg) == len(gseg) and a_covers(r0)
        cod_gj = r0["cand"][30:33]
        cod_y = r0["t"][30 + shift:33 + shift] if comparable else ""
        rows.append({"site": name, "fwd_0based": pf, "gj_codons_base_removed": " ".join(gseg[i:i + 3] for i in range(0, len(gseg), 3)),
                     "congener_codons": " ".join(yseg[i:i + 3] for i in range(0, len(yseg), 3)) if comparable else "not comparable",
                     "gj_codon_after_skip": cod_gj, "gj_residue": T3.get(cod_gj, "X"),
                     "congener_codon": cod_y or "-", "congener_residue": T3.get(cod_y, "-") if cod_y else "-",
                     "distant_references": refres.get(pf, ""),
                     "best_placements": "; ".join(f"{r['base']} at {r['offset']:+d}" for r in top),
                     "best_score": top[0]["score"], "best_mismatches": top[0]["mismatches"],
                     "scan_placement_score": r0["score"], "scan_placement_mismatches": r0["mismatches"],
                     "scan_placement_residue_matches_congener": ("yes" if comparable and T3.get(cod_y) == T3.get(cod_gj) else
                                            f"no (congener {T3.get(cod_y, '-')})" if comparable else "not comparable")})
        detail += [f"## {name}", "", "Residues: codons 21-38 of the corrected window, read in the frame of the corrected gene "
                   "(codon starts at multiples of 3; the codon at 30 is the one right after the scan's extra base). "
                   "'stops' = stop codons anywhere in the 60-nt corrected window.", "",
                   "| offset | base removed | score | mismatches | gapless | codons 21-38 | residues | stops |", "|---|---|---|---|---|---|---|---|"]
        for r in res:
            c = r["cand"]
            cods = [c[i:i + 3] for i in range(21, 39, 3)]
            aa = "".join(T3.get(x, "X") for x in cods)
            st = sum(1 for i in range(0, 58, 3) if T3.get(c[i:i + 3], "") == "*")
            detail.append(f"| {r['offset']:+d} | {r['base']} | {r['score']:.0f} | {r['mismatches']} | {r['gapless']} | {' '.join(cods)} | {aa} | {st} |")
        if comparable:
            cods = [yseg[i:i + 3] for i in range(3, 21, 3)]          # yseg starts at window index 18
            detail.append(f"| congener | | | | | {' '.join(cods)} | {''.join(T3.get(x, 'X') for x in cods)} | |")
        detail.append("")
    with open(OUT / "test1_sister_codons.tsv", "w") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (OUT / "test1_placements.md").write_text("# Test 1: candidate placements of each extra base, scored against OR703946\n\n" +
                                             "Window +-30 nt on the coding strand; local alignment, match 2, mismatch -3, gap -5/-2. "
                                             "A site is comparable when the alignment of the window (scan placement removed) is "
                                             "gapless and spans window positions 24-36 around the site.\n\n" + "\n".join(detail) + "\n")
    for r in rows:
        print(r)


# ---------------------------------------------------------------- test 2: introns
SC = "NC_001224.1"
SITE_SC = {"cob 26496": ("cob", 408), "cob 26477": ("cob", 426), "cox1 24711": ("cox1", 231), "cob 26266": ("cob", 639)}
# extra base sits between S. cerevisiae CDS nt x and x+1 (x = 3*(residue-1)), residues 137, 143, 78, 214


def _taxonomy():
    nodes = {}
    with open(ROOT / "data/raw/taxdump/nodes.dmp") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            nodes[f[0]] = f[1]
    return nodes


def _in_taxon(taxid, target, nodes):
    t = str(taxid)
    for _ in range(60):
        if t == target:
            return True
        if t not in nodes or nodes[t] == t:
            return False
        t = nodes[t]
    return False


def _cds_records():
    seqs_needed, rows = set(), []
    G = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    nodes = _taxonomy()
    for r in csv.DictReader(open(ROOT / "data/work/cds.tsv"), delimiter="\t"):
        if r["gene"] not in ("cox1", "cob"):
            continue
        a = r["accession"]
        if a not in G or not _in_taxon(G[a]["taxid"], "147537", nodes):     # Saccharomycotina
            continue
        parts = [tuple(map(int, p.split("-"))) for p in r["parts"].split(";")]
        rows.append((a, r["gene"], int(r["strand"]), parts, G[a]))
        seqs_needed.add(a)
    return rows, seqs_needed, G


def introns():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, need, G = _cds_records()
    seqs = {}
    for c in {G[a]["code"] for a in need}:
        seqs.update(S.read_fasta(ROOT / f"data/work/queries/code_{c}.fna", need))
    from Bio.Align import substitution_matrices
    pal = Align.PairwiseAligner()
    pal.mode = "global"
    pal.substitution_matrix = substitution_matrices.load("BLOSUM62")
    pal.open_gap_score, pal.extend_gap_score = -10, -0.5
    sc_prot = {}
    out_pos = defaultdict(list)          # gene -> list of (S. cerevisiae CDS position, accession, organism)
    lens = {}
    for a, gene, strand, parts, meta in rows:
        s = seqs.get(a)
        if not s:
            continue
        ordered = sorted(parts) if strand == 1 else sorted(parts, reverse=True)
        cds = "".join(s[x:y] if strand == 1 else rc(s[x:y]) for x, y in ordered)
        table = F.codon_table(int(meta["code"]))
        prot = "".join(table.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3)).rstrip("*")
        bounds, cum = [], 0
        for x, y in ordered[:-1]:
            cum += y - x
            if True:
                bounds.append(cum)                       # intron after CDS nt `cum` (only gaps > 30 nt were joins? see below)
        gaps = [((ordered[i + 1][0] - ordered[i][1]) if strand == 1 else (ordered[i][0] - ordered[i + 1][1])) for i in range(len(ordered) - 1)]
        bounds = [b for b, g in zip(bounds, gaps) if g > 30]   # real introns, not 1-2 nt frameshift joins
        prev = out_pos.get((gene, a))
        if prev is not None and len(prev[1]) >= len(bounds):     # keep the CDS with the most introns (the full gene,
            continue                                                # not an intron-ORF fusion sharing only its first exons)
        if a == SC:
            sc_prot[gene] = prot
            lens[gene] = len(cds)
        out_pos[(gene, a)] = (prot, bounds, meta["organism"])
    res = defaultdict(list)
    for (gene, a), (prot, bounds, org) in out_pos.items():
        if not bounds or gene not in sc_prot:
            continue
        aln = pal.align(prot, sc_prot[gene])[0]
        m = {}
        for (qs, qe), (ts, te) in zip(*aln.aligned):
            for i in range(qe - qs):
                m[qs + i] = ts + i
        for b in bounds:
            r, ph = (b - 1) // 3, b % 3                      # intron after residue index r (0-based), phase
            if b % 3 == 0:
                r, ph = b // 3 - 1, 3
            if r in m:
                res[gene].append((3 * m[r] + ph, a, org))
    # null shares
    lines = ["# Test 2: introns (exploratory, post hoc)", "",
             "Intron insertion positions in cox1 and cob of RefSeq budding yeasts (Saccharomycotina, CDS joins with gaps > "
             "30 nt), mapped to S. cerevisiae CDS coordinates by global protein alignment (BLOSUM62).", ""]
    def share(positions, gene, w=12):
        L = lens[gene]
        near = {j for j in range(1, L) for p in positions if abs(j - p) <= w}   # j, p: junction after CDS nt j / p
        return len(near) / (L - 1)
    from math import comb
    def p_ge3(p):
        return sum(comb(4, k) * p ** k * (1 - p) ** (4 - k) for k in (3, 4))
    sc_only = {g: [p for p, a, _ in res[g] if a == SC] for g in ("cox1", "cob")}
    allpos = {g: sorted({p for p, a, _ in res[g]}) for g in ("cox1", "cob")}
    tot = lens["cox1"] + lens["cob"] - 2
    for label, P in (("S. cerevisiae alone", sc_only), ("all RefSeq budding yeasts with cox1/cob introns", allpos)):
        for w in (12, 2):
            sh = (share(P["cox1"], "cox1", w) * (lens["cox1"] - 1) + share(P["cob"], "cob", w) * (lens["cob"] - 1)) / tot
            lines.append(f"- **{label}, within {w} nt:** {len(P['cox1'])} cox1 and {len(P['cob'])} cob distinct intron positions; "
                         f"share of coding positions within {w} nt of one: **{sh:.3f}**; P(>= 3 of 4 random sites) = "
                         f"**{p_ge3(sh):.3g}**; P(>= 1 of 4) = {1 - (1 - sh) ** 4:.3g}")
    for label, P in (("S. cerevisiae alone", sc_only), ("all RefSeq budding yeasts", allpos)):
        ds = [int(min((abs(p - x) for p in P[g]), default=-1)) for g, x in SITE_SC.values()]
        for w in (12, 2):
            k = sum(1 for d in ds if d is not None and d <= w)
            sh = (share(P["cox1"], "cox1", w) * (lens["cox1"] - 1) + share(P["cob"], "cob", w) * (lens["cob"] - 1)) / tot
            pk = sum(comb(4, j) * sh ** j * (1 - sh) ** (4 - j) for j in range(k, 5))
            lines.append(f"- **Observed, {label}, within {w} nt:** nearest distances {ds} nt (sites in the order "
                         f"{', '.join(SITE_SC)}); {k} of 4 within {w} nt; P(>= {k} of 4) = **{pk:.3g}**")
    n_sp = len({a for g in res for _, a, _ in res[g]})
    lines += ["", f"Species with mapped introns: {n_sp}.", "", "## Per site: nearest intron positions (S. cerevisiae coordinates)", "",
              "| site | extra base between S. cerevisiae CDS nt | nearest S. cerevisiae intron | species with an intron within 12 nt (distance, organism) |",
              "|---|---|---|---|"]
    for site, (gene, x) in SITE_SC.items():
        # the extra base lies at the junction after S. cerevisiae CDS nt x; an intron after nt p is |p - x| nt away
        scd = min(((abs(p - x), p) for p in sc_only[gene]), default=(None, None))
        near = sorted({(abs(p - x), a, org) for p, a, org in res[gene] if abs(p - x) <= 12})
        lines.append(f"| {site} | {x}/{x + 1} | {('after nt ' + str(scd[1]) + ', ' + str(scd[0]) + ' nt') if scd[1] else '-'} | "
                     f"{'; '.join(f'{d} nt {org} ({a})' for d, a, org in near) or 'none'} |")
    (OUT / "test2_introns.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))



# ---------------------------------------------------------------- test 2b: introns mapped onto G. jiainica
GJ = "BK073075.1"
EXTRA = {"cox1": [24710], "cob": [26495, 26476, 26265]}
SITE_NAME = {24710: "cox1 24711", 26495: "cob 26496", 26476: "cob 26477", 26265: "cob 26266"}


def _blosum_aligner():
    from Bio.Align import substitution_matrices
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -11, -1
    return a


def gj_genes():
    """G. jiainica cox1/cob coding strand with the extra bases removed; returns per gene:
    (corrected coding nt string, frame, protein, {extra fwd pos: junction index J (extra base sits after corrected nt J)})."""
    gj, _ = load()
    G, ref_family, ref_gene = S.load_meta()
    hits = [h for h in F.read_hits(ROOT / "results/scan/fungal/hits_code3.tsv") if h.q == GJ and ref_gene.get(h.s) in EXTRA]
    out = {}
    for gene, extras in EXTRA.items():
        iv = [h.fwd_interval() for h in hits if ref_gene[h.s] == gene]
        a, b = min(x for x, _ in iv) - 60, max(y for _, y in iv) + 60
        fwd = list(range(b - 1, a - 1, -1))                  # coding-strand order of fwd positions (minus strand)
        keep = [p for p in fwd if p not in extras]
        seq = "".join(rc(gj[p]) for p in keep)
        J = {}
        for e in extras:
            J[e] = sum(1 for p in fwd[:fwd.index(e)] if p not in extras) - 1
        best = None
        for f in range(3):
            prot = "".join(T3.get(seq[i:i + 3], "X") for i in range(f, len(seq) - 2, 3))
            core = prot[20:-20]
            if best is None or core.count("*") < best[2].count("*"):
                best = (f, prot, core)
        out[gene] = (seq, best[0], best[1], J)
    return out


def map_junction(spliced, k, code, target_prot, target_frame, al):
    """Map the junction after spliced[:k] onto a target protein. Returns (target nt index the junction follows,
    score, aa aligned left, aa aligned right) or None."""
    table = F.codon_table(code)
    best = None
    for f in range(3):
        prot = "".join(table.get(spliced[i:i + 3], "X") for i in range(f, len(spliced) - 2, 3))
        if prot.count("*") > 1:
            continue
        a = al.align(prot, target_prot)[0]
        if best is None or a.score > best[0]:
            best = (a.score, f, prot, a)
    if best is None:
        return None
    score, f, prot, a = best
    m = {}
    for (qs, qe), (ts, te) in zip(*a.aligned):
        for i in range(qe - qs):
            m[qs + i] = ts + i
    kk, ph = (k - f) // 3, (k - f) % 3
    left = sum(1 for q in m if q < kk)
    right = sum(1 for q in m if q >= kk)
    if ph == 0:
        if kk - 1 not in m:
            return None
        t = target_frame + 3 * m[kk - 1] + 2
    else:
        if kk not in m:
            return None
        t = target_frame + 3 * m[kk] + ph - 1
    return t, score, left, right


def _local_dp_ends(q, t, ma=2, mi=-3, g=-4):
    """Linear-gap Smith-Waterman. Returns, for each query prefix length i (1..len q), the best score of an alignment
    ending exactly at q[i-1], and the target end (exclusive) of that alignment."""
    import numpy as np
    n, m = len(q), len(t)
    tb = np.frombuffer(t.encode(), dtype=np.uint8)
    prev = np.zeros(m + 1)
    best = np.zeros(n + 1)
    arg = np.zeros(n + 1, dtype=int)
    for i in range(1, n + 1):
        sub = np.where(tb == ord(q[i - 1]), ma, mi)
        cur = np.maximum(0, np.maximum(prev[:-1] + sub, prev[1:] + g))
        cur = np.concatenate(([0.0], cur))
        for jj in range(1, m + 1):                        # horizontal gaps (sequential)
            v = cur[jj - 1] + g
            if v > cur[jj]:
                cur[jj] = v
        diag_only = np.concatenate(([0.0], prev[:-1] + sub))  # alignment must END with q[i-1] matched/mismatched
        k = int(np.argmax(diag_only))
        best[i], arg[i] = diag_only[k], k
        prev = cur
    return best, arg


def spliced_junction(q, region1, off1, region2, off2):
    """Best single-intron spliced alignment of q (G. jiainica corrected coding window) to exon1 in region1 and exon2
    in region2 (coding orientation). Returns (split i in q, relative exon1 end, exon2 start, score, tied splits)."""
    b1, e1 = _local_dp_ends(q, region1)
    b2r, e2r = _local_dp_ends(q[::-1], region2[::-1])
    n = len(q)
    scores = [(b1[i] + b2r[n - i], i) for i in range(10, n - 9)]
    top = max(scores)[0]
    ties = [i for sc, i in scores if sc == top]
    i = ties[0]
    return i, off1 + int(e1[i]), off2 + (len(region2) - int(e2r[n - i])), top, ties


def introns_gj():
    """Per-site intron junctions mapped onto G. jiainica coordinates.
    (a) RefSeq Saccharomycotina with annotated cox1/cob introns: junctions from the annotation, mapped by aligning the
        annotated protein to G. jiainica's (extra bases removed). (b) The congener and the six closest relatives
        (unannotated TPA records): junctions refined by a single-intron spliced nucleotide alignment of G. jiainica's
        coding window to the scan's intron inserts (insert ends can be shifted where an intron ORF continues the exon)."""
    OUT.mkdir(parents=True, exist_ok=True)
    genes = gj_genes()
    al = _blosum_aligner()
    rows, need, G = _cds_records()
    seqs = {}
    for c in {G[a]["code"] for a in need}:
        seqs.update(S.read_fasta(ROOT / f"data/work/queries/code_{c}.fna", need))
    nodes = _taxonomy()
    out = []
    seen = {}
    for a, gene, strand, parts, meta in rows:
        s = seqs.get(a)
        if not s:
            continue
        ordered = sorted(parts) if strand == 1 else sorted(parts, reverse=True)
        gaps = [((ordered[i + 1][0] - ordered[i][1]) if strand == 1 else (ordered[i][0] - ordered[i + 1][1]))
                for i in range(len(ordered) - 1)]
        if sum(1 for g in gaps if g > 30) < 1:
            continue
        if (a, gene) in seen and seen[(a, gene)] >= len(ordered):
            continue
        seen[(a, gene)] = len(ordered)
        cds = "".join(s[x:y] if strand == 1 else rc(s[x:y]) for x, y in ordered)
        table = F.codon_table(int(meta["code"]))
        prot = "".join(table.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))
        seq, fr, gprot, J = genes[gene]
        aln = al.align(prot, gprot)[0]
        m = {}
        for (qs, qe), (ts, te) in zip(*aln.aligned):
            for i in range(qe - qs):
                m[qs + i] = ts + i
        cum = 0
        for (x, y), gp in zip(ordered[:-1], gaps):
            cum += y - x
            if gp <= 30:
                continue
            r, ph = cum // 3, cum % 3                     # junction after `cum` CDS nt: within codon r at phase ph
            if ph == 0:
                if r - 1 not in m:
                    continue
                t = fr + 3 * m[r - 1] + 2
            else:
                if r not in m:
                    continue
                t = fr + 3 * m[r] + ph - 1
            fam = "Saccharomycetaceae" if _in_taxon(meta["taxid"], "4893", nodes) else "other Saccharomycotina"
            for e, jj in J.items():
                if abs(t - jj) <= 30:
                    out.append({"site": SITE_NAME[e], "accession": a, "organism": meta["organism"], "group": fam,
                                "source": "RefSeq annotation", "distance_nt": abs(t - jj),
                                "side": "downstream" if t > jj else ("upstream" if t < jj else "at the extra base"),
                                "intron_len": gp, "note": ""})
    # (b) congener + closest relatives: refine scan inserts
    cand = {"OR703946.1": ROOT / "results/heldout/OR703946/inserts.tsv"}
    for acc in ("BK073182.1", "BK072186.1", "BK073486.1", "BK072188.1", "BK072187.1"):
        cand[acc] = ROOT / "results/scan/fungal/inserts.tsv"
    tseq = S.read_fasta(ROOT / "data/work/queries/code_3.fna", set(cand) - {"OR703946.1"})
    tseq["OR703946.1"] = str(SeqIO.read(ROOT / "data/raw/OR703946/OR703946.1.gb", "genbank").seq).upper()
    orgs = {"OR703946.1": "Grigorovia yakushimaensis (congener)", "BK073182.1": "Kazachstania bromeliacearum",
            "BK072186.1": "Arxiozyma bovina", "BK073486.1": "Monosporozyma solicola", "BK072188.1": "Arxiozyma telluris",
            "BK072187.1": "Arxiozyma heterogenica"}
    refined = []
    for acc, path in cand.items():
        s = tseq[acc]
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["accession"] != acc or r["gene"] not in EXTRA or r["class"] != "intron":
                continue
            L, R, st = int(r["left"]), int(r["right"]), int(r["strand"])
            if R - L > 10000:
                refined.append(f"| {orgs[acc]} | {r['gene']} | scan insert {L + 1}-{R} ({R - L} nt) | not refined: longer than 10 kb (gene split across the record ends or a rearrangement, not a single intron) | | |")
                continue
            seq, fr, gprot, J = genes[r["gene"]]
            so = s if st == 1 else rc(s)                    # inserts are in oriented coordinates
            r1, o1 = so[max(0, L - 2000):L + 60], max(0, L - 2000)
            r2, o2 = so[R - 60:R + 2000], R - 60
            for e, jj in J.items():
                lo, hi = max(0, jj - 150), min(len(seq), jj + 151)
                i, x1, x2, sc, ties = spliced_junction(seq[lo:hi], r1, o1, r2, o2)
                t = lo + i - 1                              # junction after corrected G. jiainica nt t
                tie_d = sorted({abs(lo + k - 1 - jj) for k in ties})
                junc = f"exon1 ends {x1}, exon2 starts {x2 + 1} (oriented); intron {x2 - x1} nt"
                if abs(t - jj) <= 30:
                    refined.append(f"| {orgs[acc]} | {r['gene']} {SITE_NAME[e]} | scan insert {L + 1}-{R} | {junc} | "
                                   f"**{abs(t - jj)} nt** {'downstream' if t > jj else 'upstream' if t < jj else 'at the extra base'} "
                                   f"(equal-score splits give {tie_d[0]}–{tie_d[-1]} nt) | {sc:.0f} |")
    out.sort(key=lambda x: (x["site"], x["distance_nt"]))
    with open(OUT / "test2b_annotated_introns_on_gj.tsv", "w") as fh:
        w = csv.DictWriter(fh, list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    from collections import Counter
    summ = []
    for site in SITE_NAME.values():
        for grp in ("Saccharomycetaceae", "other Saccharomycotina"):
            c = Counter(x["distance_nt"] for x in out if x["site"] == site and x["group"] == grp)
            summ.append(f"| {site} | {grp} | {', '.join(f'{d} nt: {n}' for d, n in sorted(c.items())) or 'none within 30 nt'} |")
    (OUT / "test2b_refined.md").write_text("\n".join(refined) + "\n")
    print("\n".join(summ)); print(); print("\n".join(refined))
    return out



def introns_rel():
    """Exon structures of cox1/cob in the congener and the six closest relatives from miniprot (S. cerevisiae proteins,
    code 3; records rotated where a gene straddles the record ends), each junction mapped onto G. jiainica by protein
    alignment, and junctions within 40 nt of a site refined by the single-intron spliced nucleotide alignment."""
    import subprocess, tempfile
    genes = gj_genes()
    al = _blosum_aligner()
    tmp = Path(tempfile.mkdtemp(dir=ROOT / "data/explore_tmp"))
    sc = {}
    for rec in SeqIO.parse(ROOT / "data/work/ref_proteins.faa", "fasta"):
        if rec.id in ("NC_001224.1|NP_009305.1", "NC_001224.1|NP_009315.1"):
            sc["cox1" if rec.id.endswith("305.1") else "cob"] = str(rec.seq)
    orgs = {"OR703946.1": "G. yakushimaensis (congener)", "BK073182.1": "K. bromeliacearum (1)",
            "NC_004918.1": "M. servazzii (2)", "BK072186.1": "A. bovina (3)", "BK073486.1": "M. solicola (4)",
            "BK072188.1": "A. telluris (5)", "BK072187.1": "A. heterogenica (6)"}
    seqs = S.read_fasta(ROOT / "data/work/queries/code_3.fna", set(orgs) - {"OR703946.1"})
    seqs["OR703946.1"] = str(SeqIO.read(ROOT / "data/raw/OR703946/OR703946.1.gb", "genbank").seq).upper()
    lines = ["| genome | gene | rotation | exons (coding strand, 1-based) | intron junctions mapped onto G. jiainica (nearest site, distance) |",
             "|---|---|---|---|---|"]
    refined = []
    for acc, org in orgs.items():
        for gene in ("cox1", "cob"):
            (tmp / "q.faa").write_text(f">{gene}\n{sc[gene]}\n")
            for rot in (0, len(seqs[acc]) // 2):
                g = seqs[acc][rot:] + seqs[acc][:rot]
                (tmp / "t.fna").write_text(f">t\n{g}\n")
                gff = subprocess.run(["miniprot", "-T", "3", "-t", "4", "--gff", "-N", "0", str(tmp / "t.fna"),
                                      str(tmp / "q.faa")], capture_output=True, text=True).stdout
                mrna = [l.split("\t") for l in gff.splitlines() if "\tmRNA\t" in l]
                if not mrna:
                    continue
                tgt = mrna[0][8].split("Target=")[1].split()
                if int(tgt[1]) <= 3 and int(tgt[2]) >= len(sc[gene]) - 3:
                    break
            cds = [l.split("\t") for l in gff.splitlines() if "\tCDS\t" in l and "Parent=MP000001" in l]
            strand = cds[0][6]
            ex = sorted((int(c[3]) - 1, int(c[4])) for c in cds)
            if strand == "-":
                ex = ex[::-1]
            gs = g if strand == "+" else rc(g)
            n = len(g)
            exo = [(x, y) if strand == "+" else (n - y, n - x) for x, y in ex]     # oriented 0-based
            spliced = "".join(gs[x:y] for x, y in exo)
            seq, fr, gprot, J = genes[gene]
            table = F.codon_table(3)
            prot = "".join(table.get(spliced[i:i + 3], "X") for i in range(0, len(spliced) - 2, 3))
            aln = al.align(prot, gprot)[0]
            m = {}
            for (qs, qe), (ts, te) in zip(*aln.aligned):
                for i in range(qe - qs):
                    m[qs + i] = ts + i
            cum, juncs = 0, []
            for k in range(len(exo) - 1):
                cum += exo[k][1] - exo[k][0]
                r, ph = cum // 3, cum % 3
                key = r - 1 if ph == 0 else r
                if key not in m:
                    juncs.append("unmapped")
                    continue
                t = fr + 3 * m[key] + (2 if ph == 0 else ph - 1)
                near = min(J.items(), key=lambda kv: abs(t - kv[1]))
                d = abs(t - near[1])
                txt = f"{SITE_NAME[near[0]]}: {d} nt" + (" (rough)" if d <= 40 else "")
                if d <= 40:
                    lo, hi = max(0, near[1] - 150), min(len(seq), near[1] + 151)
                    r1, o1 = gs[max(0, exo[k][1] - 400):exo[k][1] + 60], max(0, exo[k][1] - 400)
                    r2, o2 = gs[exo[k + 1][0] - 60:exo[k + 1][0] + 400], exo[k + 1][0] - 60
                    i, x1, x2, sc_, ties = spliced_junction(seq[lo:hi], r1, o1, r2, o2)
                    tt = lo + i - 1
                    tie_d = sorted({abs(lo + q - 1 - near[1]) for q in ties})
                    ctx = f"`{gs[x1 - 9:x1]}|{gs[x1:x1 + 6]}…{gs[x2 - 6:x2]}|{gs[x2:x2 + 9]}`"
                    txt = (f"{SITE_NAME[near[0]]}: **{abs(tt - near[1])} nt** {'downstream' if tt > near[1] else 'upstream' if tt < near[1] else 'at the extra base'}"
                           f" (equal-score splits {tie_d[0]}–{tie_d[-1]}); refined junction {ctx}")
                juncs.append(txt)
            exs = ", ".join(f"{x + 1}-{y}" for x, y in exo)
            lines.append(f"| {org} ({acc}) | {gene} | {rot} | {exs} (oriented) | {'; '.join(juncs) or 'no intron'} |")
    (OUT / "test2b_relatives_miniprot.md").write_text(
        "# Test 2b addendum: cox1/cob exon structures in the congener and the six closest relatives (miniprot 0.18, -T 3)\n\n"
        "Coordinates are on the coding strand of the (rotated) record. Distances are in G. jiainica coordinates "
        "(extra bases removed); 0 = the intron sits exactly where G. jiainica has its extra base. 'Refined' junctions use "
        "the spliced nucleotide alignment; ties are splits with equal score (sequence alone cannot place them).\n\n"
        + "\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------------------------------------------------------------- test 2b extras: intron ORFs, Rfam, M. servazzii
TMP = ROOT / "data/explore_tmp"
INTRON_ORF = re.compile(r"laglidadg|giy|endonuc|maturase|^other:a?i\d|^other:bi\d|intron|^other:orf", re.I)


def _run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw).stdout


def intron_orfs():
    """(i) G. jiainica's whole mitogenome against intron-ORF-like reference proteins (and S. cerevisiae exon-intron fusion
    ORFs), DIAMOND --more-sensitive, frameshift-aware, range culling, E <= 1; hits outside core genes and rps3 listed.
    (ii) The congener's cob intron: in-frame ORF from the exon-1 end, DIAMOND hits over it, Rfam group I/II models."""
    TMP.mkdir(parents=True, exist_ok=True)
    ids = set()
    gene = {}
    for r in csv.DictReader(open(ROOT / "data/work/ref_proteins.tsv"), delimiter="\t"):
        gene[r["ref_id"]] = r
        if INTRON_ORF.search(r["gene"]) or (r["accession"] == "NC_001224.1" and r["gene"] in ("cob", "cox1") and int(r["length"]) > 600):
            ids.add(r["ref_id"])
    with open(TMP / "intron_orf_refs.faa", "w") as fh:
        for rec in SeqIO.parse(ROOT / "data/work/ref_proteins.faa", "fasta"):
            if rec.id in ids:
                fh.write(f">{rec.id}\n{rec.seq}\n")
    _run(["diamond", "makedb", "--quiet", "--in", str(TMP / "intron_orf_refs.faa"), "-d", str(TMP / "intron_orf_refs")])
    gj, yak = load()
    (TMP / "gj.fna").write_text(f">{GJ}\n{gj}\n")
    out = _run(["diamond", "blastx", "--quiet", "--more-sensitive", "-q", str(TMP / "gj.fna"), "-d", str(TMP / "intron_orf_refs"),
                "--query-gencode", "3", "-F", "15", "--range-culling", "-k", "5", "-e", "1", "--outfmt", "6",
                "qseqid", "sseqid", "qstart", "qend", "evalue", "bitscore", "pident"])
    G, ref_family, ref_gene = S.load_meta()
    hits = [h for h in F.read_hits(ROOT / "results/scan/fungal/hits_code3.tsv") if h.q == GJ]
    kept = F.cull(hits, G[GJ]["family_pool"], ref_family)
    spans = {}
    for h in kept:
        g = ref_gene[h.s]
        if g.startswith("other:") and g not in ("other:rps3", "other:VAR1"):
            continue                                   # count only core genes and rps3/VAR1 as "inside"
        a, b = h.fwd_interval()
        lo, hi = spans.get(g, (10 ** 9, 0))
        spans[g] = (min(lo, a), max(hi, b))
    lines = ["# Test 2b: intron-ORF remnants and the congener's cob intron", "",
             f"Reference set: {len(ids)} intron-ORF-like proteins (gene labels matching LAGLIDADG, GIY-YIG, endonuclease, maturase, "
             "AI/BI ORFs, intron, unnamed ORFs; plus S. cerevisiae exon-intron fusion proteins carried under cob/cox1).", "",
             "G. jiainica gene spans (scan alignments, forward coordinates): " + "; ".join(f"{g} {a + 1}-{b}" for g, (a, b) in sorted(spans.items(), key=lambda x: x[1])), "",
             "| hit (forward, 1-based) | reference gene | E | bits | identity | inside |", "|---|---|---|---|---|---|"]
    n_all = n_out = 0
    for l in out.splitlines():
        f = l.split("\t")
        a, b = sorted((int(f[2]), int(f[3])))
        n_all += 1
        inside = [g for g, (x, y) in spans.items() if min(b, y) - max(a - 1, x) > 0]
        if not inside:
            n_out += 1
        lines.append(f"| {a}-{b} | {gene[f[1]]['gene']} | {f[4]} | {f[5]} | {f[6]} | {', '.join(inside) or '**outside core genes and rps3**'} |")
    lines += ["", f"Hits: {n_all}; outside core genes and rps3: {n_out}.", ""]
    # congener intron
    ex = exon_struct_simple(yak, "cob")
    e1, e2 = ex[0][1], ex[1][0]                       # exon 1 end (excl.), exon 2 start (0-based, + strand)
    i0 = e1 - 9                                        # last full codon of exon 1 (CAT TGG GGT ends at 21281, group I rule)
    orf, k = "", e1 - 9 + 9
    # in-frame read-through from the codon after the last exon-1 codon under the group I rule
    lines += ["## The congener's cob intron (OR703946)", "",
              f"- miniprot exon structure: exon 1 ends at {e1}, exon 2 starts at {e2 + 1} (1-based; {e2 - e1} nt between)."]
    refined = spliced_for_congener(yak)
    lines.append(f"- spliced nucleotide alignment of G. jiainica's coding window: junction {refined}")
    T3_ = T3
    start = refined_exon1_end(yak) - 9
    prot = ""
    for j in range(start + 9, e2, 3):
        aa = T3_.get(yak[j:j + 3], "X")
        prot += aa
        if aa == "*":
            break
    lines.append(f"- open frame continuing from exon 1 into the intron (group I rule exon end): {len(prot)} codons to the first stop")
    oh = [h for h in F.read_hits(ROOT / "results/heldout/OR703946/hits_code3.tsv") if min(h.fwd_interval()[1], e2) - max(h.fwd_interval()[0], e1) > 300]
    best = {}
    for h in oh:
        g = ref_gene.get(h.s, "?")
        if g not in best or h.bits > best[g].bits:
            best[g] = h
    lines.append("- scan alignments over the intron (> 300 nt inside it): " +
                 "; ".join(f"{g} {h.s.split('|')[1]} ({gene.get(h.s, {}).get('length', '?')} aa) bits {h.bits:.0f}" for g, h in sorted(best.items(), key=lambda x: -x[1].bits)))
    fa = TMP / "or703946_cob_intron.fa"
    fa.write_text(f">intron\n{yak[e1 - 30:e2 + 30]}\n")
    for m in ("RF00028", "RF00029"):
        _run([str(ROOT / "tools/local/bin/cmsearch"), "--cpu", "4", "--max", "-E", "1000", "--tblout", str(TMP / f"{m}.tbl"),
              str(ROOT / f"data/rfam/{m}.cm"), str(fa)])
        top = [l.split() for l in open(TMP / f"{m}.tbl") if not l.startswith("#")]
        lines.append(f"- Rfam {m} ({'group I' if m == 'RF00028' else 'group II'} intron model), cmsearch --max: best score "
                     f"{top[0][14] if top else '-'}, E {top[0][15] if top else '-'}")
    (OUT / "test2b_intron_orfs.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-8:]))


def exon_struct_simple(seq, gene):
    """miniprot exon structure (+ strand records only; used for the congener's cob)."""
    import tempfile
    tmp = Path(tempfile.mkdtemp(dir=TMP))
    q = next(str(r.seq) for r in SeqIO.parse(ROOT / "data/work/ref_proteins.faa", "fasta") if r.id == ("NC_001224.1|NP_009315.1" if gene == "cob" else "NC_001224.1|NP_009305.1"))
    (tmp / "q.faa").write_text(f">q\n{q}\n")
    (tmp / "t.fna").write_text(f">t\n{seq}\n")
    gff = _run(["miniprot", "-T", "3", "-t", "2", "--gff", "-N", "0", str(tmp / "t.fna"), str(tmp / "q.faa")])
    cds = sorted((int(l.split("\t")[3]) - 1, int(l.split("\t")[4])) for l in gff.splitlines() if "\tCDS\t" in l and "Parent=MP000001" in l)
    return cds


def spliced_for_congener(yak):
    genes = gj_genes()
    seq, fr, gprot, J = genes["cob"]
    ex = exon_struct_simple(yak, "cob")
    e1, e2 = ex[0][1], ex[1][0]
    jj = J[26476]
    lo, hi = max(0, jj - 150), min(len(seq), jj + 151)
    i, x1, x2, sc, ties = spliced_junction(seq[lo:hi], yak[e1 - 400:e1 + 60], e1 - 400, yak[e2 - 60:e2 + 400], e2 - 60)
    td = sorted({abs(lo + k - 1 - jj) for k in ties})
    t = lo + i - 1
    return f"{abs(t - jj)} nt {'downstream' if t > jj else 'upstream'} of cob 26477 (equal-score splits {td[0]}-{td[-1]} nt); " \
           f"`{yak[x1 - 9:x1]}|{yak[x1:x1 + 6]}...{yak[x2 - 6:x2]}|{yak[x2:x2 + 9]}`; intron length {x2 - x1} nt " \
           f"(the tied split one base later, which follows the group I rule, gives the same length)"


def refined_exon1_end(yak):
    """Exon 1 end (0-based, exclusive) of the congener's cob under the group I rule: last exon base T, intron ends in G.
    The spliced alignment ties two splits one base apart; the rule picks the one ending ...GGT|AT...TTG|GCA."""
    ex = exon_struct_simple(yak, "cob")
    e1, e2 = ex[0][1], ex[1][0]
    for d in range(-3, 4):
        x1 = e1 + d
        if yak[x1 - 1] == "T" and yak[e2 - 1] == "G" and yak[x1 - 9:x1].endswith("GGT"):
            return x1
    return e1


def servazzii_skips():
    """M. servazzii cob (NC_004918.1, RefSeq annotation): its intron and two 1-nt skips mapped onto G. jiainica."""
    genes = gj_genes()
    al = _blosum_aligner()
    rows, need, G = _cds_records()
    s = S.read_fasta(ROOT / "data/work/queries/code_3.fna", {"NC_004918.1"})["NC_004918.1"]
    parts = next(p for a, g, st, p, m in rows if a == "NC_004918.1" and g == "cob")
    parts = sorted(parts)
    cds = "".join(s[x:y] for x, y in parts)
    prot = "".join(T3.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))
    seq, fr, gprot, J = genes["cob"]
    a = al.align(prot, gprot)[0]
    m = {}
    for (qs, qe), (ts, te) in zip(*a.aligned):
        for i in range(qe - qs):
            m[qs + i] = ts + i
    lines = ["# M. servazzii cob (NC_004918.1): annotated intron and 1-nt skips mapped onto G. jiainica", "",
             f"Annotated CDS parts (0-based, half-open): {parts}", "",
             "| junction | after M. servazzii CDS nt | skipped base context | G. jiainica coding nt it maps after | distance to cob 26496 / 26477 / 26266 |", "|---|---|---|---|---|"]
    cum = 0
    for (x, y), (x2, y2) in zip(parts, parts[1:]):
        cum += y - x
        gap = x2 - y
        r, ph = cum // 3, cum % 3
        t = fr + 3 * m[r - 1] + 2 if ph == 0 else fr + 3 * m[r] + ph - 1
        kind = "intron" if gap > 30 else f"{gap}-nt skip"
        ctx = f"`{s[y - 9:y]}[{s[y:x2]}]{s[x2:x2 + 6]}`" if gap <= 30 else ""
        lines.append(f"| {kind} | {cum} | {ctx} | {t} | " + " / ".join(str(t - J[e]) for e in (26495, 26476, 26265)) + " |")
    (OUT / "test2b_servazzii_cob.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


# ---------------------------------------------------------------- test 5: rest of the mitogenome
def mitogenome():
    """Protein genes (miniprot, S. cerevisiae proteins, code 3; relaxed for atp8), VAR1 open frame, rRNAs and RPM1
    (local alignment of S. cerevisiae sequences from NC_001224.1), scan events and stop columns, tRNA sets."""
    TMP.mkdir(parents=True, exist_ok=True)
    gj, yak = load()
    names = {"NP_009305.1": "cox1", "NP_009312.1": "atp8", "NP_009313.1": "atp6", "NP_009315.1": "cob", "NP_009319.1": "atp9",
             "NP_009320.1": "VAR1", "NP_009326.1": "cox2", "NP_009328.1": "cox3"}
    with open(TMP / "sc_mito.faa", "w") as fh:
        for r in SeqIO.parse(ROOT / "data/work/ref_proteins.faa", "fasta"):
            acc = r.id.split("|")[-1]
            if r.id.startswith("NC_001224.1") and acc in names:
                fh.write(f">{names[acc]}_{acc}\n{r.seq}\n")
    _run(["diamond", "makedb", "--quiet", "--in", str(TMP / "sc_mito.faa"), "-d", str(TMP / "sc_mito")])
    (TMP / "yak.fna").write_text(f">OR703946.1\n{yak}\n")
    rows = []
    for label, fa in (("G. jiainica", TMP / "gj.fna"), ("congener", TMP / "yak.fna")):
        if label == "G. jiainica":
            fa.write_text(f">{GJ}\n{gj}\n")
        for relaxed in (False, True):
            extra = ["--outs=0.3", "-p", "0.1", "-n", "1", "-l", "3"] if relaxed else []
            out = _run(["miniprot", "-T", "3", "-t", "4", "-N", "3", *extra, str(fa), str(TMP / "sc_mito.faa")])
            for l in out.splitlines():
                f = l.split("\t")
                tags = {t.split(":")[0]: t.split(":")[-1] for t in f[12:]}
                g = f[0].split("_")[0]
                if g == "VAR1":                        # too divergent for miniprot (relaxed runs give spurious spliced hits)
                    continue
                if any(r["genome"] == label and r["gene"] == g for r in rows):
                    continue
                rows.append({"genome": label, "gene": g, "settings": "relaxed" if relaxed else "default", "qlen": int(f[1]),
                             "q_aligned": f"{int(f[2]) + 1}-{f[3]}", "target": f"{int(f[7]) + 1}-{f[8]} ({f[4]})",
                             "frameshifts": int(tags.get("fs", 0)), "stops": int(tags.get("st", 0))})
    lines = ["# Test 5: the rest of the mitogenome (G. jiainica BK073075.1 and the congener OR703946.1)", "",
             "## Protein genes (miniprot 0.18, -T 3; S. cerevisiae proteins; best alignment per gene; 'relaxed' = --outs=0.3 -p 0.1 -n 1 -l 3; "
             "VAR1 is assessed by its open frame below, as it is too divergent for miniprot)", "",
             "| genome | gene | settings | query length | aligned | target | frameshifts | stops |", "|---|---|---|---|---|---|---|---|"]
    lines += [f"| {r['genome']} | {r['gene']} | {r['settings']} | {r['qlen']} | {r['q_aligned']} | {r['target']} | {r['frameshifts']} | {r['stops']} |"
              for r in sorted(rows, key=lambda r: (r["genome"], r["gene"]))]
    # VAR1 open frame: longest stop-to-stop frame overlapping the best VAR1/rps3 scan hit
    G, ref_family, ref_gene = S.load_meta()

    def var1(acc, seq, hits_path):
        hs = [h for h in F.read_hits(hits_path) if h.q == acc and ref_gene.get(h.s, "") in ("other:rps3", "other:VAR1")]
        h = max(hs, key=lambda h: h.bits)
        a, b = h.fwd_interval()
        best = (0, 0, 0)
        for strand in (1, -1):
            sq = seq if strand == 1 else rc(seq)
            A, B = (a, b) if strand == 1 else (len(seq) - b, len(seq) - a)
            for f in range(3):
                i = f
                while i < len(sq) - 2:
                    j = i
                    while j < len(sq) - 2 and T3.get(sq[j:j + 3], "X") != "*":
                        j += 3
                    ov = min(j, B) - max(i, A)
                    if ov > best[0]:
                        best = (ov, (j - i) // 3, strand)
                    i = j + 3
        return f"{ref_gene[h.s]} hit {a + 1}-{b} (bits {h.bits:.0f}); longest open frame over it: {best[1]} codons (strand {best[2]:+d})"
    lines += ["", "## VAR1 (rps3)", "",
              f"- G. jiainica: {var1(GJ, gj, ROOT / 'results/scan/fungal/hits_code3.tsv')}",
              f"- congener: {var1('OR703946.1', yak, ROOT / 'results/heldout/OR703946/hits_code3.tsv')}"]
    # rRNAs and RPM1
    gb = TMP / "NC_001224.1.gb"
    if not gb.exists():
        _run(["curl", "-sf", "-o", str(gb), "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001224.1&rettype=gbwithparts&retmode=text"])
    rec = SeqIO.read(gb, "genbank")
    feats = [((f.qualifiers.get("gene") or f.qualifiers.get("product"))[0], str(f.extract(rec.seq)).upper())
             for f in rec.features if f.type in ("rRNA", "ncRNA", "misc_RNA")]
    al = aligner()
    lines += ["", "## RNA genes (local alignment of S. cerevisiae sequences; match 2, mismatch -3, gap -5/-2)", "",
              "| genome | gene (S. cerevisiae length) | score | strand | coverage of the S. cerevisiae sequence | identity |", "|---|---|---|---|---|---|"]
    for label, sq in (("G. jiainica", gj), ("congener", yak)):
        for name, qq in feats:
            best = None
            for strand, t in ((1, sq), (-1, rc(sq))):
                a = al.align(qq, t)[0]
                if best is None or a.score > best[0]:
                    best = (a.score, strand, a, t)
            sc, strand, a, t = best
            n = sum(e - b for b, e in a.aligned[0])
            ident = sum(1 for (qs, qe), (ts, te) in zip(*a.aligned) for i in range(qe - qs) if qq[qs + i] == t[ts + i]) / max(1, n)
            lines.append(f"| {label} | {name} ({len(qq)}) | {sc:.0f} | {strand:+d} | {n / len(qq):.2f} | {ident:.2f} |")
    # scan events and stop columns
    ev = [r for r in csv.DictReader(open(ROOT / "results/scan/fungal/events.tsv"), delimiter="\t") if r["accession"] == GJ]
    st = [r for r in csv.DictReader(open(ROOT / "results/scan/fungal/stop_columns.tsv"), delimiter="\t") if r["accession"] == GJ]
    lines += ["", "## Every scan event and stop column in G. jiainica", "",
              "| type | gene | forward position (0-based) | shift / residue | codon | class |", "|---|---|---|---|---|---|"]
    lines += [f"| event | {r['gene']} | {r['pos_fwd']} | {r['shift']} | {r['codon']} | {r['class'] or '-'} |" for r in ev]
    lines += [f"| stop column | {r['gene']} | {r['pos']} | {r['subject_aa']} | | |" for r in st]
    # tRNAs
    def parse(path):
        out = []
        for l in open(path):
            f = l.split()
            if len(f) >= 9 and f[1].isdigit():
                out.append((f[4], f[5], float(f[8]), " ".join(f[9:])))
        return out
    a_ = parse(ROOT / "results/explore/grigorovia/trnascan_BK073075.1.txt")
    b_ = parse(ROOT / "results/heldout/OR703946/trnascan.txt")
    from collections import Counter
    ca, cb = Counter((t, x) for t, x, _, _ in a_), Counter((t, x) for t, x, _, _ in b_)
    lines += ["", "## tRNAs (tRNAscan-SE 2.0.12 -O)", "", f"G. jiainica {len(a_)} tRNAs, congener {len(b_)}; notes (pseudogene etc.): "
              f"{sum(1 for x in a_ + b_ if x[3])}.", "", "| type | anticodon | G. jiainica (scores) | congener (scores) |", "|---|---|---|---|"]
    for k in sorted(set(ca) | set(cb)):
        lines.append(f"| {k[0]} | {k[1]} | {ca.get(k, 0)} ({', '.join(f'{s:.1f}' for t, x, s, n in a_ if (t, x) == k)}) | "
                     f"{cb.get(k, 0)} ({', '.join(f'{s:.1f}' for t, x, s, n in b_ if (t, x) == k)}) |")
    sa = [s for _, _, s, _ in a_]; sb = [s for _, _, s, _ in b_]
    lines.append(f"\nScore ranges: G. jiainica {min(sa):.1f}-{max(sa):.1f}; congener {min(sb):.1f}-{max(sb):.1f}.")
    (OUT / "test5_tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def reads_strand():
    """Reads at the four sites by strand (exact 20 + base + 20 nt match, with or without the extra base), and reads spanning
    both cob sites 26496 and 26477 (19 nt apart) with both, one or neither base. SRR16974510, both read files."""
    import gzip
    gj, _ = load()
    pats = {}
    for name, p in SITES.items():
        w = gj[p - 20:p + 21]
        pats[name] = {"with": w, "without": w[:20] + w[21:]}
    a, b = 26476, 26495                               # 0-based forward positions (26477, 26496)
    joint = gj[a - 15:b + 16]
    ia, ib = 15, 15 + (b - a)
    variants = {"both bases": joint, "lacks 26477 only": joint[:ia] + joint[ia + 1:],
                "lacks 26496 only": joint[:ib] + joint[ib + 1:], "lacks both": joint[:ia] + joint[ia + 1:ib] + joint[ib + 1:]}
    comp = str.maketrans("ACGTN", "TGCAN")
    counts = defaultdict(int)
    for f in (ROOT / "data/reads/SRR16974510/SRR16974510_1.fastq.gz", ROOT / "data/reads/SRR16974510/SRR16974510_2.fastq.gz"):
        with gzip.open(f, "rt") as fh:
            for n, line in enumerate(fh):
                if n % 4 != 1:
                    continue
                r = line.strip()
                rr = r.translate(comp)[::-1]
                for name, d in pats.items():
                    for kind, w in d.items():
                        if w in r:
                            counts[(name, kind, "forward")] += 1
                        if w in rr:
                            counts[(name, kind, "reverse")] += 1
                for kind, w in variants.items():
                    if w in r or w in rr:
                        counts[("joint", kind)] += 1
    lines = ["# Reads at the four sites, by strand (SRR16974510; exact match of 20 nt + base + 20 nt, or of the 40 nt without the base)", "",
             "'forward' = the read matches the BK073075.1 forward strand; 'reverse' = its reverse complement.", "",
             "| site | with the base, forward | with, reverse | without, forward | without, reverse |", "|---|---|---|---|---|"]
    for name in SITES:
        lines.append(f"| {name} | {counts[(name, 'with', 'forward')]} | {counts[(name, 'with', 'reverse')]} | "
                     f"{counts[(name, 'without', 'forward')]} | {counts[(name, 'without', 'reverse')]} |")
    lines += ["", "## Reads spanning both cob sites (26477 and 26496; window 26462-26511)", "", "| version | reads |", "|---|---|"]
    lines += [f"| {k} | {counts[('joint', k)]} |" for k in variants]
    (OUT / "reads_strand.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"sister": sister, "introns": introns, "introns_gj": introns_gj, "introns_rel": introns_rel,
     "intron_orfs": intron_orfs, "servazzii_skips": servazzii_skips, "mitogenome": mitogenome,
     "reads_strand": reads_strand}[sys.argv[1]]()
