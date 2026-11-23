"""Two quick checks requested after the outside review (post hoc, exploratory).

  trna   tRNA-Gly genes (ARAGORN 1.2.41, yeast mitochondrial code) in the four Grigorovia mitogenomes and the outgroups
         used for dN/dS (M. servazzii, A. bovina, A. telluris) plus K. bromeliacearum; each aligned globally to the
         G. jiainica gene: identity and every difference.
  genes  every mitochondrial protein gene of G. jiainica: codons, glycine codons by type, extra bases (miniprot
         frameshifts, test 5), internal stops, and branch dN/dS for the G. jiainica branch (codeml, same CDS, codon
         alignments, tree and settings as `grigorovia_genus.py dnds`: F3x4, yeast mitochondrial code, cleandata = 1).
         VAR1 is too divergent to align (test 5); its open frame is counted only.

    uv run python src/explore/quick_checks.py trna    -> results/explore/grigorovia/h1h2/trna_gly_relatives.md
    uv run python src/explore/quick_checks.py genes   -> results/explore/grigorovia/h1h2/genus/genes_per_gene.md
"""
import re
import subprocess
import sys
import tempfile
from collections import Counter
from math import erfc, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src/explore"))
sys.path.insert(0, str(ROOT / "src"))
import grigorovia_genus as GG  # noqa: E402

H = ROOT / "results/explore/grigorovia/h1h2"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(x):
    return x.translate(COMP)[::-1]


def aragorn_gly(seq):
    """tRNA-Gly genes in one sequence: list of (sequence 5'->3', anticodon, coordinates)."""
    out = []
    with tempfile.TemporaryDirectory() as d:
        Path(d, "s.fa").write_text(f">s\n{seq}\n")
        txt = subprocess.run([str(ROOT / "tools/bin/aragorn"), "-w", "-t", "-gc3", "-l", f"{d}/s.fa"], capture_output=True,
                             text=True).stdout
    for m in re.finditer(r"tRNA-Gly\s+(c?)\[(\d+),(\d+)\]\s+\d+\s+\((\w+)\)", txt):
        c, a, b, ac = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4).upper()
        s = seq[a - 1:b]
        out.append((rc(s) if c else s, ac, f"{'c' if c else ''}[{a},{b}]"))
    return out


def nw(a, b, match=1, mis=-1, gap=-2):
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = i * gap
    for j in range(1, m + 1):
        S[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            S[i][j] = max(S[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mis), S[i - 1][j] + gap, S[i][j - 1] + gap)
    i, j, A, B = n, m, [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i][j] == S[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mis):
            A.append(a[i - 1]); B.append(b[j - 1]); i -= 1; j -= 1
        elif i > 0 and S[i][j] == S[i - 1][j] + gap:
            A.append(a[i - 1]); B.append("-"); i -= 1
        else:
            A.append("-"); B.append(b[j - 1]); j -= 1
    return "".join(reversed(A)), "".join(reversed(B))


def trna():
    seqs = GG.mito_seqs()
    og = GG.S.read_fasta(ROOT / "data/work/queries/code_3.fna", {"NC_004918.1", "BK072186.1", "BK072188.1", "BK073182.1"})
    names = {"NC_004918.1": "M. servazzii", "BK072186.1": "A. bovina", "BK072188.1": "A. telluris",
             "BK073182.1": "K. bromeliacearum"}
    genomes = {}
    for sp, v in seqs.items():
        genomes[f"G. {sp}"] = v if isinstance(v, dict) else {sp: v}
    for a, n in names.items():
        genomes[n] = {a: og[a]}
    found = {}
    for g, d in genomes.items():
        hits = []
        for cid, s in d.items():
            for t in aragorn_gly(s):
                hits.append((cid,) + t)
        found[g] = hits
    ref = found["G. jiainica"][0][1]
    L = ["# tRNA-Gly genes in G. jiainica and its relatives (ARAGORN 1.2.41, yeast mitochondrial code; post hoc)", "",
         f"G. jiainica tRNA-Gly ({len(ref)} nt, anticodon {found['G. jiainica'][0][2]}): `{ref}`", "",
         "| genome | tRNA-Gly genes | anticodon | length | identity to G. jiainica | differences (G. jiainica position: its base > other) |",
         "|---|---|---|---|---|---|"]
    for g, hits in found.items():
        if not hits:
            L.append(f"| {g} | 0 | | | | |")
            continue
        for cid, s, ac, co in hits:
            A, B = nw(ref, s)
            pos, diffs, same, cols = 0, [], 0, 0
            for x, y in zip(A, B):
                if x != "-":
                    pos += 1
                if x != "-" and y != "-":
                    cols += 1
                    same += x == y
                if x != y:
                    diffs.append(f"{pos}{'' if x != '-' else '+'}:{x}>{y}")
            ident = same / max(len(A), 1)
            L.append(f"| {g} | {len(hits)} | {ac} | {len(s)} | {ident:.1%} ({same}/{len(A)} columns) | "
                     f"{', '.join(diffs) if diffs else 'none'} |")
    L += ["", "Positions count along the G. jiainica gene (5'->3'); '+' marks a base the relative has and G. jiainica lacks."]
    # positions where G. jiainica differs from every relative, and pairwise differences among the four Grigorovia
    per = {}
    for g, hits in found.items():
        if g == "G. jiainica" or not hits:
            continue
        A, B = nw(ref, hits[0][1])
        pos, d = 0, set()
        for x, y in zip(A, B):
            if x != "-":
                pos += 1
                if y != x:
                    d.add(pos)
        per[g] = d
    allrel = set.intersection(*per.values())
    cong = set.intersection(*(v for k, v in per.items() if k.startswith("G. ")))
    L += ["", f"Positions where G. jiainica differs from all {len(per)} relatives: {len(allrel)} ({', '.join(map(str, sorted(allrel)))}); "
          f"from all three congeners: {len(cong)} ({', '.join(map(str, sorted(cong)))}).", "",
          "## Pairwise differences among the four Grigorovia tRNA-Gly genes (substitutions / end indels, global alignment)", "",
          "| | " + " | ".join(g for g in found if g.startswith("G. ")) + " |", "|---|" + "---|" * sum(g.startswith("G. ") for g in found)]
    gl = [g for g in found if g.startswith("G. ") and found[g]]
    for g1 in gl:
        row = []
        for g2 in gl:
            A, B = nw(found[g1][0][1], found[g2][0][1])
            sub = sum(1 for x, y in zip(A, B) if x != y and "-" not in (x, y))
            ind = sum(1 for x, y in zip(A, B) if "-" in (x, y))
            row.append("—" if g1 == g2 else f"{sub} / {ind}")
        L.append(f"| {g1} | " + " | ".join(row) + " |")
    # stem pairing at the changed positions, on G. jiainica's tRNAscan-SE structure
    ss = (ROOT / "results/explore/grigorovia/trnascan_ss_BK073075.1.txt").read_text()
    m = re.search(r"Type: Gly.*?\n.*?\nSeq: (\S+)\nStr: (\S+)", ss, re.S)
    sq, st = m.group(1), m.group(2)
    assert sq == ref, "tRNAscan-SE and ARAGORN sequences differ"
    stack, pair = [], {}
    for i, ch in enumerate(st):
        if ch == ">":
            stack.append(i)
        elif ch == "<":
            j = stack.pop()
            pair[i], pair[j] = j, i
    okp = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}
    L += ["", "## Stem pairs at the changed positions (on G. jiainica's tRNAscan-SE structure)", "",
          f"Structure: `{st}`", "",
          "| genome | stem pairs changed (G. jiainica pair / relative's bases at the same two positions) |", "|---|---|"]
    for g, hits in found.items():
        if g == "G. jiainica" or not hits:
            continue
        A, B = nw(ref, hits[0][1])
        rel, pos = {}, 0
        for x, y in zip(A, B):
            if x != "-":
                rel[pos] = y
                pos += 1
        done, items = set(), []
        for i in sorted(pair):
            j = pair[i]
            if i in done or (rel.get(i) == ref[i] and rel.get(j) == ref[j]):
                continue
            done |= {i, j}
            a, b = min(i, j), max(i, j)
            jp, rp = (ref[a], ref[b]), (rel.get(a, "-"), rel.get(b, "-"))
            kind = ("compensatory (both changed, both pair)" if rp[0] != jp[0] and rp[1] != jp[1] and jp in okp and rp in okp
                    else "pairs in both" if jp in okp and rp in okp else "pairs in G. jiainica only" if jp in okp
                    else "pairs in the relative only" if rp in okp else "pairs in neither")
            items.append(f"{a + 1}/{b + 1}: {jp[0]}-{jp[1]} / {rp[0]}-{rp[1]} ({kind})")
        L.append(f"| {g} | {'; '.join(items) if items else 'none'} |")
    (H / "trna_gly_relatives.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def genes():
    chi2_sf1 = lambda x: erfc(sqrt(max(x, 0) / 2))  # noqa: E731
    cds, notes = GG.cds_sets()
    taxa = sorted({t for g in cds.values() for t in g})
    tree_fg = (H / "genus/dnds/run_cox1_cob/M2_Gjiainica/tree.nwk").read_text().strip()
    tree0 = tree_fg.replace(" #1", "")
    t5 = {}
    for line in (H / "test5_tables.md").read_text().splitlines():
        f = [x.strip() for x in line.strip("|").split("|")]
        if len(f) == 8 and f[0] == "G. jiainica":
            t5[f[1]] = (f[6], f[7])
    D = H / "genus/dnds_per_gene"
    D.mkdir(parents=True, exist_ok=True)
    L = ["# G. jiainica's mitochondrial protein genes one by one (post hoc, after the outside review)", "",
         f"CDS, codon alignments and tree as in `dnds.md` (tree `{tree_fg}`); codeml F3x4, yeast mitochondrial code, "
         "cleandata = 1. Extra bases are the miniprot frameshifts of test 5 (G. jiainica's four sites); they are excluded "
         "from the CDS. Stops: internal stops in the CDS as used (test 5 reports none in any gene).", "",
         "| gene | codons | GGT | GGA | GGC | GGG | extra bases | internal stops | omega G. jiainica | omega background | "
         "different from background (2ΔlnL, p) | neutral rejected (2ΔlnL, p) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    table = GG.F.codon_table(3)
    for gene in GG.GENES:
        if len(cds[gene]) < len(taxa):
            L.append(f"| {gene} | missing in {set(taxa) - set(cds[gene])} | | | | | | | | | | |")
            continue
        pa, ca = GG.codon_align(cds[gene])
        Lc = len(next(iter(ca.values())))
        phy = D / f"{gene}.phy"
        phy.write_text(f" {len(taxa)} {Lc}\n" + "".join(f"{t}  {ca[t]}\n" for t in taxa))
        l0, w0, _ = GG.codeml(str(phy), tree0, 0, wd=D / f"run_{gene}/M0")
        l2, w2, _ = GG.codeml(str(phy), tree_fg, 2, wd=D / f"run_{gene}/M2")
        l1, _, _ = GG.codeml(str(phy), tree_fg, 2, fix_omega=1, omega=1, wd=D / f"run_{gene}/M2fix")
        s = cds[gene]["Gjiainica"]
        cod = [s[i:i + 3] for i in range(0, len(s) - 2, 3)]
        c = Counter(cod)
        stops = sum(1 for x in cod if x == "NNN" or table.get(x) == "*")
        fs, st = t5.get(gene, ("?", "?"))
        a, b = 2 * (l2 - l0), 2 * (l2 - l1)
        L.append(f"| {gene} | {len(cod)} | {c['GGT']} | {c['GGA']} | {c['GGC']} | {c['GGG']} | {fs} | {stops} | "
                 f"{w2[1] if len(w2) > 1 else w2[0]:.3f} | {w2[0]:.3f} | {a:.2f}, {chi2_sf1(a):.2g} | {b:.2f}, {chi2_sf1(b):.2g} |")
    # VAR1: the open frame (stop to stop, + strand) that covers most of the rps3 hit 19206-19868 (test 5), counted only
    gj = GG.mito_seqs()["jiainica"]
    hs, he = 19206 - 1, 19868
    best, cover = "", 0
    for fr in range(3):
        starts = [i for i in range(fr, len(gj) - 2, 3) if table.get(gj[i:i + 3]) == "*"]
        for x, y in zip(starts, starts[1:]):
            o = max(0, min(y, he) - max(x + 3, hs))
            if o > cover:
                cover, best = o, gj[x + 3:y]
    vc = Counter(best[i:i + 3] for i in range(0, len(best) - 2, 3))
    L.append(f"| VAR1 (open frame over the rps3 hit, + strand, stop to stop) | {len(best) // 3} | {vc['GGT']} | {vc['GGA']} | {vc['GGC']} | "
             f"{vc['GGG']} | n/a | 0 in the frame | not estimated (too divergent to align) | | | |")
    L += ["", "Notes from the CDS extraction: " + ("; ".join(n for n in notes if n.startswith("Gjiainica")) or "none for G. jiainica") + "."]
    (H / "genus/genes_per_gene.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    {"trna": trna, "genes": genes}[sys.argv[1]]()
