"""Grigorovia jiainica follow-up checks that were first run as one-off commands (exploratory, post-freeze).
Consolidated here so that every file in candidate/funnel/ can be regenerated. Run after code/explore/grigorovia.py a b c.

    uv run python code/explore/grigorovia_followup.py relatives   # closest relatives by protein identity (DIAMOND)
    uv run python code/explore/grigorovia_followup.py usage       # in-frame AGG and GGN codon use, Grigorovia + relatives
    uv run python code/explore/grigorovia_followup.py trnas       # ARAGORN and tRNAscan-SE -O on Grigorovia + relatives
    uv run python code/explore/grigorovia_followup.py gly         # tRNA-Gly loops, base 37, anticodon read support
    uv run python code/explore/grigorovia_followup.py family      # tRNA-Gly anticodons in all Saccharomycetaceae
    uv run python code/explore/grigorovia_followup.py relsites    # do the relatives carry the extra bases? (step e)
    uv run python code/explore/grigorovia_followup.py heldout     # OR703946 (G. yakushimaensis): frozen scan + sites + tRNA

Tools: DIAMOND, ARAGORN (tools/bin/aragorn), tRNAscan-SE (tools/local/bin); see scripts/setup_tools.sh.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import flags as FL  # noqa: E402
import fscore as F  # noqa: E402
import grigorovia as Gr  # noqa: E402
import scan as S  # noqa: E402
import second_assembly_check as C  # noqa: E402

ROOT = F.ROOT
OUT = Gr.OUT
ACC = Gr.ACC
COMP = Gr.COMP
TOOLS = {"PATH": f"{ROOT / 'tools/local/bin'}:{os.environ.get('PATH', '')}"}
CORE = FL.CORE


def meta():
    return {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}


def family_members(G):
    return [a for a, r in G.items() if r["family_pool_name"] == "Saccharomycetaceae"
            and r["set"] in ("refseq_main", "fungal_tpa") and int(r["length"]) > 0]


def seqs_of(G, accs):
    out = {}
    for c in {G[a]["code"] for a in accs}:
        out.update(S.read_fasta(ROOT / f"data/work/queries/code_{c}.fna", set(accs)))
    return out


def relatives():
    """Rebuild Grigorovia's core proteins from its DIAMOND alignments; rank Saccharomycetaceae genomes by blastx identity."""
    G, ref_family, ref_gene, code, seq, hits, kept = Gr.load()
    refs = Gr.ref_seqs({h.s for h in kept})
    with open(OUT / "grig_core_proteins.faa", "w") as out:
        for h in kept:
            g = ref_gene.get(h.s, "?")
            if g not in CORE:
                continue
            subj, sp, prot = refs[h.s], h.sstart, []
            for tok in F.TOKEN.findall(h.btop):
                if tok.isdigit():
                    n = int(tok)
                    prot.append(subj[sp - 1:sp - 1 + n])
                    sp += n
                    continue
                qa, sa = tok[0], tok[1]
                if qa in "\\/":
                    continue
                if qa == "-":
                    sp += 1
                    continue
                prot.append(qa)
                if sa != "-":
                    sp += 1
            out.write(f">{g}\n{''.join(prot).replace('*', 'X')}\n")
    M = meta()
    fam = [a for a in family_members(M) if a != ACC]
    by = defaultdict(list)
    for a in fam:
        by[M[a]["code"]].append(a)
    subprocess.run(["diamond", "makedb", "--in", str(OUT / "grig_core_proteins.faa"), "-d", str(OUT / "grig_core"), "--quiet"], check=True)
    res = OUT / "fam_vs_grig.tsv"
    res.unlink(missing_ok=True)
    for c, accs in by.items():
        s = seqs_of(M, accs)
        fa = OUT / f"fam_code{c}.fna"
        fa.write_text("".join(f">{a}\n{s[a]}\n" for a in accs))
        with open(res, "a") as fh:
            subprocess.run(["diamond", "blastx", "--sensitive", "-d", str(OUT / "grig_core"), "-q", str(fa), "--query-gencode", c,
                            "-F", "15", "-k", "20", "--max-hsps", "1", "-e", "1e-5", "--quiet", "--threads", "8",
                            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "bitscore"], stdout=fh, check=True)
    best = defaultdict(dict)
    for q, s_, pid, L, bits in csv.reader(open(res), delimiter="\t"):
        if s_ not in best[q] or float(bits) > best[q][s_][2]:
            best[q][s_] = (float(pid), int(L), float(bits))
    rows = sorted(((sum(p * L for p, L, _ in d.values()) / sum(L for _, L, _ in d.values()), len(d), q)
                   for q, d in best.items() if len(d) >= 4), reverse=True)
    json.dump([q for *_, q in rows[:6]], open(OUT / "closest_relatives.json", "w"))
    for w, n, q in rows[:12]:
        print(f"{w:.1f}% over {n} genes | {q} {M[q]['organism']}")


def _kept_columns(G, ref_family, ref_gene, accs, seqs):
    hits = defaultdict(list)
    for path in ("results/scan/fungal/hits_code3.tsv", "results/scan/main/hits_code3.tsv"):
        for h in F.read_hits(ROOT / path):
            if h.q in accs:
                hits[h.q].append(h)
    out = {}
    for a in accs:
        cnt = Counter()
        for h in F.cull(hits[a], G[a]["family_pool"], ref_family):
            if ref_gene.get(h.s) not in CORE:
                continue
            so = seqs[a] if h.strand == 1 else seqs[a].translate(COMP)[::-1]
            for _, c, _, _ in Gr.columns(h, so, F.codon_table(3)):
                cnt[c] += 1
        out[a] = cnt
    return out


def usage():
    G, ref_family, ref_gene = S.load_meta()
    rel = json.load(open(OUT / "closest_relatives.json"))
    accs = [ACC] + rel
    seqs = seqs_of(G, accs)
    cols = _kept_columns(G, ref_family, ref_gene, set(accs), seqs)
    lines = ["# In-frame AGG and glycine codon use (kept core alignments): Grigorovia and its six closest relatives", "",
             "| genome | aligned codons | AGG | AGG per 1,000 | AGA | AGA per 1,000 | GGA | GGT | GGG | GGC |", "|---|---|---|---|---|---|---|---|---|---|"]
    for a in accs:
        c = cols[a]
        tot = sum(c.values())
        lines.append(f"| {G[a]['organism']} ({a}) | {tot} | {c['AGG']} | {1000 * c['AGG'] / tot:.2f} | {c['AGA']} | {1000 * c['AGA'] / tot:.1f} | "
                     f"{c['GGA']} | {c['GGT']} | {c['GGG']} | {c['GGC']} |")
    (OUT / "agg_usage_relatives.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def _aragorn(fa):
    return subprocess.run([str(ROOT / "tools/bin/aragorn"), "-t", "-gcyeast", "-c", "-w", str(fa)], capture_output=True, text=True).stdout


def _trnascan(fa, out, ss=None):
    cmd = ["tRNAscan-SE", "-O", "-q", "--forceow", "--thread", "4", "-o", str(out)] + (["-f", str(ss)] if ss else []) + [str(fa)]
    subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, env={**os.environ, **TOOLS}, check=True)


def trnas():
    G = meta()
    rel = json.load(open(OUT / "closest_relatives.json"))
    accs = [ACC] + rel
    seqs = seqs_of(G, accs)
    STD = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}
    lines = ["# tRNAs: ARAGORN 1.2.41 (-t -gcyeast -c) and tRNAscan-SE 2.0.12 (-O, Infernal 1.1.5)", "",
             "| genome | ARAGORN tRNAs | ARAGORN Arg / Gly | tRNAscan-SE tRNAs | tRNAscan-SE Arg / Gly | Pro, Cys (tRNAscan-SE) |", "|---|---|---|---|---|---|"]
    for a in accs:
        fa = OUT / f"{a}.fna"
        fa.write_text(f">{a}\n{seqs[a]}\n")
        ar = re.findall(r"tRNA-(\w+)\s+c?\[[\d,]+\]\s+\d+\s+\(([acgt]{3})\)", _aragorn(fa))
        (OUT / f"aragorn_{a}.txt").write_text("\n".join(f"{x}({y})" for x, y in ar) + "\n")
        _trnascan(fa, OUT / f"trnascan_{a}.txt", OUT / f"trnascan_ss_{a}.txt")
        ts = [l.split() for l in open(OUT / f"trnascan_{a}.txt") if l.split() and l.split()[0] == a]
        iso = Counter(r[4] for r in ts)
        lines.append(f"| {G[a]['organism']} ({a}) | {len(ar)} | {','.join(y for x, y in ar if x == 'Arg')} / {','.join(y for x, y in ar if x == 'Gly')} | "
                     f"{len(ts)} | {','.join(r[5] for r in ts if r[4] == 'Arg')} / {','.join(r[5] for r in ts if r[4] == 'Gly')} | "
                     f"{iso.get('Pro', 0)}, {iso.get('Cys', 0)} |")
    lines += ["", "Anticodons 5'->3' as DNA. TCT reads AGA and, by wobble, AGG. tRNAscan-SE names the UCA-anticodon tRNA by the "
              "standard code, so 'Trp' appears missing in every genome; that is a labelling artefact."]
    (OUT / "step_d_trnas.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def gly():
    G = meta()
    accs = [ACC, "BK073182.1", "NC_004918.1", "BK072186.1"]
    lines = ["# tRNA-Gly structure (tRNAscan-SE -O secondary structure; ARAGORN anticodon)", "",
             "| genome | anticodon | score | loop (32-38) | loop length | base 37 |", "|---|---|---|---|---|---|"]
    for a in accs:
        t = (OUT / f"trnascan_ss_{a}.txt").read_text()
        for blk in t.strip().split("\n\n"):
            if "Type: Gly" not in blk:
                continue
            ac, p34 = re.search(r"Anticodon: (\w+) at (\d+)-\d+", blk).groups()
            score = re.search(r"Score: ([\d.]+)", blk).group(1)
            sq = re.search(r"Seq: (\S+)", blk).group(1)
            st = re.search(r"Str: (\S+)", blk).group(1)
            i = int(p34) - 1
            lo, hi = i, i + 2
            while lo > 0 and st[lo - 1] == ".":
                lo -= 1
            while hi < len(st) - 1 and st[hi + 1] == ".":
                hi += 1
            lines.append(f"| {G[a]['organism']} ({a}) | {ac} | {score} | {sq[lo:i]}-{sq[i:i + 3]}-{sq[i + 3:hi + 1]} | {hi - lo + 1} | {sq[i + 3]} |")
    # read support at Grigorovia tRNA-Gly anticodon position 34 (minus-strand tRNA; anticodon at 16416-16414)
    s = S.read_fasta(ROOT / "data/work/queries/code_3.fna", {ACC})[ACC]
    p = 16416 - 1
    pats = {b: s[p - 20:p] + b + s[p + 1:p + 21] for b in "ACGT"}
    rc = lambda x: x.translate(COMP)[::-1]  # noqa: E731
    cnt = Counter()
    for fq in sorted((ROOT / "data/reads/SRR16974510").glob("*.fastq.gz")):
        with gzip.open(fq, "rt") as fh:
            for k, line in enumerate(fh):
                if k % 4 == 1:
                    for b, pat in pats.items():
                        if pat in line or rc(pat) in line:
                            cnt[b] += 1
    lines += ["", "Read support at Grigorovia tRNA-Gly anticodon position 34 (SRR16974510; exact 20+1+20 nt, both strands): "
              + ", ".join(f"{rc(b)}CC {cnt[b]}" for b in "CAGT") + "."]
    (OUT / "trna_gly_structure.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def family():
    G = meta()
    fam = family_members(G)
    seqs = seqs_of(G, fam)
    ev = defaultdict(list)
    for d in ("results/scan/main", "results/scan/fungal"):
        for e in csv.DictReader(open(ROOT / d / "events.tsv"), delimiter="\t"):
            if e["accession"] in seqs:
                ev[e["accession"]].append(e)
    rows = []
    W = OUT / "fam_trna"
    W.mkdir(exist_ok=True)
    for a in fam:
        fa = W / f"{a}.fna"
        fa.write_text(f">{a}\n{seqs[a]}\n")
        gly = re.findall(r"tRNA-Gly\s+c?\[[\d,]+\]\s+\d+\s+\(([acgt]{3})\)", _aragorn(fa))
        agg = [e for e in ev[a] if int(e["shift"]) == 1 and "AGG" in (e["codon"], e["codon_left"])]
        rows.append({"accession": a, "organism": G[a]["organism"], "gly_aragorn": ",".join(gly) or "-",
                     "agg_shift_events": len(agg), "agg_clean": sum(e["class"] == "clean" for e in agg)})
    json.dump(rows, open(OUT / "fam_trna_gly.json", "w"), indent=1)
    print(Counter(r["gly_aragorn"] for r in rows), [r for r in rows if "gcc" in r["gly_aragorn"]])


def _carries(win, other, al):
    """Local alignment of a +-60 window (extra base at index 60) to another genome, both strands. Returns
    (aligned bases, identity, 'lacks'/'carries'); 'lacks' = a gap in the other genome within 3 nt of index 60."""
    best = None
    for t in (other + other[:200], (other + other[:200]).translate(COMP)[::-1]):
        a = al.align(win, t)[0]
        if best is None or a.score > best.score:
            best = a
    qb, tb = best.aligned
    cov = sum(e - s for s, e in qb)
    ident = sum(1 for (qs, qe), (ts, te) in zip(qb, tb) for i in range(qe - qs) if win[qs + i] == best.query[ts + i]) / max(1, cov)
    gaps = [(qe1, (qs2 - qe1) - (ts2 - te1)) for (qs1, qe1), (ts1, te1), (qs2, qe2), (ts2, te2) in zip(qb, tb, qb[1:], tb[1:])
            if (qs2 - qe1) > (ts2 - te1)]
    lacks = any(abs(p - 60) <= 3 and n % 3 for p, n in gaps)
    return cov, ident, ("lacks" if lacks else "carries")


SITES = (24710, 26495, 26476, 26265)   # 0-based forward positions of the four extra bases in BK073075.1


def relsites():
    G = meta()
    rel = json.load(open(OUT / "closest_relatives.json"))
    seqs = seqs_of(G, [ACC] + rel)
    q = seqs[ACC]
    al = C.aligner()
    lines = ["# Step e: do the six closest relatives carry the extra base? (+-60 nt windows, local alignment; comparable = "
             ">= 90 aligned nt at >= 70% identity; automatic call, check gaps placed in repeats by eye)", "",
             "| site (1-based) | " + " | ".join(G[a]["organism"] for a in rel) + " |", "|---|" + "---|" * len(rel)]
    for p in SITES:
        cells = []
        for a in rel:
            cov, ident, call = _carries(q[p - 60:p + 61], seqs[a], al)
            cells.append(f"{call} ({ident:.2f})" if cov >= 90 and ident >= 0.70 else f"n.c. ({ident:.2f})")
        lines.append(f"| {p + 1} | " + " | ".join(cells) + " |")
    (OUT / "step_e_relatives.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def heldout():
    from Bio import SeqIO
    rec = SeqIO.read(ROOT / "data/raw/OR703946/OR703946.1.gb", "genbank")
    yak = str(rec.seq).upper()
    G = meta()
    out = ROOT / "results/heldout/OR703946"
    out.mkdir(parents=True, exist_ok=True)
    res = S.scan(["OR703946.1"], out, 8, seqs={"OR703946.1": yak}, code_of={"OR703946.1": 3},
                 family_of={"OR703946.1": G[ACC]["family_pool"]})
    S.write_tables(res, out)
    r = res[0]
    null = FL.Null(ROOT / "results/controls/plant_null/null.json")
    stops = {c for c, x in F.codon_table(3).items() if x == "*"}
    cl = FL.clean_core(r["events"], null.split_votes)
    lines = ["# Held-out relative: Grigorovia yakushimaensis OR703946.1 (frozen pipeline, unchanged)", "",
             f"- record: {rec.description}; {len(yak)} bp; CDS features: {sum(f.type == 'CDS' for f in rec.features)}; "
             f"comment: {rec.annotations.get('comment', '')[:200].replace(chr(10), ' ')}",
             f"- core genes aligned: {r['summary']['genes_hit']}; breaks: {dict(Counter(e['class'] for e in r['events']))}; "
             f"clean core shifts: {len(cl)}; flags: {sum(t['flag'] for t in FL.genome_tests(r['events'], 3, r['summary']['gc'], null, stops))}; "
             f"clean +1 at AGG: {sum(1 for e in cl if int(e['shift']) == 1 and 'AGG' in FL.votes(e, null.cap))}"]
    q = S.read_fasta(ROOT / "data/work/queries/code_3.fna", {ACC})[ACC]
    al = C.aligner()
    for p in SITES:
        cov, ident, call = _carries(q[p - 60:p + 61], yak, al)
        lines.append(f"- site {p + 1}: {cov} nt aligned at {ident:.2f}; automatic call: {call}")
    fa = out / "OR703946.1.fna"
    fa.write_text(f">OR703946.1\n{yak}\n")
    gly = re.findall(r"tRNA-Gly\s+c?\[[\d,]+\]\s+\d+\s+\(([acgt]{3})\)", _aragorn(fa))
    lines.append(f"- tRNA-Gly (ARAGORN): {', '.join(gly)}")
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def figure():
    """Codons, code-3 reading and S. cerevisiae COB residues around the cob double site, for code/figures.py."""
    G, ref_family, ref_gene, code, seq, hits, kept = Gr.load()
    table = F.codon_table(3)
    so = seq.translate(COMP)[::-1]
    L = len(seq)
    o1, o2 = L - 1 - 26495, L - 1 - 26476
    h = next(h for h in hits if h.s.startswith("NC_001224.1|") and Gr.S.load_meta()[2].get(h.s) == "cob" and h.strand == -1)
    cols = {c[0]: c for c in Gr.columns(h, so, table)}
    sc = Gr.ref_seqs({h.s})[h.s]
    plan = [("codon", o1 - 21 + 3 * i) for i in range(7)] + [("skip", o1)] + [("codon", o1 + 1 + 3 * i) for i in range(6)] \
        + [("skip", o2)] + [("codon", o2 + 1 + 3 * i) for i in range(4)]
    items = []
    for kind, p in plan:
        if kind == "skip":
            items.append({"kind": "skip", "base": so[p]})
        else:
            c = so[p:p + 3]
            col = cols.get(p)
            items.append({"kind": "codon", "codon": c, "aa": table.get(c, "X"), "ref": sc[col[3] - 1] if col else None})
    old = so[o1:o1 + 30]
    data = {"genome": ACC, "species": "Grigorovia jiainica", "strain": "NRRL Y-48843", "gene": "cob", "sites_1based": [26496, 26477],
            "reference": f"Saccharomyces cerevisiae S288C COB ({h.s.split('|')[1]})", "items": items,
            "without_skip": [(old[i:i + 3], table.get(old[i:i + 3], "X")) for i in range(0, 30, 3)],
            "reads_with_without": {"26496": [1718, 0], "26477": [1733, 0]}}
    # anticodons (5'->3', DNA letters) of the tRNAs encoded by the Grigorovia mitogenome, from tRNAscan-SE -O (step d)
    trna = defaultdict(list)
    for line in open(OUT / f"trnascan_{ACC}.txt"):
        f = line.split()
        if f and f[0] == ACC:
            trna[f[4]].append(f[5])
    data["trna_anticodons"] = {k: sorted(set(v)) for k, v in trna.items()}
    data["trna_note"] = ("tRNAscan-SE 2.0.12 -O on BK073075.1; 'Sup' rows are named by the standard code "
                         "(TCA reads UGA = Trp in code 3)")
    json.dump(data, open(OUT / "figure_site.json", "w"), indent=1)
    print(json.dumps(data)[:300])


if __name__ == "__main__":
    {"figure": figure, "relatives": relatives, "usage": usage, "trnas": trnas, "gly": gly, "family": family, "relsites": relsites,
     "heldout": heldout}[sys.argv[1]]()
