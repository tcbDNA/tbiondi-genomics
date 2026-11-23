"""Exploratory funnel for the Grigorovia jiainica +1 AGG flag (PREREG log, user decision after the scan report).
Not confirmatory. Steps a-e here; f (reads) and g (literature) separately.

    uv run python code/explore/grigorovia.py a|b|c|d|e
"""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fscore as F  # noqa: E402
import flags as FL  # noqa: E402
import scan as S  # noqa: E402

ROOT = F.ROOT
ACC = "BK073075.1"
OUT = ROOT / "results/explore/grigorovia"
COMP = str.maketrans("ACGTN", "TGCAN")


def load():
    G, ref_family, ref_gene = S.load_meta()
    code = int(G[ACC]["code"])
    seq = S.read_fasta(ROOT / f"data/work/queries/code_{code}.fna", {ACC})[ACC]
    hits = [h for h in F.read_hits(ROOT / f"results/scan/fungal/hits_code{code}.tsv") if h.q == ACC]
    kept = F.cull(hits, G[ACC]["family_pool"], ref_family)
    return G, ref_family, ref_gene, code, seq, hits, kept


def ref_seqs(ids: set[str]) -> dict[str, str]:
    out, cur, buf = {}, None, []
    with open(ROOT / "data/work/ref_final.faa") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in ids:
                    out[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur in ids:
        out[cur] = "".join(buf)
    return out


def columns(h: F.Hit, seq_o: str, table: dict):
    """Aligned columns of a hit: (query_codon_start_oriented, query_codon, query_aa, subject_pos_1based) for 'aa'
    columns; frameshift ops move the cursor as in walk()."""
    cur, sp = h.ostart, h.sstart
    out = []
    for tok in F.TOKEN.findall(h.btop):
        if tok.isdigit():
            for _ in range(int(tok)):
                c = seq_o[cur:cur + 3]
                out.append((cur, c, table.get(c, "X"), sp))
                cur += 3
                sp += 1
            continue
        qa, sa = tok[0], tok[1]
        if qa in "\\/":
            cur += 1 if qa == "\\" else -1
        elif qa == "-":
            sp += 1
        elif sa == "-":
            cur += 3
        else:
            c = seq_o[cur:cur + 3]
            out.append((cur, c, table.get(c, "X"), sp))
            cur += 3
            sp += 1
    return out


def step_a():
    G, ref_family, ref_gene, code, seq, hits, kept = load()
    table = F.codon_table(code)
    ev = [e for e in csv.DictReader(open(ROOT / "results/scan/fungal/events.tsv"), delimiter="\t") if e["accession"] == ACC]
    lines = [f"# Step a: sites (Grigorovia jiainica {ACC}, code {code}; AGG = {table['AGG']}, CTG = {table['CTG']})", ""]
    rows = []
    for e in ev:
        strand = int(e["strand"])
        so = seq if strand == 1 else seq.translate(COMP)[::-1]
        ep = F.Event(ACC, e["subject"], strand, int(e["pos_oriented"]), int(e["shift"]), int(e["raw"]), 0, 0)
        tr = F.trigger(ep, so, table)
        f0 = ep.pos % 3
        base = ep.pos if ep.shift == 1 else ep.pos - 1
        right, left = tr["tpos"], tr["tpos"] - tr["window"]
        win_codons = sorted({so[t - ((t - f0) % 3): t - ((t - f0) % 3) + 3] for t in range(left, right + 1)})
        cs_r = right - ((right - f0) % 3)
        nxt = so[cs_r + 3:cs_r + 6]
        ctx_lo = base - 30
        ctx = so[ctx_lo:base] + "[" + so[base] + "]" + so[base + 1:base + 31]
        # old-frame codon split of the context, starting at a codon boundary
        s0 = ctx_lo + ((f0 - ctx_lo) % 3)
        oldframe = " ".join(so[i:i + 3] for i in range(s0, base + 31, 3))
        # shifted frame: from the base after the skip, to the next shift in the same hit or the hit end
        h = next((h for h in kept if h.s == e["subject"] and h.strand == strand and h.ostart <= ep.pos <= h.oend), None)
        stops_after, end_note = "", ""
        if h:
            later = sorted(int(o["pos_oriented"]) for o in ev if o["subject"] == e["subject"] and int(o["strand"]) == strand
                           and int(o["pos_oriented"]) > ep.pos)
            stop_at = min([h.oend] + later)
            start = base + 1 if ep.shift == 1 else base
            sh = so[start:stop_at]
            aa = "".join(table.get(sh[i:i + 3], "X") for i in range(0, len(sh) - 2, 3))
            stops_after = aa.count("*")
            end_note = (f"shifted frame read {len(aa)} codons to {'the next shift' if later and later[0] < h.oend else 'the hit end'}"
                        f" with {stops_after} stop(s); hit covers subject {h.sstart}-{h.send} of {h.slen}")
        really = "AGG" in (e["codon"], e["codon_left"]) and int(e["window"]) <= 4
        rows.append({"gene": e["gene"], "pos_fwd": e["pos_fwd"], "class": e["class"], "codons_L_R": f"{e['codon_left']}/{e['codon']}",
                     "window": e["window"], "codons_in_window": ",".join(win_codons), "next_codon": nxt,
                     "CTG_in_window": int("CTG" in win_codons), "AGG_in_window": int("AGG" in win_codons), "frame_note": end_note})
        lines += [f"## {e['gene']} @ {e['pos_fwd']} (+{e['shift']}, {e['class']})",
                  f"- placement window {tr['window']} nt; old-frame codons in window: {', '.join(win_codons)}; codon after R: {nxt}; "
                  f"CTG in window: {'yes' if 'CTG' in win_codons else 'no'}; AGG in window: {'yes' if 'AGG' in win_codons else 'no'}",
                  f"- ±30 nt (coding strand, skipped base in brackets): `{ctx}`",
                  f"- old frame: `{oldframe}`",
                  f"- {end_note}", ""]
    core_agg = [r for r in rows if r["class"] == "clean" and r["AGG_in_window"] and r["gene"] in FL.CORE]
    only_agg = [r for r in core_agg if not r["CTG_in_window"]]
    lines += [f"**Clean core-gene shifts with AGG in the placement window: {len(core_agg)}; of these without a CTG in the window: "
              f"{len(only_agg)}.**"]
    (OUT / "step_a_sites.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__" and sys.argv[1] == "a":
    step_a()


def step_b():
    """Residue the reference proteins (all hits covering the site, own family included) align to the codon read after
    the skipped base (GGN if the skipped base is the A of AGG)."""
    G, ref_family, ref_gene, code, seq, hits, kept = load()
    table = F.codon_table(code)
    ev = [e for e in csv.DictReader(open(ROOT / "results/scan/fungal/events.tsv"), delimiter="\t")
          if e["accession"] == ACC and e["class"] == "clean" and e["gene"] in FL.CORE]
    meta = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    lines = ["# Step b: residue facing the codon read after the skipped base", ""]
    gly_sites = 0
    for e in ev:
        strand = int(e["strand"])
        so = seq if strand == 1 else seq.translate(COMP)[::-1]
        base = int(e["pos_oriented"])
        after = so[base + 1:base + 4]
        res = Counter()
        detail = []
        for h in hits:
            if h.strand != strand or not (h.ostart <= base <= h.oend):
                continue
            cols = columns(h, so, table)
            # the aligned column whose query codon starts within 2 nt after the skipped base
            col = next((c for c in cols if base + 1 <= c[0] <= base + 3), None)
            if col is None:
                continue
            subj = ref_seqs({h.s})[h.s] if not hasattr(step_b, "_r") else step_b._r.get(h.s)
            aa = subj[col[3] - 1] if subj and 0 < col[3] <= len(subj) else "?"
            fam = meta.get(h.s.split("|")[0], {}).get("family_pool_name", "?")
            res[aa] += 1
            detail.append(f"{aa} ({fam}, {ref_gene.get(h.s)}, bits {h.bits:.0f})")
        top = res.most_common(1)[0][0] if res else "?"
        gly = res.get("G", 0) > sum(res.values()) / 2
        gly_sites += gly
        lines += [f"## {e['gene']} @ {e['pos_fwd']}: old-frame codons {e['codon_left']}/{e['codon']}; codon read after the skip: "
                  f"{after} = {table.get(after, 'X')}",
                  f"- reference residues at that position ({sum(res.values())} hits): {dict(res.most_common())}; majority Gly: "
                  f"{'yes' if gly else 'no'}", "- " + "; ".join(detail), ""]
    lines.append(f"**Sites where most reference proteins have glycine: {gly_sites} of {len(ev)}.**")
    (OUT / "step_b_residue.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def step_c():
    """In-frame AGG codons elsewhere in the kept alignments, and what the references have there."""
    G, ref_family, ref_gene, code, seq, hits, kept = load()
    table = F.codon_table(code)
    refs = ref_seqs({h.s for h in kept})
    n = Counter()
    res = Counter()
    lines = ["# Step c: AGG codons read in frame in the kept alignments", ""]
    allc = Counter()
    for h in kept:
        so = seq if h.strand == 1 else seq.translate(COMP)[::-1]
        g = ref_gene.get(h.s, "?")
        for cur, c, aa, sp in columns(h, so, table):
            allc[c] += 1
            if c == "AGG":
                s_aa = refs[h.s][sp - 1] if 0 < sp <= len(refs[h.s]) else "?"
                n[g] += 1
                res[s_aa] += 1
                lines.append(f"- {g}: AGG at oriented {cur} (strand {h.strand}) aligned to reference residue {s_aa}")
    tot = sum(allc.values())
    lines += ["", f"**In-frame AGG codons in kept alignments: {sum(n.values())} of {tot} aligned codons** "
                  f"(by gene: {dict(n)}; reference residues there: {dict(res)}). For comparison, in-frame AGA: {allc['AGA']}, "
                  f"CGN: {sum(v for k, v in allc.items() if k.startswith('CG'))}."]
    (OUT / "step_c_inframe_agg.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-3:]))


if __name__ == "__main__" and sys.argv[1] in ("b", "c"):
    {"b": step_b, "c": step_c}[sys.argv[1]]()


def step_f():
    """Reads from the assembly's own SRA run (SRR16974510): exact matches of the assembly context with the extra base
    (20 nt + base + 20 nt) versus without it (20 + 20), both strands."""
    import gzip
    G, ref_family, ref_gene, code, seq, hits, kept = load()
    ev = [e for e in csv.DictReader(open(ROOT / "results/scan/fungal/events.tsv"), delimiter="\t")
          if e["accession"] == ACC and e["class"] == "clean" and e["gene"] in FL.CORE]
    pats = {}
    for e in ev:
        p = int(e["pos_fwd"])
        w, wo = seq[p - 20:p + 21], seq[p - 20:p] + seq[p + 1:p + 21]
        rc = lambda x: x.translate(COMP)[::-1]  # noqa: E731
        pats[e["gene"] + "@" + e["pos_fwd"]] = {"with": (w, rc(w)), "without": (wo, rc(wo))}
    counts = {k: Counter() for k in pats}
    nreads = 0
    for fq in sorted((ROOT / "data/reads/SRR16974510").glob("*.fastq.gz")):
        with gzip.open(fq, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                nreads += 1
                for k, d in pats.items():
                    for kind, (a, b) in d.items():
                        if a in line or b in line:
                            counts[k][kind] += 1
    lines = [f"# Step f: reads (SRR16974510, {nreads} reads, both files)", "",
             "| shift | reads with the extra base | reads without it | share with |", "|---|---|---|---|"]
    ok = 0
    for k, c in counts.items():
        tot = c["with"] + c["without"]
        share = c["with"] / tot if tot else float("nan")
        ok += tot >= 10 and share >= 0.9
        lines.append(f"| {k} | {c['with']} | {c['without']} | {share:.3f} |")
    lines += ["", f"**Shifts where >= 10 spanning reads exist and >= 90% carry the extra base: {ok} of {len(counts)}.**"]
    (OUT / "step_f_reads.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__" and sys.argv[1] == "f":
    step_f()
