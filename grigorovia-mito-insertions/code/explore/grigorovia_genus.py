"""Test 6 (exploratory): mitogenomes across the four Grigorovia species.
G. transvaalensis and G. humatica are assembled de novo from public WGS reads (SPAdes 4.3.0, run separately); this
script picks their mitochondrial contigs, confirms them by read mapping, and compares the four species.

    uv run python code/explore/grigorovia_genus.py pick      # mito contigs from the SPAdes assemblies
    uv run python code/explore/grigorovia_genus.py confirm   # read mapping: depth, spanning reads at key positions
    uv run python code/explore/grigorovia_genus.py compare   # frozen scan, cox1/cob exon structures, tRNA-Gly
    uv run python code/explore/grigorovia_genus.py dnds      # codon alignments, tree, codeml branch models
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fscore as F  # noqa: E402
import scan as S  # noqa: E402

ROOT = F.ROOT
OUT = ROOT / "results/explore/grigorovia/h1h2/genus"
# G. transvaalensis: strain CLIB 448 (ERR9706995). The type strain NRRL Y-17245 (SRR16974369 = SRR6475934) has
# essentially no mtDNA reads (4 of 200,000 hit mitochondrial proteins; rho0-like), so it cannot be used here.
ASM = {"humatica": ROOT / "data/assembly/humatica_SRR16974597", "transvaalensis": ROOT / "data/assembly/transvaalensis_CLIB448_ERR9706995"}
READS = {"humatica": [ROOT / f"data/reads/genus/SRR16974597/SRR16974597_{i}.fastq" for i in (1, 2)],
         "transvaalensis": [ROOT / f"data/reads/genus/ERR9706995/sub_{i}.fastq" for i in (1, 2)]}
TMP = ROOT / "data/explore_tmp"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw).stdout


def mito_seqs():
    """All four Grigorovia mitogenomes: deposited (jiainica, yakushimaensis) and assembled (the other two)."""
    seqs = {"jiainica": S.read_fasta(ROOT / "data/work/queries/code_3.fna", {"BK073075.1"})["BK073075.1"],
            "yakushimaensis": str(SeqIO.read(ROOT / "data/raw/OR703946/OR703946.1.gb", "genbank").seq).upper()}
    for sp in ASM:
        f = OUT / f"{sp}_mito.fna"
        if f.exists():
            recs = list(SeqIO.parse(f, "fasta"))
            seqs[sp] = {r.id: str(r.seq).upper() for r in recs}
    return seqs


def pick():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Mitochondrial contigs picked from the SPAdes assemblies", "",
             "Rule: a contig is mitochondrial if it has a DIAMOND blastx hit (code 3, e <= 1e-10) to an S. cerevisiae "
             "mitochondrial protein, or a minimap2 (asm20) alignment of >= 300 nt to the congener mitogenome OR703946, "
             "or GC < 25% (yeast mtDNA is AT-rich; nuclear contigs here are 30-40%), and k-mer coverage >= 3x the assembly's "
             "length-weighted median (mtDNA is multi-copy).", ""]
    for sp, d in {**ASM, "transvaalensis_NRRL_Y-17245": ROOT / "data/assembly/transvaalensis_SRR16974369"}.items():
        contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(d / "contigs.fasta", "fasta")}
        cov = {k: float(k.split("_cov_")[1]) for k in contigs}
        # length-weighted median coverage
        tot = sum(len(v) for v in contigs.values())
        acc = 0
        med = None
        for k in sorted(contigs, key=lambda k: cov[k]):
            acc += len(contigs[k])
            if acc >= tot / 2:
                med = cov[k]
                break
        dm = run(["diamond", "blastx", "--quiet", "--sensitive", "-p", "4", "--query-gencode", "3", "-F", "15", "-e", "1e-10",
                  "--range-culling", "-k", "1", "-q", str(d / "contigs.fasta"), "-d", str(TMP / "sc_mito"), "--outfmt", "6", "qseqid", "sseqid"])
        hit_genes = defaultdict(set)
        for l in dm.splitlines():
            q, s_ = l.split("\t")
            hit_genes[q].add(s_.split("_")[0])
        mm = run(["minimap2", "-x", "asm20", str(ROOT / "results/heldout/OR703946/OR703946.1.fna"), str(d / "contigs.fasta")])
        aln = Counter()
        for l in mm.splitlines():
            f = l.split("\t")
            aln[f[0]] += int(f[9])
        gc = {k: (v.count("G") + v.count("C")) / len(v) for k, v in contigs.items()}
        mt_like = lambda k: hit_genes.get(k) or aln[k] >= 300 or (gc[k] < 0.25 and len(contigs[k]) >= 300)
        keep = [k for k in contigs if mt_like(k) and cov[k] >= 3 * med]
        weak = [k for k in contigs if mt_like(k) and cov[k] < 3 * med and len(contigs[k]) >= 500]
        with open(OUT / (f"{sp}_mito.fna" if sp in ASM else f"{sp}_mito_unused.fna"), "w") as fh:
            for k in sorted(keep, key=lambda k: -len(contigs[k])):
                fh.write(f">{sp}_{k}\n{contigs[k]}\n")
        lines += [f"## G. {sp}", "", f"- assembly: {len(contigs)} contigs, {tot:,} nt; length-weighted median k-mer coverage {med:.1f}",
                  f"- mitochondrial contigs kept: {len(keep)}, total {sum(len(contigs[k]) for k in keep):,} nt", ""]
        lines += [f"  - {k}: {len(contigs[k]):,} nt, coverage {cov[k]:.0f}, genes {','.join(sorted(hit_genes.get(k, []))) or '-'}, "
                  f"aligned to OR703946 {aln[k]:,} nt, GC {gc[k]:.3f}" for k in sorted(keep, key=lambda k: -len(contigs[k]))]
        lines += ["", f"- matching but low coverage, >= 500 nt (possible nuclear copies, not used): {len(weak)}"]
        lines += [f"  - {k}: {len(contigs[k]):,} nt, coverage {cov[k]:.1f}, genes {','.join(sorted(hit_genes.get(k, []))) or '-'}" for k in weak[:10]]
        lines.append("")
    (OUT / "pick.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def confirm(species=None):
    """Map all reads back to the picked contigs; report depth and whether contig ends overlap (circular)."""
    lines = ["# Read-mapping confirmation of the assembled mitogenomes", ""]
    for sp in (species or list(ASM)):
        ref = OUT / f"{sp}_mito.fna"
        bam = TMP / f"{sp}_mito.bam"
        sam = subprocess.Popen(["minimap2", "-t", "8", "-ax", "sr", str(ref), *map(str, READS[sp])], stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "sort", "-@", "4", "-o", str(bam), "-"], stdin=sam.stdout, check=True, capture_output=True)
        sam.wait()
        run(["samtools", "index", str(bam)])
        depth = defaultdict(list)
        for l in run(["samtools", "depth", "-a", "-Q", "20", str(bam)]).splitlines():
            c, p, d = l.split("\t")
            depth[c].append(int(d))
        flag = run(["samtools", "flagstat", str(bam)]).splitlines()
        mapped = next((l for l in flag if "primary mapped" in l), "")
        lines += [f"## G. {sp}", "", f"- reads mapped (primary): {mapped.strip()}", "",
                  "| contig | length | median depth (MAPQ >= 20) | positions with depth >= 10 | min depth | ends overlap (circular)? |",
                  "|---|---|---|---|---|---|"]
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(ref, "fasta")}
        for c, ds in depth.items():
            s = seqs[c]
            ov = max((k for k in range(20, 200) if s[:k] == s[-k:]), default=0)
            srt = sorted(ds)
            lines.append(f"| {c} | {len(s):,} | {srt[len(srt) // 2]} | {sum(1 for d in ds if d >= 10) / len(ds):.3f} | {srt[0]} | "
                         f"{'yes, ' + str(ov) + ' nt' if ov else 'no'} |")
        lines.append("")
    (OUT / f"confirm_{'_'.join(species or list(ASM))}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))



# ------------------------------------------------------------------ compare
SC_PROT = {}


def sc_prot():
    if not SC_PROT:
        names = {"NP_009305.1": "cox1", "NP_009315.1": "cob", "NP_009326.1": "cox2", "NP_009328.1": "cox3",
                 "NP_009313.1": "atp6", "NP_009312.1": "atp8", "NP_009319.1": "atp9"}
        for r in SeqIO.parse(ROOT / "data/work/ref_proteins.faa", "fasta"):
            a = r.id.split("|")[-1]
            if r.id.startswith("NC_001224.1") and a in names:
                SC_PROT[names[a]] = str(r.seq)
    return SC_PROT


def cds_from_cigar(gs, tstart, cigar):
    """Coding sequence along miniprot's alignment on the oriented target `gs` (0-based start): M/D codons kept,
    F (frameshift: 1-2 target bases skipped) and N/U/V (introns) removed; G (partial codon) kept as-is."""
    out, t = [], tstart
    for n, op in re.findall(r"(\d+)([MIDFGNUV])", cigar):
        n = int(n)
        if op in "MD":
            out.append(gs[t:t + 3 * n]); t += 3 * n
        elif op == "I":
            pass
        elif op == "F":
            t += n
        elif op == "G":
            out.append(gs[t:t + n]); t += n
        elif op == "N":
            t += n
        elif op == "U":
            out.append(gs[t:t + 1]); t += n + 1
        elif op == "V":
            out.append(gs[t:t + 2]); t += n + 2
    return "".join(out)


def exon_structure(seq, gene, relaxed=False):
    """miniprot (code 3) best alignment of the S. cerevisiae protein; tries a half rotation if the gene straddles the ends.
    Returns (oriented sequence, oriented 0-based exons, rotation, strand, query start, query end) or None."""
    import tempfile
    tmp = Path(tempfile.mkdtemp(dir=TMP))
    q = sc_prot()[gene]
    (tmp / "q.faa").write_text(f">{gene}\n{q}\n")
    best = None
    for rot in (0, len(seq) // 2):
        g = seq[rot:] + seq[:rot]
        (tmp / "t.fna").write_text(f">t\n{g}\n")
        extra = ["--outs=0.3", "-p", "0.1", "-n", "1", "-l", "3"] if relaxed else []
        gff = run(["miniprot", "-T", "3", "-t", "2", "--gff", "-N", "0", *extra, str(tmp / "t.fna"), str(tmp / "q.faa")])
        mrna = [l.split("\t") for l in gff.splitlines() if "\tmRNA\t" in l]
        paf = [l.split("\t") for l in gff.splitlines() if l.startswith("##PAF")]
        if not mrna:
            continue
        cig = next(t[5:] for t in paf[0] if t.startswith("cg:Z:"))
        pstrand, pts = paf[0][5], int(paf[0][8])
        tgt = mrna[0][8].split("Target=")[1].split()
        qs, qe = int(tgt[1]), int(tgt[2])
        cds = [l.split("\t") for l in gff.splitlines() if "\tCDS\t" in l and "Parent=MP000001" in l]
        strand = cds[0][6]
        ex = sorted((int(c[3]) - 1, int(c[4])) for c in cds)
        n = len(g)
        gs = g if strand == "+" else rc(g)
        exo = [(x, y) if strand == "+" else (n - y, n - x) for x, y in ex]
        exo.sort()
        n_ = len(g)
        tstart = pts if pstrand == "+" else n_ - int(paf[0][9])
        cand = (gs, exo, rot, strand, qs, qe, cds_from_cigar(gs, tstart, cig))
        if best is None or (qe - qs) > (best[5] - best[4]):
            best = cand
        if qs <= 3 and qe >= len(q) - 3:
            break
    return best


def carries_cds(win, cds, al):
    """Extra base at win[60]. Returns 'lacks'/'carries'/'not spanned' with aligned length and identity."""
    if not cds:
        return "gene not found"
    best = None
    for t in (cds, rc(cds)):
        a = al.align(win, t)[0]
        if best is None or a.score > best[0].score:
            best = (a, t)
    a, t = best
    qb, tb = a.aligned
    cov = sum(e - s_ for s_, e in qb)
    ident = sum(1 for (qs, qe), (ts, te) in zip(qb, tb) for i in range(qe - qs) if win[qs + i] == t[ts + i]) / max(1, cov)
    if not (qb[0][0] <= 55 and qb[-1][1] >= 66):
        return f"not spanned ({cov} nt)"
    # net length difference across query 45..75 (robust to where the aligner puts a gap inside a repeat)
    m = {}
    for (qs, qe), (ts, te) in zip(qb, tb):
        for i in range(qe - qs):
            m[qs + i] = ts + i
    lo = next((x for x in range(45, 56) if x in m), None)
    hi = next((x for x in range(75, 65, -1) if x in m), None)
    if lo is None or hi is None:
        return f"not spanned ({cov} nt)"
    net = (hi - lo) - (m[hi] - m[lo])          # bases G. jiainica has in excess over the other species here
    call = "lacks" if net % 3 == 1 else ("same length" if net % 3 == 0 else f"other ({net:+d})")
    return f"{call} ({cov} nt, {ident:.2f})"


def compare():
    import grigorovia_h1h2 as H
    import grigorovia_followup as GF
    import flags as FL
    seqs = mito_seqs()
    lines = ["# Test 6: the four Grigorovia mitogenomes compared", ""]
    # 1. frozen scan on the assembled contigs
    G = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    fam = G["BK073075.1"]["family_pool"]
    null = FL.Null(ROOT / "results/controls/plant_null/null.json")
    stops = {c for c, x in F.codon_table(3).items() if x == "*"}
    lines += ["## Frozen scan pipeline (unchanged) on each mitogenome", "",
              "| species | record | length | core genes aligned | clean core shifts | +1 at AGG (clean) | flags | stop columns in core genes |", "|---|---|---|---|---|---|---|---|"]
    flat = {}
    for sp, v in seqs.items():
        d = v if isinstance(v, dict) else {("BK073075.1" if sp == "jiainica" else "OR703946.1"): v}
        flat[sp] = d
        out = OUT / f"scan_{sp}"
        res = S.scan(list(d), out, 8, seqs=d, code_of={k: 3 for k in d}, family_of={k: fam for k in d})
        S.write_tables(res, out)
        ev = [e for r in res for e in r["events"]]
        cl = FL.clean_core(ev, null.split_votes)
        genes = sorted({g for r in res for g in r["summary"]["genes_hit"].split(",") if g and not g.startswith("other")})
        agg = sum(1 for e in cl if int(e["shift"]) == 1 and "AGG" in FL.votes(e, null.cap))
        nflag = sum(t["flag"] for r in res for t in FL.genome_tests(r["events"], 3, r["summary"]["gc"], null, stops))
        st = sum(1 for r in res for x in r["stops"] if not x["gene"].startswith("other"))
        lines.append(f"| G. {sp} | {', '.join(d)[:60]} | {sum(len(x) for x in d.values()):,} | {', '.join(genes)} | {len(cl)} | {agg} | {nflag} | {st} |")
        (OUT / f"clean_shifts_{sp}.json").write_text(json.dumps([{k: e[k] for k in ('gene', 'pos_fwd', 'shift', 'codon', 'class')} for e in cl], indent=1))
    # 2. the four G. jiainica sites in each species, against each species' spliced cox1/cob CDS (introns removed)
    q = flat["jiainica"]["BK073075.1"]
    al = GF.C.aligner()
    site_gene = {24710: "cox1", 26495: "cob", 26476: "cob", 26265: "cob"}
    cds_of = {}
    for sp, d in flat.items():
        for gene in ("cox1", "cob"):
            best = None
            for k, v in d.items():
                es = exon_structure(v, gene)
                if es and (best is None or es[5] - es[4] > best[5] - best[4]):
                    best = es
            cds_of[(sp, gene)] = best[6] if best else ""
    lines += ["", "## Do the other species carry G. jiainica's extra bases?", "",
              "G. jiainica window of +-60 nt aligned locally (both strands) to each species' spliced CDS of the same gene (miniprot "
              "exon structure, introns removed). Across 15 nt either side of the extra base, 'lacks' = G. jiainica has one base more "
              "than the other species (net difference 1, mod 3); 'same length' = no net difference, i.e. the other species has a "
              "base there too; 'not spanned' = the alignment does not cover that stretch.", "",
              "| site | " + " | ".join(f"G. {sp}" for sp in seqs if sp != "jiainica") + " |", "|---|" + "---|" * (len(seqs) - 1)]
    for p in GF.SITES:
        win = q[p - 60:p + 61]
        row = []
        for sp in seqs:
            if sp == "jiainica":
                continue
            row.append(carries_cds(win, cds_of[(sp, site_gene[p])], al))
        lines.append(f"| {p + 1} | " + " | ".join(row) + " |")
    # 3. cox1/cob exon structures, junctions mapped onto G. jiainica
    genes = H.gj_genes()
    bal = H._blosum_aligner()
    lines += ["", "## cox1 and cob introns (miniprot exon structures; junctions mapped onto G. jiainica, extra bases removed)", "",
              "| species | gene | exons (oriented) | introns: position after G. jiainica coding nt (nearest site, distance) |", "|---|---|---|---|"]
    table = F.codon_table(3)
    for sp, d in flat.items():
        for gene in ("cox1", "cob"):
            hits = []
            for k, v in d.items():
                es = exon_structure(v, gene)
                if es:
                    hits.append((es[5] - es[4], k, es))
            if not hits:
                lines.append(f"| G. {sp} | {gene} | not found | |")
                continue
            _, k, (gs, exo, rot, strand, qs, qe, _c) = max(hits, key=lambda h: h[0])
            seq, fr, gprot, J = genes[gene]
            spliced = "".join(gs[x:y] for x, y in exo)
            prot = "".join(table.get(spliced[i:i + 3], "X") for i in range(0, len(spliced) - 2, 3))
            aln = bal.align(prot, gprot)[0]
            m = {}
            for (a0, a1), (b0, b1) in zip(*aln.aligned):
                for i in range(a1 - a0):
                    m[a0 + i] = b0 + i
            cum, txt = 0, []
            for j in range(len(exo) - 1):
                cum += exo[j][1] - exo[j][0]
                r, ph = cum // 3, cum % 3
                key = r - 1 if ph == 0 else r
                if key not in m:
                    txt.append("unmapped")
                    continue
                t = fr + 3 * m[key] + (2 if ph == 0 else ph - 1)
                near = min(J.items(), key=lambda kv: abs(t - kv[1]))
                dist = abs(t - near[1])
                if dist <= 40:                               # refine with the spliced nucleotide alignment
                    jj = near[1]
                    lo, hi = max(0, jj - 150), min(len(seq), jj + 151)
                    e1, e2 = exo[j][1], exo[j + 1][0]
                    i_, x1, x2, sc_, ties = H.spliced_junction(seq[lo:hi], gs[max(0, e1 - 400):e1 + 60], max(0, e1 - 400),
                                                               gs[e2 - 60:e2 + 400], e2 - 60)
                    tt = lo + i_ - 1
                    td = sorted({abs(lo + q - 1 - jj) for q in ties})
                    txt.append(f"{H.SITE_NAME[near[0]]}: {abs(tt - jj)} nt (refined; equal-score splits {td[0]}-{td[-1]} nt)")
                else:
                    txt.append(f"after {t} ({H.SITE_NAME[near[0]]}: {dist} nt)")
            lines.append(f"| G. {sp} | {gene} | {k}: {', '.join(f'{x + 1}-{y}' for x, y in exo)} (protein {qs}-{qe}) | {'; '.join(txt) or 'no intron'} |")
    # 4. tRNA-Gly
    lines += ["", "## tRNA-Gly (tRNAscan-SE -O and ARAGORN)", "", "| species | tRNAscan-SE tRNAs | tRNA-Gly anticodons (tRNAscan-SE) | ARAGORN tRNA-Gly |", "|---|---|---|---|"]
    for sp, d in flat.items():
        fa = OUT / f"{sp}_all.fna"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in d.items()))
        GF._trnascan(fa, OUT / f"trnascan_{sp}.txt", OUT / f"trnascan_ss_{sp}.txt")
        ts = [l.split() for l in open(OUT / f"trnascan_{sp}.txt") if len(l.split()) > 8 and l.split()[1].isdigit()]
        gly = [r[5] for r in ts if r[4] == "Gly"]
        ar = re.findall(r"tRNA-Gly\s+c?\[[\d,]+\]\s+\d+\s+\(([acgt]{3})\)", GF._aragorn(fa))
        lines.append(f"| G. {sp} | {len(ts)} | {', '.join(gly) or 'none'} | {', '.join(ar) or 'none'} |")
    (OUT / "compare.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))



# ------------------------------------------------------------------ dN/dS
GENES = ["cox1", "cob", "cox2", "cox3", "atp6", "atp8", "atp9"]
OUTGROUPS = {"Mservazzii": "NC_004918.1", "Abovina": "BK072186.1", "Atelluris": "BK072188.1"}
LABEL = {"jiainica": "Gjiainica", "yakushimaensis": "Gyakushima", "humatica": "Ghumatica", "transvaalensis": "Gtransvaal"}


def cds_sets():
    seqs = mito_seqs()
    out_seqs = {LABEL[sp]: (v if isinstance(v, dict) else {sp: v}) for sp, v in seqs.items()}
    og = S.read_fasta(ROOT / "data/work/queries/code_3.fna", set(OUTGROUPS.values()))
    for k, a in OUTGROUPS.items():
        out_seqs[k] = {a: og[a]}
    table = F.codon_table(3)
    cds = defaultdict(dict)
    notes = []
    for tax, d in out_seqs.items():
        for gene in GENES:
            best = None
            for k, v in d.items():
                es = exon_structure(v, gene, relaxed=(gene == "atp8")) or exon_structure(v, gene, relaxed=True)
                if es and (best is None or es[5] - es[4] > best[5] - best[4]):
                    best = es
            if not best:
                notes.append(f"{tax} {gene}: not found")
                continue
            gs, exo, rot, strand, qs, qe, sp = best
            # trim to whole codons, drop a terminal stop
            sp = sp[: len(sp) - len(sp) % 3]
            if table.get(sp[-3:], "") == "*":
                sp = sp[:-3]
            inner = [i for i in range(0, len(sp) - 2, 3) if table.get(sp[i:i + 3], "X") == "*"]
            short = [y - x for x, y in exo if y - x < 30]
            if inner:
                notes.append(f"{tax} {gene}: {len(inner)} internal stop(s) in the spliced CDS; those codons set to NNN")
                for i in inner:
                    sp = sp[:i] + "NNN" + sp[i + 3:]
            if len(exo) > 1:
                gaps = [exo[j + 1][0] - exo[j][1] for j in range(len(exo) - 1)]
                notes.append(f"{tax} {gene}: {len(exo)} segments, gaps between them {gaps} nt (introns, or skipped bases where < 5 nt)")
            cds[gene][tax] = sp
    return cds, notes


def codon_align(cds_g):
    import tempfile
    table = F.codon_table(3)
    tmp = Path(tempfile.mkdtemp(dir=TMP))
    prot = {t: "".join(table.get(s[i:i + 3], "X") if "N" not in s[i:i + 3] else "X" for i in range(0, len(s), 3)) for t, s in cds_g.items()}
    (tmp / "p.faa").write_text("".join(f">{t}\n{p}\n" for t, p in prot.items()))
    aln = run(["mafft", "--quiet", "--localpair", "--maxiterate", "1000", str(tmp / "p.faa")])
    pa = {r.id: str(r.seq) for r in SeqIO.parse(__import__("io").StringIO(aln), "fasta")}
    ca = {}
    for t, a in pa.items():
        s, i, out = cds_g[t], 0, []
        for c in a:
            if c == "-":
                out.append("---")
            else:
                out.append(s[i:i + 3])
                i += 3
        ca[t] = "".join(out)
    return pa, ca


def nj_tree(pa_concat):
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    msa = MultipleSeqAlignment([SeqRecord(Seq(s), id=t) for t, s in pa_concat.items()])
    dm = DistanceCalculator("identity").get_distance(msa)
    tree = DistanceTreeConstructor().nj(dm)
    return tree, dm


def newick_plain(tree, fg=None):
    def rec(cl):
        if cl.is_terminal():
            return cl.name + (" #1" if cl.name == fg else "")
        return "(" + ",".join(rec(c) for c in cl.clades) + ")"
    return rec(tree.root) + ";"


def codeml(phy, tree_str, model, fix_omega=0, omega=0.2, wd=None, codonfreq=2):
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "tree.nwk").write_text(tree_str + "\n")
    ctl = f"""seqfile = {phy}
treefile = tree.nwk
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = {codonfreq}
clock = 0
model = {model}
NSsites = 0
icode = 2
fix_kappa = 0
kappa = 2
fix_omega = {fix_omega}
omega = {omega}
cleandata = 1
getSE = 0
RateAncestor = 0
"""
    (wd / "codeml.ctl").write_text(ctl)
    subprocess.run(["codeml", "codeml.ctl"], cwd=wd, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    out = (wd / "out.txt").read_text()
    lnl = float(re.search(r"lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s+(-?[\d.]+)", out).group(1))
    w = re.search(r"w \(dN/dS\) for branches:\s+([\d.\s]+)", out)
    ws = [float(x) for x in w.group(1).split()] if w else [float(re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", out).group(1))]
    rows = {}
    m = re.search(r"dN & dS for each branch\s+branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS\s+([\s\S]+?)\n\n", out)
    if m:
        for l in m.group(1).strip().splitlines():
            f = l.split()
            rows[f[0]] = {"t": float(f[1]), "N": float(f[2]), "S": float(f[3]), "w": float(f[4]), "dN": float(f[5]),
                          "dS": float(f[6]), "NdN": float(f[7]), "SdS": float(f[8])}
    return lnl, ws, rows


def dnds():
    from math import erfc, sqrt
    chi2_sf1 = lambda x: erfc(sqrt(x / 2))          # chi-square survival function, 1 df
    D = OUT / "dnds"
    D.mkdir(parents=True, exist_ok=True)
    cds, notes = cds_sets()
    taxa = sorted({t for g in cds.values() for t in g})
    aligned = {}
    for gene in GENES:
        if len(cds[gene]) < len(taxa):
            notes.append(f"{gene}: missing in {set(taxa) - set(cds[gene])}; gene left out")
            continue
        aligned[gene] = codon_align(cds[gene])
    # concatenated protein alignment for the tree
    pa_cat = {t: "".join(aligned[g][0][t] for g in aligned) for t in taxa}
    tree, dm = nj_tree(pa_cat)
    sets = {"cox1+cob": ["cox1", "cob"], "cox2+cox3+atp6+atp8+atp9": [g for g in ("cox2", "cox3", "atp6", "atp8", "atp9") if g in aligned],
            "all seven": list(aligned)}
    for g in ("cox1", "cob"):
        sets[g] = [g]
    order = taxa
    lines = ["# Branch dN/dS (codeml, PAML 4.10.10; yeast mitochondrial code, icode = 2; F3x4; cleandata = 1)", "",
             f"Taxa: {', '.join(order)}. CDS from miniprot exon structures (S. cerevisiae proteins); G. jiainica's four extra bases "
             "fall in the gaps between miniprot segments and are therefore excluded. Codon alignments: MAFFT L-INS-i on proteins, "
             "back-translated. Tree: neighbour joining on concatenated protein p-distances (unrooted), "
             f"`{newick_plain(tree)}`.", "",
             "Tests per Grigorovia lineage (its terminal branch as foreground):",
             "- **A: different from background** — two-ratio model (model = 2) against one ratio (M0), 1 df.",
             "- **B: still constrained** — foreground omega free against foreground omega fixed at 1 (neutral), 1 df. "
             "Rejecting omega = 1 with omega < 1 means purifying selection on that branch.",
             "- **Power**: expected substitutions on the branch, N·dN and S·dS (two-ratio model). With few synonymous "
             "changes, omega on that branch cannot be estimated well and neither test can reject.", ""]
    res_json = {}
    for sname, gl in sets.items():
        if not gl:
            continue
        ca = {t: "".join(aligned[g][1][t] for g in gl) for t in order}
        L = len(next(iter(ca.values())))
        phy = D / f"{sname.replace('+', '_').replace(' ', '_')}.phy"
        phy.write_text(f" {len(order)} {L}\n" + "".join(f"{t}  {ca[t]}\n" for t in order))
        base = D / f"run_{phy.stem}"
        l0, w0, rows0 = codeml(str(phy), newick_plain(tree), 0, wd=base / "M0")
        lines += [f"## {sname} ({L // 3} codons before gap removal)", "", f"One-ratio M0: omega = {w0[0]:.3f}, lnL = {l0:.2f}", "",
                  "| lineage | omega background | omega branch | A: 2ΔlnL, p | B (omega = 1): 2ΔlnL, p | N·dN | S·dS | dS (branch) |",
                  "|---|---|---|---|---|---|---|---|"]
        for sp in ("jiainica", "yakushimaensis", "humatica", "transvaalensis"):
            fg = LABEL[sp]
            t2 = newick_plain(tree, fg)
            l2, w2, rows2 = codeml(str(phy), t2, 2, wd=base / f"M2_{fg}")
            l1, w1, _ = codeml(str(phy), t2, 2, fix_omega=1, omega=1, wd=base / f"M2fix_{fg}")
            idx = order.index(fg) + 1
            br = next((v for k, v in rows2.items() if k.split("..")[1] == str(idx)), None)
            A = 2 * (l2 - l0)
            B = 2 * (l2 - l1)
            pA, pB = chi2_sf1(max(A, 0)), chi2_sf1(max(B, 0))
            lines.append(f"| G. {sp} | {w2[0]:.3f} | {w2[1]:.3f} | {A:.2f}, {pA:.2g} | {B:.2f}, {pB:.2g} | "
                         f"{br['NdN']:.1f} | {br['SdS']:.1f} | {br['dS']:.3f} |" if br else f"| G. {sp} | {w2} | | | | | | |")
            res_json[f"{sname}|{sp}"] = {"w_bg": w2[0], "w_fg": w2[1], "A": A, "pA": pA, "B": B, "pB": pB, "branch": br, "lnL0": l0, "lnL2": l2, "lnL2fix": l1}
        lines.append("")
    lines += ["## Notes", ""] + [f"- {n}" for n in notes]
    (OUT / "dnds.md").write_text("\n".join(lines) + "\n")
    (OUT / "dnds.json").write_text(json.dumps(res_json, indent=1))
    print("\n".join(lines))



# ------------------------------------------------------------------ assembly, mtDNA reads, baiting
RUNS = {"G. jiainica NRRL Y-48843 (SRR16974510)": ROOT / "data/reads/SRR16974510/SRR16974510_1.fastq.gz",
        "G. humatica NRRL Y-48839 (SRR16974597)": ROOT / "data/reads/genus/SRR16974597/SRR16974597_1.fastq",
        "G. transvaalensis NRRL Y-17245, type strain (SRR16974369 = SRR6475934)": ROOT / "data/reads/genus/SRR16974369/SRR16974369_1.fastq",
        "G. transvaalensis CLIB 448 (ERR9706995)": ROOT / "data/reads/genus/ERR9706995/ERR9706995_1.fastq"}
SPADES = ROOT / "tools/spades_dl/SPAdes-4.3.0-Darwin/bin/spades.py"


def assemble():
    """Whole-genome SPAdes assemblies (skipped if present). Commands as run on 2026-11-22."""
    R = ROOT / "data/reads/genus"
    jobs = [("humatica_SRR16974597", [R / "SRR16974597/SRR16974597_1.fastq", R / "SRR16974597/SRR16974597_2.fastq"], "21,33,55,77", None),
            ("transvaalensis_SRR16974369", [R / "SRR16974369/sub_1.fastq", R / "SRR16974369/sub_2.fastq"], "21,33,55,77", R / "SRR16974369/SRR16974369"),
            ("transvaalensis_CLIB448_ERR9706995", [R / "ERR9706995/sub_1.fastq", R / "ERR9706995/sub_2.fastq"], "21,33,55,77,99", R / "ERR9706995/ERR9706995")]
    for name, (r1, r2), ks, full in jobs:
        d = ROOT / "data/assembly" / name
        if (d / "contigs.fasta").exists():
            continue
        if full is not None and not r1.exists():          # first 2,000,000 pairs
            for i, f in ((1, r1), (2, r2)):
                with open(f"{full}_{i}.fastq") as fh, open(f, "w") as out:
                    for n, line in enumerate(fh):
                        if n >= 8_000_000:
                            break
                        out.write(line)
        run(["python3", str(SPADES), "-1", str(r1), "-2", str(r2), "-o", str(d), "-t", "10", "-m", "28", "--only-assembler", "-k", ks])


def mtreads():
    """Share of the first 200,000 reads (read 1) that hit an S. cerevisiae mitochondrial protein (DIAMOND --sensitive, code 3)."""
    lines = ["# mtDNA reads per run (first 200,000 reads of read 1; DIAMOND blastx --sensitive, code 3, E <= 1e-5, "
             "against the eight S. cerevisiae mitochondrial proteins)", "", "| run | mean read length | reads hitting | share |", "|---|---|---|---|"]
    for label, fq in RUNS.items():
        sub = TMP / f"mt200k_{abs(hash(label)) % 10 ** 8}.fq"
        op = __import__("gzip").open if str(fq).endswith(".gz") else open
        with op(fq, "rt") as fh, open(sub, "w") as out:
            ln = []
            for n, line in enumerate(fh):
                if n >= 800_000:
                    break
                out.write(line)
                if n % 4 == 1:
                    ln.append(len(line.strip()))
        hits = run(["diamond", "blastx", "--quiet", "--sensitive", "-p", "8", "--query-gencode", "3", "-e", "1e-5", "-k", "1",
                    "-q", str(sub), "-d", str(TMP / "sc_mito"), "--outfmt", "6", "qseqid"])
        n = len(set(hits.split()))
        lines.append(f"| {label} | {sum(ln) / len(ln):.0f} | {n:,} | {n / 200000:.3%} |")
    (OUT / "mtreads.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def bait(rounds=2):
    """G. transvaalensis CLIB 448: the whole-genome assembly's mitochondrial contigs were fragmented, so read pairs matching
    them (minimap2 -x sr, >= 90% identity over >= 80 nt, either mate) are pulled from all 6.0 million pairs and assembled
    de novo (SPAdes with error correction); contigs >= 500 nt with GC < 30% seed the next round."""
    R = ROOT / "data/reads/genus/ERR9706995"
    T = TMP / "bait_tra"
    T.mkdir(parents=True, exist_ok=True)
    whole = OUT / "transvaalensis_mito.fna"
    (OUT / "transvaalensis_mito_wholeassembly.fna").write_text(whole.read_text())
    (T / "seed0.fna").write_text(whole.read_text())
    lines = ["# G. transvaalensis CLIB 448: mitochondrial read baiting and reassembly", "",
             f"Seed 0: the whole-genome assembly's mitochondrial contigs ({sum(1 for _ in SeqIO.parse(T / 'seed0.fna', 'fasta'))} contigs).", ""]
    for rd in range(1, rounds + 1):
        paf = run(["minimap2", "-t", "10", "-x", "sr", "--secondary=no", str(T / f"seed{rd - 1}.fna"), str(R / "ERR9706995_1.fastq"), str(R / "ERR9706995_2.fastq")])
        ids = set()
        for l in paf.splitlines():
            f = l.split("\t")
            if int(f[10]) >= 80 and int(f[9]) / int(f[10]) >= 0.9:
                ids.add(re.sub(r"/[12]$", "", f[0]))
        for i in (1, 2):
            with open(R / f"ERR9706995_{i}.fastq") as fh, open(T / f"r{rd}_{i}.fq", "w") as out:
                keep = False
                for n, line in enumerate(fh):
                    if n % 4 == 0:
                        keep = line[1:].split()[0] in ids
                    if keep:
                        out.write(line)
        d = T / f"spades{rd}"
        run(["python3", str(SPADES), "-1", str(T / f"r{rd}_1.fq"), "-2", str(T / f"r{rd}_2.fq"), "-o", str(d), "-t", "10", "-m", "28",
             "-k", "21,33,55,77,99,127"])
        recs = [r for r in SeqIO.parse(d / "contigs.fasta", "fasta") if len(r.seq) >= 500]
        mt = [r for r in recs if (str(r.seq).count("G") + str(r.seq).count("C")) / len(r.seq) < 0.30]
        SeqIO.write(mt, T / f"seed{rd}.fna", "fasta")
        lines.append(f"- round {rd}: {len(ids):,} read pairs; contigs >= 500 nt: {len(recs)}; GC < 30%: {len(mt)}, total "
                     f"{sum(len(r.seq) for r in mt):,} nt, largest {max(len(r.seq) for r in mt):,} nt")
    with open(whole, "w") as fh:
        for r in SeqIO.parse(T / f"seed{rounds}.fna", "fasta"):
            fh.write(f">transvaalensis_bait_{r.id}\n{r.seq}\n")
    lines.append(f"\nFinal G. transvaalensis mitochondrial contigs: round {rounds} ({whole.name}).")
    (OUT / "bait.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def genes_in_contigs():
    """Completeness of the seven core genes in each assembled mitogenome (miniprot exon structures)."""
    lines = ["# Core genes in the assembled mitogenomes (miniprot, S. cerevisiae proteins, code 3)", "",
             "| species | gene | contig | protein aligned | exons | CDS length |", "|---|---|---|---|---|---|"]
    for sp in ASM:
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(OUT / f"{sp}_mito.fna", "fasta")}
        for gene in GENES:
            best = None
            for k, v in seqs.items():
                es = exon_structure(v, gene) or exon_structure(v, gene, relaxed=True)
                if es and (best is None or es[5] - es[4] > best[1][5] - best[1][4]):
                    best = (k, es)
            if best:
                k, es = best
                lines.append(f"| G. {sp} | {gene} | {k} | {es[4]}-{es[5]} of {len(sc_prot()[gene])} | {len(es[1])} | {len(es[6])} |")
            else:
                lines.append(f"| G. {sp} | {gene} | not found | | | |")
    (OUT / "genes_in_contigs.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def profile():
    """Profile-likelihood 95% upper bound for omega on the G. jiainica branch (two-ratio model, foreground omega fixed)."""
    D = OUT / "dnds"
    res = json.load(open(OUT / "dnds.json"))
    tree = (D / "run_cox1_cob" / "M2_Gjiainica" / "tree.nwk").read_text().strip()
    out, lines = {}, ["# Upper bound on omega for the G. jiainica branch (profile likelihood, 2 dlnL = 3.84)", "",
                      "Largest neutral share = (upper bound - background) / (1 - background): the share of the branch that could have "
                      "evolved at omega = 1 if the rest evolved at the background omega.", "",
                      "| gene set | omega branch | omega background | profile (omega fixed: 2 dlnL) | 95% upper bound | largest neutral share |", "|---|---|---|---|---|---|"]
    for sname, stem in (("cox1+cob", "cox1_cob"), ("cox1", "cox1"), ("cob", "cob"), ("all seven", "all_seven"),
                        ("cox2+cox3+atp6+atp8+atp9", "cox2_cox3_atp6_atp8_atp9")):
        r = res[f"{sname}|jiainica"]
        l2, wbg, wfg = r["lnL2"], r["w_bg"], r["w_fg"]
        prof = []
        for w in (0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.75):
            if w <= wfg:
                continue
            l, ws, _ = codeml(str(D / f"{stem}.phy"), tree, 2, fix_omega=1, omega=w, wd=D / f"profile_{stem}" / f"w{w}")
            prof.append((w, 2 * (l2 - l)))
            if prof[-1][1] > 3.84:
                break
        up, prev = None, (wfg, 0.0)
        for w, d in prof:
            if d >= 3.84:
                up = prev[0] + (3.84 - prev[1]) * (w - prev[0]) / (d - prev[1])
                break
            prev = (w, d)
        share = max(0.0, (up - wbg) / (1 - wbg)) if up is not None else None
        out[sname] = {"w_fg": wfg, "w_bg": wbg, "profile": prof, "upper95": up, "max_neutral_share": share}
        lines.append(f"| {sname} | {wfg:.3f} | {wbg:.3f} | {', '.join(f'{w}: {d:.2f}' for w, d in prof)} | "
                     f"{up:.3f} | {share:.1%} |" if up is not None else f"| {sname} | {wfg:.3f} | {wbg:.3f} | {prof} | not reached | |")
    (OUT / "dnds_profile_jiainica.json").write_text(json.dumps(out, indent=1))
    (OUT / "dnds_profile.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def refine():
    """cob intron next to the TGG site in the two congeners that have it: junction refined by the spliced alignment."""
    import grigorovia_h1h2 as H
    genes = H.gj_genes()
    seq, fr, gprot, J = genes["cob"]
    sp = mito_seqs()
    lines = ["# cob intron junctions near the G. jiainica cob sites, refined (spliced nucleotide alignment of G. jiainica's coding window)", "",
             "| species | site | distance | equal-score splits | junction sequence |", "|---|---|---|---|---|"]
    targets = {"yakushimaensis": sp["yakushimaensis"]}
    targets.update({f"humatica ({k})": v for k, v in sp["humatica"].items()})
    targets.update({f"transvaalensis ({k})": v for k, v in sp["transvaalensis"].items()})
    for name, s in targets.items():
        es = exon_structure(s, "cob")
        if not es or len(es[1]) < 2:
            continue
        gs, exo = es[0], es[1]
        for j in range(len(exo) - 1):
            e1, e2 = exo[j][1], exo[j + 1][0]
            for e, jj in J.items():
                lo, hi = max(0, jj - 150), min(len(seq), jj + 151)
                i, x1, x2, sc, ties = H.spliced_junction(seq[lo:hi], gs[max(0, e1 - 400):e1 + 60], max(0, e1 - 400), gs[e2 - 60:e2 + 400], e2 - 60)
                t = lo + i - 1
                td = sorted({abs(lo + k - 1 - jj) for k in ties})
                if abs(t - jj) <= 30:
                    lines.append(f"| G. {name.split(' (')[0]} | {H.SITE_NAME[e]} | {abs(t - jj)} nt {'downstream' if t > jj else 'upstream'} | "
                                 f"{td[0]}-{td[-1]} nt | `{gs[x1 - 9:x1]}|{gs[x1:x1 + 6]}...{gs[x2 - 6:x2]}|{gs[x2:x2 + 9]}` |")
    lines.append("\nG. transvaalensis cob has no intron (see compare.md), so it has no row.")
    (OUT / "refine_cob_junction.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))



def downloads():
    """Sizes of everything downloaded for the H1/H2 tests (files on disk)."""
    items = [("reads", "G. humatica SRR16974597 (.sra, NCBI AWS mirror)", ROOT / "data/reads/genus/SRR16974597/SRR16974597.sra"),
             ("reads", "G. transvaalensis NRRL Y-17245 SRR16974369 (.sra)", ROOT / "data/reads/genus/SRR16974369/SRR16974369.sra"),
             ("reads", "G. transvaalensis CLIB 448 ERR9706995 (.sra)", ROOT / "data/reads/genus/ERR9706995/ERR9706995.sra")]
    for acc in ("GCA_030571715.1", "GCA_003709265.2", "GCA_030462875.1", "GCA_003708445.3", "GCA_964275055.1", "GCF_000304475.1"):
        items.append(("nuclear assembly", f"{acc} (NCBI Datasets zip)", ROOT / f"data/nuclear/{acc}/{acc}.zip"))
    for f in sorted((ROOT / "data/busco_lineages").glob("*.tar.gz")):
        items.append(("BUSCO lineage", f.name, f))
    for f in sorted((ROOT / "data/rfam").glob("*.cm")):
        items.append(("Rfam model", f.name, f))
    lines = ["# Downloads for the H1/H2 follow-up (sizes of the files on disk)", "", "| kind | item | bytes |", "|---|---|---|"]
    tot = {}
    for kind, label, f in items:
        n = f.stat().st_size
        tot[kind] = tot.get(kind, 0) + n
        lines.append(f"| {kind} | {label} | {n:,} |")
    lines += ["", "| kind | total |", "|---|---|"] + [f"| {k} | {v / 1e9:.2f} GB |" for k, v in tot.items()]
    lines.append(f"| all | {sum(tot.values()) / 1e9:.2f} GB |")
    lines.append("\nAlso used, already on disk: SRR16974510 (G. jiainica reads, 0.28 GB, downloaded for step f). The UniProt panel "
                 "(49 proteins) and NC_001224.1 (GenBank) are small. A first ENA download of ERR9706995 was abandoned after 56 MB "
                 "for the faster AWS mirror.")
    (OUT / "downloads.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))



def _prune_phy(src, dst, keep):
    lines = src.read_text().splitlines()
    rows = [l for l in lines[1:] if l.split()[0] in keep]
    L = len(rows[0].split()[1])
    dst.write_text(f" {len(rows)} {L}\n" + "\n".join(rows) + "\n")


def dnds_sensitivity():
    """Robustness of the G. jiainica branch result (referee request): codon-frequency models F1x4 (1), F3x4 (2), F61 (3),
    and F3x4 without the outgroups; for cox1+cob and cob alone. Same alignments and tree topology."""
    from math import erfc, sqrt
    chi = lambda x: erfc(sqrt(max(x, 0) / 2))
    D = OUT / "dnds"
    tree_full = (D / "run_cox1_cob" / "M2_Gjiainica" / "tree.nwk").read_text().strip()
    ing = ["Gjiainica", "Gyakushima", "Ghumatica", "Gtransvaal"]
    tree_in = "(Gyakushima,(Gtransvaal,Ghumatica),Gjiainica #1);"          # unrooted: trifurcation at the base
    settings = [("F1x4", 1, None), ("F3x4", 2, None), ("F61", 3, None), ("F3x4, no outgroups", 2, ing)]
    lines = ["# Sensitivity of the G. jiainica branch test to the codon model and the outgroups", "",
             "Two-ratio branch model (G. jiainica foreground) against one ratio (test A) and against foreground omega = 1 (test B); "
             "95% upper bound on the foreground omega from the profile likelihood; largest neutral share = (bound - background) / "
             "(1 - background), a share of synonymous change along the branch, not of time.", "",
             "| gene set | setting | lnL (two-ratio) | omega branch | omega background | A: p | B (omega = 1): p | dS branch | N·dS | 95% upper bound | largest neutral share |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    out = {}
    for stem, sname in (("cox1_cob", "cox1+cob"), ("cob", "cob"), ("cox1", "cox1")):
        for label, cf, keep in settings:
            phy = D / f"{stem}.phy"
            tree_fg, tree_0 = tree_full, tree_full.replace(" #1", "")
            if keep:
                phy2 = D / f"{stem}_ingroup.phy"
                _prune_phy(phy, phy2, set(keep))
                phy, tree_fg, tree_0 = phy2, tree_in, tree_in.replace(" #1", "")
            wd = D / f"sens_{stem}_{label.replace(' ', '_').replace(',', '')}"
            l0, w0, _ = codeml(str(phy), tree_0, 0, wd=wd / "M0", codonfreq=cf)
            l2, w2, rows2 = codeml(str(phy), tree_fg, 2, wd=wd / "M2", codonfreq=cf)
            l1, _, _ = codeml(str(phy), tree_fg, 2, fix_omega=1, omega=1, wd=wd / "M2fix", codonfreq=cf)
            idx = (keep or ["Abovina", "Atelluris", "Gjiainica", "Ghumatica", "Gtransvaal", "Gyakushima", "Mservazzii"]).index("Gjiainica") + 1
            order = [l.split()[0] for l in phy.read_text().splitlines()[1:]]
            idx = order.index("Gjiainica") + 1
            br = next((v for k, v in rows2.items() if k.split("..")[1] == str(idx)), None)
            prof, up, prev = [], None, (w2[1], 0.0)
            for w in (0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.5):
                if w <= w2[1]:
                    continue
                l, _, _ = codeml(str(phy), tree_fg, 2, fix_omega=1, omega=w, wd=wd / f"w{w}", codonfreq=cf)
                d = 2 * (l2 - l)
                prof.append((w, d))
                if d >= 3.84:
                    up = prev[0] + (3.84 - prev[1]) * (w - prev[0]) / (d - prev[1])
                    break
                prev = (w, d)
            share = max(0.0, (up - w2[0]) / (1 - w2[0])) if up is not None else None
            out[f"{sname}|{label}"] = {"lnL2": l2, "w_fg": w2[1], "w_bg": w2[0], "pA": chi(2 * (l2 - l0)), "pB": chi(2 * (l2 - l1)),
                                       "branch": br, "upper95": up, "share": share, "profile": prof}
            lines.append(f"| {sname} | {label} | {l2:.2f} | {w2[1]:.3f} | {w2[0]:.3f} | {chi(2 * (l2 - l0)):.2g} | {chi(2 * (l2 - l1)):.2g} | "
                         f"{br['dS']:.3f} | {br['N'] * br['dS']:.0f} | " + (f"{up:.3f} | {share:.1%} |" if up is not None else "not reached | |"))
    (OUT / "dnds_sensitivity.json").write_text(json.dumps(out, indent=1))
    (OUT / "dnds_sensitivity.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def ggc():
    """GGC codons in G. jiainica's seven core genes (extra bases removed) and what every taxon has at the four sites.
    Also GC at third codon positions."""
    import grigorovia_h1h2 as H
    D = OUT / "dnds"
    table = F.codon_table(3)
    genes = H.gj_genes()
    lines = ["# GGC codons in G. jiainica's core genes, and the residues at the four sites across the genus (post hoc)", ""]
    # G. jiainica coding sequences: cox1/cob from gj_genes (extra bases removed), others from the codon alignment
    aln = {}
    for g in GENES:
        rows = [l.split() for l in (D / f"{g}.phy").read_text().splitlines()[1:]] if (D / f"{g}.phy").exists() else None
        if rows is None:
            # single-gene alignments exist only for cox1/cob; take others from the concatenation order
            continue
        aln[g] = {t: s for t, s in rows}
    cat = {t: s for t, s in (l.split() for l in (D / "all_seven.phy").read_text().splitlines()[1:])}
    gj = cat["Gjiainica"].replace("-", "")
    cods = [gj[i:i + 3] for i in range(0, len(gj) - 2, 3)]
    n_ggc = sum(1 for c in cods if c == "GGC")
    gc3 = sum(1 for c in cods if c[2] in "GC") / len(cods)
    lines += [f"- G. jiainica, seven core genes (cox1, cob, cox2, cox3, atp6, atp8, atp9; extra bases removed): {len(cods):,} codons, "
              f"**{n_ggc} GGC**; GC at third codon positions {gc3:.1%}.",
              f"- Gly codons: " + ", ".join(f"{c} {sum(1 for x in cods if x == c)}" for c in ("GGT", "GGA", "GGG", "GGC")) + "."]
    # where the GGC codons are: codon right after each extra base in the corrected cox1/cob
    site_rows = []
    for gene in ("cox1", "cob"):
        seq, fr, prot, J = genes[gene]
        for e, jj in sorted(J.items()):
            # codon starting right after the extra base (jj is the corrected index the base follows)
            c = seq[jj + 1:jj + 4]
            # find this codon in the aligned G. jiainica sequence of this gene: match the 30-nt context
            ctx = seq[jj + 1 - 15:jj + 1 + 18]
            a = aln[gene]
            gj_al = a["Gjiainica"]
            ungapped = gj_al.replace("-", "")
            k = ungapped.find(ctx)
            col = None
            if k >= 0:
                target = k + 15                           # ungapped index of the codon's first base
                n = -1
                for ci in range(0, len(gj_al), 3):
                    cod = gj_al[ci:ci + 3]
                    if cod != "---":
                        n += 3
                    if n >= target:
                        col = ci
                        break
            row = {"site": H.SITE_NAME[e], "gj_codon_after_base": c}
            if col is not None:
                for t, s_ in a.items():
                    cc, nx = s_[col:col + 3], s_[col + 3:col + 6]
                    f_ = lambda x: f"{x} ({table.get(x, '-') if '-' not in x else 'gap'})"
                    row[t] = f"{f_(cc)}, then {f_(nx)}"
            site_rows.append(row)
    taxa = ["Gjiainica", "Gyakushima", "Ghumatica", "Gtransvaal", "Mservazzii", "Abovina", "Atelluris"]
    lines += ["", "## The codon right after each extra base (scan placement) and the next codon, in every taxon (codon alignments, `codon_alignments/`)", "",
              "At cob 26266 the scan places the extra base as the A of ACT, so GGC is the second codon; placing it as the A of AGG "
              "instead makes GGC the first, with the same protein. Gaps at cob 26477 in G. yakushimaensis and G. humatica are codons "
              "lost at their intron junction in the CDS extraction (removed by codeml's cleandata).", "",
              "| site | " + " | ".join(taxa) + " |", "|---|" + "---|" * len(taxa)]
    for r in site_rows:
        lines.append(f"| {r['site']} | " + " | ".join(r.get(t, "?") for t in taxa) + " |")
    nsite = sum(1 for r in site_rows if "GGC" in r.get("Gjiainica", ""))
    lines += ["", f"GGC codons in G. jiainica's seven core genes: {n_ggc}; of them at the extra-base sites (first or second codon "
              f"after the base): {nsite}."]
    (OUT / "ggc_sites.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"assemble": assemble, "mtreads": mtreads, "pick": pick, "bait": bait, "confirm": confirm, "genes": genes_in_contigs,
     "compare": compare, "dnds": dnds, "profile": profile, "refine": refine, "downloads": downloads, "dnds_sensitivity": dnds_sensitivity, "ggc": ggc}[sys.argv[1]]()
