"""Test 3 (exploratory): nuclear genes for respiration and for cox1/cob translation in Grigorovia and two controls.

    uv run python code/explore/nuclear_genes.py fetch         # nuclear assemblies (NCBI Datasets), UniProt panel
    uv run python code/explore/nuclear_genes.py completeness  # compleasm, one reference protein per BUSCO gene
    uv run python code/explore/nuclear_genes.py pass1         # S. cerevisiae queries, miniprot default + sensitive
    uv run python code/explore/nuclear_genes.py close         # relatives' translations as queries (miniprot, DIAMOND)
    uv run python code/explore/nuclear_genes.py congener      # G. yakushimaensis proteins as queries
    uv run python code/explore/nuclear_genes.py reads         # G. jiainica reads: polish, re-call, read hits, PET122
    uv run python code/explore/nuclear_genes.py het           # heterozygosity; heterozygous indels/stops; haplotypes
    uv run python code/explore/nuclear_genes.py rates         # identity to the congener, activators vs subunits
Outputs: results/explore/grigorovia/h1h2/nuclear/ and .../completeness/. Large intermediates: data/nuclear/ (ignored).
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import re
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from math import log
from pathlib import Path
from statistics import median

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fscore as F  # noqa: E402
import nuclear_calls as N  # noqa: E402
import read_haplotypes as RH  # noqa: E402

ROOT = F.ROOT
OUT = ROOT / "results/explore/grigorovia/h1h2/nuclear"
CMP = ROOT / "results/explore/grigorovia/h1h2/completeness"
NUC = ROOT / "data/nuclear"
TMP = ROOT / "data/explore_tmp/nuclear"
ASM = {"GCA_030571715.1": "G. jiainica", "GCA_003709265.2": "G. yakushimaensis", "GCA_030462875.1": "G. humatica",
       "GCA_003708445.3": "G. transvaalensis", "GCA_964275055.1": "M. servazzii (control)", "GCF_000304475.1": "K. africana (control)"}
GRIG = ["GCA_030571715.1", "GCA_003709265.2", "GCA_030462875.1", "GCA_003708445.3"]
GJ, YAK = "GCA_030571715.1", "GCA_003709265.2"
READS = [ROOT / f"data/reads/SRR16974510/SRR16974510_{i}.fastq.gz" for i in (1, 2)]
GENES = ["PET309", "MSS51", "CBS1", "CBS2", "CBP1", "COR1", "QCR2", "CYT1", "RIP1", "QCR6", "QCR7", "QCR8", "QCR9", "QCR10",
         "COX4", "COX5A", "COX5B", "COX6", "COX7", "COX8", "COX9", "COX12", "COX13", "COX10", "COX15", "CYC1", "CYC7", "CYC3",
         "TUF1", "MEF1", "MRPL3", "MRPL4", "MRPL7", "MRPL10", "MRPL15", "MRPL16", "MRPS5", "MRPS9", "MRPS17", "MRP4",
         "CBP3", "CBP6", "CBP2", "PET111", "PET54", "PET122", "PET494", "AEP2", "ATP22"]
ACTIVATORS = ["PET309", "MSS51", "CBS2", "CBP1", "PET111", "PET494", "PET122", "ATP22", "CBP3", "CBP6"]
SENS = ["-n", "1", "-l", "4", "--outs=0.3", "-p", "0.1"]
PANEL = OUT.parent / "sc_nuclear_panel.faa"


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw).stdout


def genome(acc):
    return NUC / acc / "genome.fna"


def fetch():
    for acc in ASM:
        d = NUC / acc
        d.mkdir(parents=True, exist_ok=True)
        if genome(acc).exists():
            continue
        z = d / f"{acc}.zip"
        run(["curl", "-sfL", "-o", str(z), f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"])
        with zipfile.ZipFile(z) as zf:
            name = next(n for n in zf.namelist() if n.endswith("_genomic.fna"))
            genome(acc).write_bytes(zf.read(name))
    if not PANEL.exists():
        q = " OR ".join(f"gene_exact:{g}" for g in GENES)
        fa = run(["curl", "-sfG", "https://rest.uniprot.org/uniprotkb/stream", "--data-urlencode",
                  f"query=({q}) AND organism_id:559292 AND reviewed:true", "--data-urlencode", "format=fasta"])
        recs = [r for r in SeqIO.parse(io.StringIO(fa), "fasta") if r.description.split("GN=")[1].split()[0] in GENES]
        for r in recs:
            r.id = r.description.split("GN=")[1].split()[0] + "|" + r.id.split("|")[1]
            r.description = ""
        SeqIO.write(recs, PANEL, "fasta")
    print({a: sum(1 for _ in SeqIO.parse(genome(a), "fasta")) for a in ASM}, "panel:", sum(1 for _ in SeqIO.parse(PANEL, "fasta")))


def _nodes():
    nodes = {}
    for line in open((ROOT / "protocol/taxonomy_extract/nodes.dmp" if (ROOT / "protocol/taxonomy_extract/nodes.dmp").exists() else ROOT / "data/raw/taxdump/nodes.dmp")):
        f = [x.strip() for x in line.split("|")]
        nodes[f[0]] = f[1]
    return nodes


def completeness():
    """compleasm 0.2.10 with saccharomycetes_odb12 genes, one reference protein per gene (S. cerevisiae where the lineage
    has one, else the longest Saccharomycetaceae sequence). The full set (315,346 references) ran ~4 h per genome."""
    CMP.mkdir(parents=True, exist_ok=True)
    full = ROOT / "data/busco_lineages"
    lib = ROOT / "data/busco_lineages_sc"
    ref = lib / "saccharomycetes_odb12" / "refseq_db.faa.gz"
    if not (full / "saccharomycetes_odb12.done").exists():
        run(["uv", "run", "python", str(ROOT / "tools/compleasm/compleasm.py"), "download", "saccharomycetes", "-L", str(full)])
    if not ref.exists():
        (lib / "saccharomycetes_odb12").mkdir(parents=True, exist_ok=True)
        for f in ("dataset.cfg", "scores_cutoff"):
            (lib / "saccharomycetes_odb12" / f).write_bytes((full / "saccharomycetes_odb12" / f).read_bytes())
        h = lib / "saccharomycetes_odb12" / "hmms"
        if not h.exists():
            h.symlink_to(full / "saccharomycetes_odb12" / "hmms")
        for f in full.glob("file_versions.tsv*"):
            (lib / f.name).write_bytes(f.read_bytes())
        nodes = _nodes()

        def sac(t):
            for _ in range(60):
                if t == "4893":
                    return True
                if t not in nodes or nodes[t] == t:
                    return False
                t = nodes[t]
            return False
        seqs = defaultdict(dict)
        cur = None
        with gzip.open(full / "saccharomycetes_odb12" / "refseq_db.faa.gz", "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    hd = line[1:].strip()
                    b, tax = hd.split("_")[0], hd.split("_")[1]
                    cur = seqs[b].setdefault(tax, [hd, ""]) if tax not in seqs[b] else None
                elif cur is not None:
                    cur[1] += line.strip()
        src = Counter()
        with gzip.open(ref, "wt") as out:
            for b, d in seqs.items():
                if "559292" in d or "4932" in d:
                    pick, k = d.get("559292") or d.get("4932"), "S. cerevisiae"
                else:
                    cand = [v for t, v in d.items() if sac(t)] or list(d.values())
                    pick, k = max(cand, key=lambda v: len(v[1])), "other Saccharomycetaceae"
                src[k] += 1
                out.write(f">{pick[0]}\n{pick[1]}\n")
        (CMP / "query_sources.json").write_text(json.dumps(dict(src)))
        (lib / "saccharomycetes_odb12.done").touch()
    rows = {}
    for acc in ASM:
        od = ROOT / f"results/explore/grigorovia/h1h2/compleasm1/{acc}"
        run(["uv", "run", "python", str(ROOT / "tools/compleasm/compleasm.py"), "run", "-a", str(genome(acc)), "-o", str(od),
             "-l", "saccharomycetes", "-L", str(lib), "-t", "8"])
        (CMP / f"{acc}_summary.txt").write_text((od / "summary.txt").read_text())
        (CMP / f"{acc}_full_table.tsv").write_text((od / "saccharomycetes_odb12" / "full_table.tsv").read_text())
        rows[acc] = dict(re.findall(r"^([SDFIM]):([\d.]+)%", (od / "summary.txt").read_text(), re.M))
    # background: BUSCO genes complete in M. servazzii but missing in all four Grigorovia
    st = defaultdict(dict)
    for acc in ASM:
        for r in csv.reader(open(CMP / f"{acc}_full_table.tsv"), delimiter="\t"):
            if r[0].startswith(("Gene", "#")):
                continue
            st[r[0]][acc] = min(st[r[0]].get(acc, "Z"), r[1][0])
    ser = [b for b in st if st[b].get("GCA_964275055.1", "M") in "SD"]
    miss4 = [b for b in ser if all(st[b].get(a, "M") == "M" for a in GRIG)]
    lines = ["# Completeness (compleasm 0.2.10, saccharomycetes_odb12, one reference protein per gene)", "",
             f"Query sources: {(CMP / 'query_sources.json').read_text()}", "",
             "| assembly | species | complete single | duplicated | fragmented | missing | complete (S + D) |", "|---|---|---|---|---|---|---|"]
    for acc, r in rows.items():
        lines.append(f"| {acc} | {ASM[acc]} | {r['S']}% | {r['D']}% | {r['F']}% | {r['M']}% | {float(r['S']) + float(r['D']):.1f}% |")
    lines += ["", f"Genus-wide background: of {len(ser)} BUSCO genes complete in M. servazzii, {len(miss4)} "
              f"({len(miss4) / len(ser):.1%}) are missing from all four Grigorovia assemblies."]
    (CMP / "completeness.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def _calls(gff):
    return N.parse(gff)


def pass1():
    OUT.mkdir(parents=True, exist_ok=True)
    qs = [r.id for r in SeqIO.parse(PANEL, "fasta")]
    res = {}
    for acc in ASM:
        for mode, extra in (("default", []), ("sensitive", SENS)):
            g = OUT / f"{acc}.pass1_{mode}.gff"
            g.write_text(run(["miniprot", "-t", "8", "--gff", "-N", "2", *extra, str(genome(acc)), str(PANEL)]))
            res[(acc, mode)] = _calls(g)
    lines = ["# Pass 1: S. cerevisiae queries (miniprot 0.18, standard code). Call = default run, or the sensitive run (*) where the default found nothing", "",
             "Rules (PREREG): present >= 70% coverage, no frameshift or stop; broken = frameshift or stop; absent < 30%; else partial.", "",
             "| gene | " + " | ".join(ASM[a] for a in ASM) + " |", "|---|" + "---|" * len(ASM)]
    tsv = []
    for q in qs:
        row = []
        for acc in ASM:
            c1, r1 = N.call(res[(acc, "default")].get(q, []))
            c2, r2 = N.call(res[(acc, "sensitive")].get(q, []))
            c, r = (c1, r1) if c1 != "absent" else (c2 + "*", r2)
            row.append(f"{c} {r['cov']:.2f}" + (f" fs{r['fs']} st{r['st']}" if r['fs'] or r['st'] else "") if r else c)
            tsv.append((q.split("|")[0], acc, c, f"{r['cov']:.3f}" if r else "", r["fs"] if r else "", r["st"] if r else ""))
        lines.append(f"| {q.split('|')[0]} | " + " | ".join(row) + " |")
    (OUT / "pass1_calls.md").write_text("\n".join(lines) + "\n")
    with open(OUT / "pass1_calls.tsv", "w") as fh:
        fh.write("gene\tassembly\tcall\tcoverage\tframeshifts\tstops\n" + "".join("\t".join(map(str, t)) + "\n" for t in tsv))
    print("\n".join(lines))


def _trans(acc, faa):
    """miniprot --trans translations (sensitive settings); returns {query: (coverage, protein)} for the best hit."""
    out = run(["miniprot", "-t", "8", "--trans", "-N", "0", *SENS, str(genome(acc)), str(faa)])
    res, q = {}, None
    for l in out.splitlines():
        if l.startswith("##STA"):
            if q and q[0] not in res:
                res[q[0]] = (q[1], q[2], l.split("\t")[1].rstrip("*"))
            q = None
        elif not l.startswith("#"):
            f = l.split("\t")
            tags = {t.split(":")[0]: t.split(":")[-1] for t in f[12:]}
            q = (f[0], (int(f[3]) - int(f[2])) / int(f[1]), tags.get("fs") == "0" and tags.get("st") == "0")
    return res


def _diamond(acc_or_path, db, out, extra=()):
    q = acc_or_path if isinstance(acc_or_path, Path) else genome(acc_or_path)
    out.write_text(run(["diamond", "blastx", "--quiet", "--ultra-sensitive", "-F", "15", "--range-culling", "-k", "3", "-e", "1e-3",
                        "-q", str(q), "-d", str(db), "--outfmt", "6", "qseqid", "sseqid", "qstart", "qend", "sstart", "send",
                        "slen", "evalue", "bitscore", "pident", "btop", *extra]))


def _best_diamond(path, min_cov=0.0):
    best = {}
    for r in csv.reader(open(path), delimiter="\t"):
        g = r[1].split("|")[1]
        cov = (int(r[5]) - int(r[4]) + 1) / int(r[6])
        if cov < min_cov:
            continue
        if g not in best or float(r[8]) > best[g]["bits"]:
            best[g] = {"bits": float(r[8]), "E": float(r[7]), "cov": cov, "id": float(r[9]), "contig": r[0], "qs": int(r[2]),
                       "qe": int(r[3]), "fs": r[10].count("/") + r[10].count("\\"), "stops": r[10].count("*"), "query": r[1].split("|")[0]}
    return best


def close():
    TMP.mkdir(parents=True, exist_ok=True)
    sub = ["PET309", "CBP1", "CBS2", "CBS1", "PET494", "PET54", "PET122", "CBP2", "ATP22", "AEP2", "PET111", "COX9", "QCR10", "MRPL15", "COR1", "MRPL4"]
    recs = [r for r in SeqIO.parse(PANEL, "fasta") if r.id.split("|")[0] in sub]
    SeqIO.write(recs, TMP / "panel_sub.faa", "fasta")
    with open(OUT / "close_queries.faa", "w") as fh:
        for acc in ("GCA_964275055.1", "GCF_000304475.1", YAK, "GCA_030462875.1"):
            for q, (cov, ok, p) in _trans(acc, TMP / "panel_sub.faa").items():
                fh.write(f">{acc}|{q}|cov{cov:.2f}\n{p}\n")
        for r in recs:
            fh.write(f">SC|{r.id}\n{r.seq}\n")
    run(["diamond", "makedb", "--quiet", "--in", str(OUT / "close_queries.faa"), "-d", str(TMP / "close_queries")])
    lines = ["# Relatives' translations as queries (M. servazzii, K. africana, G. yakushimaensis, G. humatica; and S. cerevisiae)", "",
             "| gene | " + " | ".join(ASM[a] for a in GRIG) + " |", "|---|" + "---|" * 4]
    best = {}
    for acc in GRIG:
        _diamond(acc, TMP / "close_queries", OUT / f"{acc}.diamond_close.tsv")
        best[acc] = _best_diamond(OUT / f"{acc}.diamond_close.tsv")
    for g in sub:
        cells = []
        for acc in GRIG:
            b = best[acc].get(g)
            cells.append(f"cov {b['cov']:.2f}, id {b['id']:.0f}%, E {b['E']:.0e} (query {b['query'][:8]})" if b else "none")
        lines.append(f"| {g} | " + " | ".join(cells) + " |")
    (OUT / "close_queries.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def congener():
    """G. yakushimaensis proteins: miniprot translations of S. cerevisiae panel hits (coverage >= 0.7, intact), otherwise the
    longest open frame under the best DIAMOND hit of the relatives' queries. Then searched in the other three Grigorovia."""
    contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(genome(YAK), "fasta")}
    tr = _trans(YAK, PANEL)
    prot = {q.split("|")[0]: ("miniprot", p) for q, (cov, ok, p) in tr.items() if cov >= 0.7 and ok}
    dm = _best_diamond(OUT / f"{YAK}.diamond_close.tsv")
    for g, b in dm.items():
        if g in prot:
            continue
        a, e, st = (b["qs"] - 1, b["qe"], "+") if b["qs"] < b["qe"] else (b["qe"] - 1, b["qs"], "-")
        prot[g] = (f"open frame under DIAMOND hit {b['contig']}:{b['qs']}-{b['qe']}", N.longest_orf_over(contigs[b["contig"]], a, e, st))
    with open(OUT / "yak_queries.faa", "w") as fh:
        for g, (src, p) in prot.items():
            fh.write(f">YAK|{g}\n{p}\n")
    (OUT / "yak_queries_sources.tsv").write_text("gene\tsource\tlength\n" + "".join(f"{g}\t{s}\t{len(p)}\n" for g, (s, p) in prot.items()))
    run(["diamond", "makedb", "--quiet", "--in", str(OUT / "yak_queries.faa"), "-d", str(TMP / "yakq")])
    qs = [f"YAK|{g}" for g in prot]
    res = {}
    for acc in GRIG:
        g = OUT / f"{acc}.yakq.gff"
        g.write_text(run(["miniprot", "-t", "8", "--gff", "-N", "2", *SENS, str(genome(acc)), str(OUT / "yak_queries.faa")]))
        res[acc] = _calls(g)
        _diamond(acc, TMP / "yakq", OUT / f"{acc}.diamond_yakq.tsv")
    lines = ["# Congener (G. yakushimaensis) proteins as queries: miniprot calls (sensitive settings)", "",
             f"Congener proteins: {len(prot)}; none for {', '.join(x for x in GENES if x not in prot)}.", "",
             "| gene | " + " | ".join(ASM[a] for a in GRIG) + " |", "|---|" + "---|" * 4]
    tsv = []
    for q in qs:
        cells = []
        for acc in GRIG:
            c, r = N.call(res[acc].get(q, []))
            cells.append(f"{c} {r['cov']:.2f}" + (f" fs{r['fs']} st{r['st']}" if r['fs'] or r['st'] else "") if r else c)
            tsv.append((q.split("|")[1], acc, c, f"{r['cov']:.3f}" if r else "", r["fs"] if r else "", r["st"] if r else "",
                        f"{r['contig']}:{r['start']}-{r['end']}" if r else ""))
        lines.append(f"| {q.split('|')[1]} | " + " | ".join(cells) + " |")
    (OUT / "congener_query_calls.md").write_text("\n".join(lines) + "\n")
    with open(OUT / "congener_query_calls.tsv", "w") as fh:
        fh.write("gene\tassembly\tcall\tcoverage\tframeshifts\tstops\tlocus\n" + "".join("\t".join(map(str, t)) + "\n" for t in tsv))
    print("\n".join(lines))


def reads():
    """G. jiainica (same BioSample as the assembly): map reads, read consensus, re-call; DIAMOND of reads against congener
    proteins of the flagged genes; local SPAdes assembly where the gene is missing from the assembly."""
    W = NUC / GJ
    bam = W / "reads.bam"
    if not bam.exists():
        mm = subprocess.Popen(["minimap2", "-t", "8", "-ax", "sr", str(genome(GJ)), *map(str, READS)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "sort", "-@", "4", "-o", str(bam), "-"], stdin=mm.stdout, check=True, capture_output=True)
        mm.wait()
        run(["samtools", "index", str(bam)])
    flag = run(["samtools", "flagstat", str(bam)])
    pol = W / "genome.polished.fna"
    if not pol.exists():
        run(["samtools", "consensus", "-@", "4", "-a", "--show-del", "no", "--show-ins", "yes", "-f", "fasta", str(bam), "-o", str(pol)])
    g = OUT / f"{GJ}.polished.yakq.gff"
    g.write_text(run(["miniprot", "-t", "8", "--gff", "-N", "2", *SENS, str(pol), str(OUT / "yak_queries.faa")]))
    o, p = _calls(OUT / f"{GJ}.yakq.gff"), _calls(g)
    _diamond(pol, TMP / "yakq", OUT / f"{GJ}.polished.diamond_yakq.tsv")
    d0, d1 = _best_diamond(OUT / f"{GJ}.diamond_yakq.tsv", 0.5), _best_diamond(OUT / f"{GJ}.polished.diamond_yakq.tsv", 0.5)
    lines = ["# G. jiainica: assembly against read consensus", "", f"Reads mapped: {[l for l in flag.splitlines() if 'primary mapped' in l][0].strip()}", "",
             "| gene | miniprot, assembly | miniprot, read consensus | DIAMOND, assembly (cov, id, frameshifts, stops) | DIAMOND, read consensus |", "|---|---|---|---|---|"]
    flagged = []
    for q in sorted(set(o) | set(p) | {f"YAK|{x}" for x in ("COX9", "PET111", "PET122")}):
        c1, r1 = N.call(o.get(q, []))
        c2, r2 = N.call(p.get(q, []))
        gname = q.split("|")[1]
        if c1 == "present" and c2 == "present":
            continue
        flagged.append(gname)
        f = lambda c, r: f"{c} {r['cov']:.2f} fs{r['fs']} st{r['st']}" if r else c
        dd = lambda b: f"{b['cov']:.2f}, {b['id']:.1f}%, {b['fs']}, {b['stops']} ({b['contig']}:{b['qs']}-{b['qe']})" if b else "none"
        lines.append(f"| {gname} | {f(c1, r1)} | {f(c2, r2)} | {dd(d0.get(gname))} | {dd(d1.get(gname))} |")
    # read-level hits for flagged genes (plus three present genes for scale)
    keep = set(flagged) | {"COX12", "QCR9", "COX8", "MSS51", "CBS2"}
    SeqIO.write([r for r in SeqIO.parse(OUT / "yak_queries.faa", "fasta") if r.id.split("|")[1] in keep], TMP / "readcheck.faa", "fasta")
    run(["diamond", "makedb", "--quiet", "--in", str(TMP / "readcheck.faa"), "-d", str(TMP / "readcheck")])
    hits = defaultdict(list)
    ids = defaultdict(set)
    for fq in READS:
        for l in run(["diamond", "blastx", "--quiet", "--sensitive", "-k", "1", "-e", "1e-5", "-q", str(fq), "-d", str(TMP / "readcheck"),
                      "--outfmt", "6", "qseqid", "sseqid", "pident", "-p", "8"]).splitlines():
            r, s_, pid = l.split("\t")
            hits[s_.split("|")[1]].append(float(pid))
            ids[s_.split("|")[1]].add(re.sub(r"/[12]$", "", r))
    plen = {r.id.split("|")[1]: len(r.seq) for r in SeqIO.parse(TMP / "readcheck.faa", "fasta")}
    lines += ["", "## Reads hitting the congener proteins (DIAMOND --sensitive, E <= 1e-5; both read files)", "",
              "| gene | protein length | reads | mean identity |", "|---|---|---|---|"]
    lines += [f"| {g} | {plen[g]} | {len(v)} | {sum(v) / len(v):.1f}% |" for g, v in sorted(hits.items())]
    # genes missing from both the assembly and the consensus: local assembly of the read pairs hitting them
    missing = [g for g in flagged if g not in d0 and g not in d1]
    for g in missing:
        d = TMP / f"local_{g}"
        d.mkdir(parents=True, exist_ok=True)
        (d / "ids.txt").write_text("\n".join(sorted(ids[g])) + "\n")
        for i, fq in enumerate(READS, 1):
            with gzip.open(fq, "rt") as fh, open(d / f"r{i}.fq", "w") as out:
                keepr = False
                for n, line in enumerate(fh):
                    if n % 4 == 0:
                        keepr = line[1:].split()[0] in ids[g]
                    if keepr:
                        out.write(line)
        run([str(ROOT / "tools/spades_dl/SPAdes-4.3.0-Darwin/bin/spades.py"), "-1", str(d / "r1.fq"), "-2", str(d / "r2.fq"),
             "-o", str(d / "spades"), "-t", "4", "--only-assembler"])
        rr = run(["diamond", "blastx", "--quiet", "--ultra-sensitive", "-F", "15", "-e", "1e-3", "-q", str(d / "spades/contigs.fasta"),
                  "-d", str(TMP / "yakq"), "--outfmt", "6", "qseqid", "sseqid", "qlen", "sstart", "send", "slen", "pident", "evalue", "btop"])
        top = [l.split("\t") for l in rr.splitlines() if l.split("\t")[1].endswith("|" + g)]
        if top:
            t = max(top, key=lambda x: int(x[4]) - int(x[3]))
            lines += ["", f"**{g}: missing from the assembly; assembled from the {len(ids[g])} read pairs that hit it (SPAdes 4.3.0):** "
                      f"contig {t[2]} nt, aligns to congener residues {t[3]}-{t[4]} of {t[5]} at {t[6]}% identity, frameshifts "
                      f"{t[8].count('/') + t[8].count(chr(92))}, stops {t[8].count('*')}."]
    (OUT / "reads_check.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def _pileup(bam, ref, region):
    out = run(["samtools", "mpileup", "-q", "20", "-Q", "20", "-B", "-f", str(ref), "-r", region, str(bam)])
    for line in out.splitlines():
        f = line.split("\t")
        s = re.sub(r"\^.", "", f[4]).replace("$", "")
        bases, indels, i = [], Counter(), 0
        while i < len(s):
            ch = s[i]
            if ch in "+-":
                j = i + 1
                while j < len(s) and s[j].isdigit():
                    j += 1
                k = int(s[i + 1:j])
                indels[ch + s[j:j + k].upper()] += 1
                i = j + k
                continue
            bases.append(ch)
            i += 1
        yield int(f[1]), f[2].upper(), bases, indels


def het():
    W = NUC / GJ
    bam, ref = W / "reads.bam", genome(GJ)
    run(["samtools", "faidx", str(ref)])
    fai = sorted((l.split("\t") for l in open(str(ref) + ".fai")), key=lambda f: -int(f[1]))[:25]
    n = h = alt90 = 0
    afs = Counter()
    for c, L, *_ in fai:
        for pos, rb, bases, ind in _pileup(bam, ref, f"{c}:1-{L}"):
            ref_n = sum(1 for b in bases if b in ".,")
            alt = Counter(b.upper() for b in bases if b.upper() in "ACGT")
            tot = ref_n + sum(alt.values())
            if tot < 20 or len(bases) > 120:
                continue
            n += 1
            if alt:
                x, k = alt.most_common(1)[0]
                afs[round(k / tot, 1)] += 1
                if 0.25 <= k / tot <= 0.75:
                    h += 1
                elif k / tot > 0.9:
                    alt90 += 1
    span = sum(int(f[1]) for f in fai)
    lines = ["# G. jiainica nuclear heterozygosity and heterozygous variants in the panel genes", "",
             f"Genome-wide (25 largest contigs, {span:,} nt; sites with 20-120 reads at MAPQ >= 20, base quality >= 20): {n:,} sites.",
             f"- heterozygous-like (second base in 25-75% of reads): {h:,} = {1000 * h / n:.2f} per kb",
             f"- reads disagree with the assembly (> 90% another base): {alt90:,} = {1000 * alt90 / n:.3f} per kb",
             f"- alternative-base fraction, sites with any: {dict(sorted(afs.items()))}", ""]
    # per gene: heterozygous indels >= 25% of reads, heterozygous stop-gains (one-base view)
    best = _best_diamond(OUT / f"{GJ}.diamond_yakq.tsv", 0.5)
    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(ref, "fasta")}
    lines += ["## Panel genes (locus from the best DIAMOND hit of the congener protein); pileup view", "",
              "| gene | locus | heterozygous SNVs | indels in >= 25% of reads (not multiple of 3) | in-frame indels | stop-gains (single-base view) |", "|---|---|---|---|---|---|"]
    spots = []
    for g, b in sorted(best.items()):
        c, qs, qe = b["contig"], b["qs"], b["qe"]
        strand = 1 if qs < qe else -1
        a, e = min(qs, qe), max(qs, qe)
        s = seqs[c]
        fs, inf, stops, nh = [], [], [], 0
        for pos, rb, bases, ind in _pileup(bam, ref, f"{c}:{a}-{e}"):
            tot = len(bases)
            if tot < 15:
                continue
            for x, k in ind.items():
                if k / tot >= 0.25:
                    (fs if (len(x) - 1) % 3 else inf).append(f"{pos}{x}({k}/{tot})")
                    spots.append((g, c, pos))
            alt = Counter(x.upper() for x in bases if x.upper() in "ACGT")
            if alt:
                x, k = alt.most_common(1)[0]
                if 0.25 <= k / tot <= 0.75:
                    nh += 1
                    off = (pos - a) % 3 if strand == 1 else (e - pos) % 3
                    if strand == 1:
                        cs = pos - off
                        cod, i = s[cs - 1:cs + 2], pos - cs
                        xx = x
                    else:
                        ce = pos + off
                        cod, i = str(Seq(s[ce - 3:ce]).reverse_complement()), ce - pos
                        xx = str(Seq(x).reverse_complement())
                    alt_cod = cod[:i] + xx + cod[i + 1:]
                    if len(cod) == 3 and str(Seq(alt_cod).translate()) == "*" and str(Seq(cod).translate()) != "*":
                        stops.append(f"{pos}{rb}>{x}({k}/{tot})")
                        spots.append((g, c, pos))
        lines.append(f"| {g} | {c}:{a}-{e} | {nh} | {', '.join(fs) or '-'} | {', '.join(inf) or '-'} | {', '.join(stops) or '-'} |")
    # each flagged spot: read haplotypes from alignments, translated in the gene frame
    lines += ["", "## Flagged spots: reads grouped by the exact sequence they carry (alignment-based; code/explore/read_haplotypes.py)", "",
              "| gene | window | read versions (reads) | length vs assembly | translation in the gene frame | reads matching the assembly |", "|---|---|---|---|---|---|"]
    done = set()
    for g, c, pos in spots:
        key = (g, c, pos // 60)
        if key in done:
            continue
        done.add(key)
        b = best[g]
        strand = 1 if b["qs"] < b["qe"] else -1
        a, e = min(b["qs"], b["qe"]), max(b["qs"], b["qe"])
        w0, w1 = pos - 12, pos + 28
        hap = RH.haplotypes(str(bam), c, w0, w1)
        refw = seqs[c][w0 - 1:w1]
        tot = sum(hap.values())
        versions = [(sq, k) for sq, k in hap.most_common() if k >= 2]
        tr = []
        for sq, k in versions:
            if strand == 1:
                off = (3 - (w0 - a) % 3) % 3
                t = sq[off:]
            else:
                off = (3 - (e - w1) % 3) % 3
                t = str(Seq(sq).reverse_complement())[off:]
            t = t[:len(t) - len(t) % 3]
            tr.append(f"{str(Seq(t).translate())} ({k})")
        lines.append(f"| {g} | {c}:{w0}-{w1} | {' / '.join(str(k) for _, k in versions)} (of {tot}) | "
                     f"{' / '.join(f'{len(sq) - len(refw):+d}' for sq, _ in versions)} | {'; '.join(tr)} | {hap.get(refw, 0)} |")
    (OUT / "heterozygosity.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def rates():
    files = {"G. jiainica": OUT / f"{GJ}.polished.diamond_yakq.tsv", "G. humatica": OUT / "GCA_030462875.1.diamond_yakq.tsv",
             "G. transvaalensis": OUT / "GCA_003708445.3.diamond_yakq.tsv"}
    best = {sp: _best_diamond(f, 0.5) for sp, f in files.items()}
    genes = sorted(set().union(*[set(b) for b in best.values()]))
    struct = [g for g in genes if g not in ACTIVATORS]
    lines = ["# Identity to the G. yakushimaensis protein (DIAMOND, best hit covering >= 50%); G. jiainica = read-consensus assembly", "",
             "| gene | group | " + " | ".join(files) + " |", "|---|---|" + "---|" * len(files)]
    for g in ACTIVATORS + struct:
        if g not in genes:
            continue
        lines.append(f"| {g} | {'activator/assembly factor' if g in ACTIVATORS else 'subunit, heme, cyt c, translation'} | " +
                     " | ".join(f"{best[sp][g]['id']:.1f}" if g in best[sp] else "-" for sp in files) + " |")
    summ = {}
    for label, gs in (("activators/assembly factors", ACTIVATORS), ("subunits, heme, cyt c, translation", struct)):
        summ[label] = {sp: (median([best[sp][g]["id"] for g in gs if g in best[sp]]), sum(1 for g in gs if g in best[sp])) for sp in files}
    lines += ["", "| group | " + " | ".join(files) + " |", "|---|" + "---|" * len(files)]
    for label, d in summ.items():
        lines.append(f"| {label} (median %, genes) | " + " | ".join(f"{m:.1f} ({k})" for m, k in d.values()) + " |")
    a, s_ = summ["activators/assembly factors"], summ["subunits, heme, cyt c, translation"]
    ratio = lambda d, sp: -log(d[sp][0] / 100) / -log(d["G. humatica"][0] / 100)
    lines += ["", f"Distance ratio (-ln identity) G. jiainica / G. humatica: activators {ratio(a, 'G. jiainica'):.2f}, "
              f"subunits {ratio(s_, 'G. jiainica'):.2f}."]
    (OUT / "rates.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-6:]))


def rbh():
    """Reciprocal best hits (referee request): G. jiainica's and the congener's proteins for the activators and a few subunits,
    searched against the whole reviewed S. cerevisiae S288C proteome (UniProt); is the expected gene the best hit, and by
    how much over the next?"""
    TMP.mkdir(parents=True, exist_ok=True)
    prot_fa = TMP / "sc_proteome.faa"
    if not prot_fa.exists():
        prot_fa.write_text(run(["curl", "-sfG", "https://rest.uniprot.org/uniprotkb/stream", "--data-urlencode",
                                "query=organism_id:559292 AND reviewed:true", "--data-urlencode", "format=fasta"]))
    n_sc = sum(1 for _ in SeqIO.parse(prot_fa, "fasta"))
    run(["diamond", "makedb", "--quiet", "--in", str(prot_fa), "-d", str(TMP / "sc_proteome")])
    genes = ["PET309", "MSS51", "CBS2", "CBP1", "PET111", "PET122", "PET494", "ATP22", "CBP3", "CBP6", "COR1", "COX9", "QCR10", "MEF1"]
    pol = {r.id: str(r.seq).upper() for r in SeqIO.parse(NUC / GJ / "genome.polished.fna", "fasta")}
    best = _best_diamond(OUT / f"{GJ}.polished.diamond_yakq.tsv", 0.5)
    queries = []
    for g in genes:
        b = best.get(g)
        if b:
            a, e, st = (b["qs"] - 1, b["qe"], "+") if b["qs"] < b["qe"] else (b["qe"] - 1, b["qs"], "-")
            queries.append((f"jiainica|{g}", N.longest_orf_over(pol[b["contig"]], a, e, st)))
    loc = TMP / "local_PET122" / "spades" / "contigs.fasta"
    if loc.exists():
        c = max(SeqIO.parse(loc, "fasta"), key=lambda r: len(r.seq))
        sq = str(c.seq).upper()
        best_orf = max((N.longest_orf_over(sq, 0, len(sq), st) for st in "+-"), key=len)
        queries.append(("jiainica|PET122", best_orf))
    for r in SeqIO.parse(OUT / "yak_queries.faa", "fasta"):
        if r.id.split("|")[1] in genes:
            queries.append((f"yakushimaensis|{r.id.split('|')[1]}", str(r.seq)))
    (TMP / "rbh_queries.faa").write_text("".join(f">{k}\n{v}\n" for k, v in queries))
    out = run(["diamond", "blastp", "--quiet", "--ultra-sensitive", "-k", "5", "-e", "1e-3", "-q", str(TMP / "rbh_queries.faa"),
               "-d", str(TMP / "sc_proteome"), "--outfmt", "6", "qseqid", "sseqid", "stitle", "bitscore", "evalue", "pident", "qcovhsp"])
    hits = defaultdict(list)
    for l in out.splitlines():
        q, sid, title, bits, ev, pid, qcov = l.split("\t")
        gn = re.search(r"GN=(\S+)", title)
        name = gn.group(1) if gn else sid
        if name not in [h[0] for h in hits[q]]:
            hits[q].append((name, float(bits), float(ev), float(pid), float(qcov)))
    lines = [f"# Reciprocal best hits: proteins searched against the reviewed S. cerevisiae S288C proteome ({n_sc:,} proteins; DIAMOND "
             "blastp --ultra-sensitive)", "",
             "G. jiainica proteins: the longest open frame (from its first Met, so it can be shorter than the gene) under the best congener-query hit in the read-consensus assembly "
             "(PET122: the locally assembled contig). Congener proteins: as used for the gene calls.", "",
             "| query | protein length | best hit (bits, E, identity) | second-best gene (bits) | expected gene best? |", "|---|---|---|---|---|"]
    qlen = dict(queries)
    for q, _ in queries:
        h = hits.get(q, [])
        g = q.split("|")[1]
        if not h:
            lines.append(f"| {q} | {len(qlen[q])} | no hit (E <= 1e-3) | | no hit |")
            continue
        b1 = h[0]
        b2 = h[1] if len(h) > 1 else None
        lines.append(f"| {q} | {len(qlen[q])} | {b1[0]} ({b1[1]:.0f}, {b1[2]:.0e}, {b1[3]:.0f}%) | "
                     f"{f'{b2[0]} ({b2[1]:.0f})' if b2 else '-'} | {'yes' if b1[0] == g else 'no'} |")
    (OUT / "rbh.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"fetch": fetch, "completeness": completeness, "pass1": pass1, "close": close, "congener": congener, "reads": reads,
     "het": het, "rates": rates, "rbh": rbh}[sys.argv[1]]()
