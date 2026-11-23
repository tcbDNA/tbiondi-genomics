"""Descriptive context for every scan flag, and the scoring of the six predictions.

ADDED AFTER THE FREEZE (PREREG log entry "added after the freeze (user request)"). Context only: it does not
change which flags are strong.

    uv run python code/flag_context.py --scan results/scan/main results/scan/fungal --analysis results/scan/analysis

Writes <analysis>/flag_context.tsv/.md, many_shift_tally.tsv/.md, completeness.md, lineage_sensitivity.tsv/.md,
paired_break_tally.tsv/.md and predictions.md. Uses NCBI E-utilities (public
records only) for SRA platform and other GenBank mitochondrial records of a flagged taxon.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as A  # noqa: E402
import flags as FL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NAMED = {"NC_010769.1": "Aphrocallistes vastus", "NC_068742.1": "Perkinsus marinus", "NC_068743.1": "Perkinsus olseni",
         "NC_036748.1": "Huiozyma sinensis", "NC_004918.1": "Monosporozyma servazzii"}
SPONGES = {"NC_027419.1": "Oopsacas minuta", "NC_028054.1": "Vazella pourtalesii", "NC_072601.1": "Caulophacus iocasicus"}
CTG_GENERA = {"Kazachstania", "Nakaseomyces", "Huiozyma", "Monosporozyma"}
# known insertion/deletion editing lineages (PREREG log, many-shift tally entry), with citations
CALCARONEA = {"Sycettidae", "Leucosoleniidae", "Grantiidae", "Heteropiidae", "Amphoriscidae", "Jenkinidae", "Lelapiidae"}
YANG2017 = {"NC_034794.1", "NC_034798.1"}


def known_editor(m: dict) -> str:
    if m["order_pool_name"] == "Physariida":
        return "myxomycete insertion editing (Mahendran et al. 1991; Gott et al. 2005; Bundschuh et al. 2011)"
    if m["accession"] in YANG2017:
        return "heteroloboseid RNA editing (Yang et al. 2017, Genome Biol Evol 9:1161)"
    if m["class_name"] == "Kinetoplastea":
        return "kinetoplastid U insertion/deletion editing (Benne et al. 1986)"
    if m["class_name"] == "Diplonemea":
        return "diplonemid U insertion editing (Moreira et al. 2016)"
    if m["family_pool_name"] in CALCARONEA:
        return "calcaronean poly(U) insertion editing (Lavrov et al. 2016)"
    return ""


def get(url: str, tries: int = 5) -> str:
    for i in range(tries):
        try:
            time.sleep(0.4)
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception:  # noqa: BLE001
            time.sleep(3 * (i + 1))
    return ""


def gb_comments(accs: set[str]) -> dict[str, str]:
    """COMMENT blocks of wanted records from the local RefSeq gbff and fungal TPA flat files."""
    out = {}
    for path in (ROOT / "data/raw/refseq_mito/mitochondrion.1.genomic.gbff.gz",
                 ROOT / "data/raw/fungal_tpa/BK072095-BK074789.gb.gz"):
        cur, buf, in_c = None, [], False
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("VERSION"):
                    cur = line.split()[1]
                elif line.startswith("COMMENT"):
                    in_c, buf = cur in accs, [line[12:].rstrip()]
                elif in_c and line.startswith("            "):
                    buf.append(line[12:].rstrip())
                elif in_c:
                    out[cur] = " ".join(buf)
                    in_c = False
    return out


def platform(acc: str, comment: str, srr_of: dict) -> dict:
    tech = re.search(r"Sequencing [Tt]echnology\s*:+\s*([^:#]+?)(?:\s{2,}|##|$|\s[A-Z][a-z]+ [A-Za-z]+\s*::)", comment)
    meth = re.search(r"Assembly Method\s*:+\s*([^:#]+?)(?:\s{2,}|##|$|\s[A-Z][a-z]+ [A-Za-z]+\s*::)", comment)
    src = re.search(r"(?:identical to|derived from) ([A-Z]{1,2}_?\d+)", comment)
    srr = re.search(r"entry ((?:SRR|ERR|DRR)\d+)", comment) or (re.match(r"((?:SRR|ERR|DRR)\d+)", srr_of.get(acc, "")))
    r = {"platform": tech.group(1).strip() if tech else "", "assembly_method": meth.group(1).strip() if meth else "",
         "source_accession": src.group(1) if src else "", "sra_run": srr.group(1) if srr else ""}
    if acc.startswith("BK"):
        r["assembly_method"] = r["assembly_method"] or "Sarhan et al. 2025 pipeline (TPA assembly from the SRA run)"
        if r["sra_run"]:
            ids = re.findall(r"<Id>(\d+)</Id>", get(EU + "esearch.fcgi?db=sra&term=" + r["sra_run"]))
            if ids:
                info = get(EU + f"efetch.fcgi?db=sra&id={ids[0]}&rettype=runinfo&retmode=text")
                rows = list(csv.DictReader(info.splitlines()))
                if rows:
                    r["platform"] = f"{rows[0].get('Platform', '')} {rows[0].get('Model', '')} ({rows[0].get('LibraryLayout', '')})".strip()
    return r


def other_assemblies(acc: str, taxid: str, length: int, source: str, scan_same_taxon: list[str],
                     exclude: set[str] = frozenset()) -> str:
    term = f"txid{taxid}[Organism:exp] AND mitochondrion[filter]"
    ids = re.findall(r"<Id>(\d+)</Id>", get(EU + "esearch.fcgi?db=nuccore&retmax=200&term=" + urllib.parse.quote(term)))
    found = []
    if ids:
        js = get(EU + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(ids[:200]))
        try:
            res = json.loads(js)["result"]
        except Exception:  # noqa: BLE001
            res = {}
        for uid in res.get("uids", []):
            d = res[uid]
            av = d.get("accessionversion", "")
            if av.split(".")[0] in (acc.split(".")[0], source) or av.split(".")[0] in exclude or int(d.get("slen", 0)) < 0.8 * length:
                continue
            found.append(f"{av} ({d.get('slen')} bp)")
    parts = []
    if scan_same_taxon:
        parts.append("in scan, same species: " + ", ".join(scan_same_taxon))
    parts.append(f"other GenBank mitochondrial records of the species (>= length threshold; scan records and their sources "
                 f"excluded): {len(found)}"
                 + (" [" + "; ".join(found[:6]) + (" ..." if len(found) > 6 else "") + "]" if found else ""))
    return "; ".join(parts)


# Papers per known lineage and gene, from PRIOR_ART.md only (reading depth in brackets: A abstract, K keyword
# windows, F-agent full text read by an agent, C cited second-hand). No new searching (PREREG log, site-source column).
def site_paper(m: dict, gene: str, codon: str) -> str:
    fam, order, cls = m["family_pool_name"], m["order_pool_name"], m["class_name"]
    if cls in ("Aves",) or order == "Testudines":
        if gene == "nad3":
            return "nad3-174: Mindell 1998 (C); Russell & Beckenbach 2008 (F-agent); Andreu-Sanchez 2021 (F-agent)"
        if order == "Testudines" and gene in ("nad4l", "nad4"):
            return f"turtle CTN-AGN {gene} shift: Russell & Beckenbach 2008 (F-agent; exact position not checked against the paper)"
        if fam == "Podicipedidae" and gene == "nad1":
            return "no paper in PRIOR_ART (grebe nad1 is listed there as annotation only)"
        return "not in PRIOR_ART (not searched)"
    if cls == "Hexactinellida":
        return "glass-sponge TGG: Rosengarten 2008 (F-agent); Haen 2014 (A)"
    if fam == "Perkinsidae":
        return "Perkinsus AGG: Masuda 2010 (A); Zhang 2011 (K); Gornik 2022 (F-agent)"
    if fam == "Saccharomycetaceae" and codon == "CTG":
        return "CTGY yeasts: Langkjaer 2003 (A); Szabooova 2018 (A); Wolters 2023 (F-agent)".replace("Szabooova", "Szaboova")
    if fam in ("Anomiidae",):
        return "Heteranomia TAG: Lin et al. 2026 bioRxiv (K)"
    if fam == "Formicidae" and gene == "cob":
        return "Polyrhachis cob: Beckenbach 2005 (A)"
    if fam == "Ostreidae" and gene == "cob":
        return "Crassostrea cob: Milbury & Gaffney 2005 (A)"
    if fam == "Veneridae":
        return "venerid nad4/cob: Wang 2020 (A)"
    if fam == "Phaeodactylaceae" and gene == "cox1":
        return "Phaeodactylum cox1: Oudot-Le Secq & Green 2011 (A)"
    if fam in ("Rhizophydiaceae", "Spizellomycetaceae") and gene == "cob":
        return "no paper found (PRIOR_ART: annotation only)"
    if fam == "Dipodascaceae":
        return "Magnusiomyces bypassing: Lang 2014 (A) (bypass, not a 1-nt shift)"
    return "not in PRIOR_ART (not searched)"


def record_annotation(acc: str, gene: str, pos: int, junctions: list[dict]) -> str:
    g = gene.upper().replace("NAD", "ND").replace("COB", "CYTB")
    cand = []
    for j in junctions:
        jg = j["gene"].upper().replace("NAD", "ND").replace("COB", "CYTB")
        if j["accession"] == acc and (jg == g or g in jg or jg in g) and abs(int(j["genome_pos"]) - 1 - pos) <= 150:
            cand.append((abs(int(j["genome_pos"]) - 1 - pos), j))
    if not cand:
        return "not annotated in this record"
    off, j = min(cand, key=lambda x: x[0])
    note = j.get("note", "")
    return (f"annotated in this record ({j['gene']} at {j['genome_pos']}, {off} nt from the site"
            f"{'; note: ' + note[:50] if note else ''})")


# ================= additions logged in PREREG before the analysis ran (completeness, units, real-data
# sensitivity, own chance, independent assemblies) =================
import random as _random  # noqa: E402
import units as U  # noqa: E402
import scan as S  # noqa: E402
import fscore as F  # noqa: E402

LINEAGES = [
    ("CTG yeasts", "+1", "CTG", lambda m: m["organism"].split()[0] in CTG_GENERA),
    ("glass sponges", "+1", "TGG", lambda m: m["class_name"] == "Hexactinellida"),
    ("Perkinsus", "+1", "AGG", lambda m: m["family_pool_name"] == "Perkinsidae"),
    ("birds", "+1", "AGT", lambda m: m["class_name"] == "Aves"),
    ("turtles", "+1", "AGT", lambda m: m["order_pool_name"] == "Testudines"),
]
PAPER_GENE_SITES = {  # PRIOR_ART section 1, gene-level only
    "NC_036748.1": ["cox3", "atp6", "cob"],            # K. (Huiozyma) sinensis, Szaboova 2018
    "NC_004918.1": ["atp6", "cox2", "cox3", "cob"],    # K. (Monosporozyma) servazzii, 5 sites
    "NC_010769.1": ["cox3", "nad6"],                   # A. vastus, Rosengarten 2008
}


def completeness(meta, gs) -> list[str]:
    exp = Counter((r["set"], r["code"]) for r in meta.values()
                  if r["set"] in ("refseq_main", "fungal_tpa") and int(r["code"] or 0) > 0 and int(r["length"]) > 0)
    sc = Counter((meta[a]["set"], meta[a]["code"]) for a in gs)
    nocore = Counter((meta[a]["set"], meta[a]["code"]) for a, g in gs.items() if not (set(g["genes_hit"].split(",")) & FL.CORE))
    lines = ["# Completeness (per batch)", "", "| set | code | expected | scored | not scored | scored, no core gene aligned |",
             "|---|---|---|---|---|---|"]
    for k in sorted(exp, key=lambda k: (k[0], int(k[1]))):
        lines.append(f"| {k[0]} | {k[1]} | {exp[k]} | {sc[k]} | {exp[k] - sc[k]} | {nocore[k]} |")
    lines.append(f"| all | | {sum(exp.values())} | {sum(sc.values())} | {sum(exp.values()) - sum(sc.values())} | {sum(nocore.values())} |")
    lines += ["", "A genome without a result row is \"not scored\", never 0 shifts."]
    return lines


def own_chance(n: int, d: int, code: int, gc: float, null, plants, draws: int = 20000, seed: int = 1) -> float:
    """Share of genomes with n random clean shifts (direction d, drawn from the planted null's clean shifts of
    the same stratum) that the identical test flags at any codon."""
    key = null.key(code, gc, d)
    parts = key.split("|")
    cut = null.gc_cut.get(str(code))
    pool = []
    for p in plants:
        if p["clean_h"] != "1" or int(p["shift"]) != d:
            continue
        if parts[0] != "ALL" and p["code"] != parts[0]:
            continue
        if parts[1] in ("lo", "hi") and cut is not None and (("lo" if float(p["gc"]) < cut else "hi") != parts[1]):
            continue
        v = FL.votes({"window": p["window_h"], "codon": p["codon_R_h"], "codon_left": p["codon_L_h"]}, null.cap)
        pool.append((p["gene"], v))
    if not pool or n < 2:
        return 0.0
    stops = {c for c, x in F.codon_table(code).items() if x == "*"}
    q = {c: null.q(code, gc, d, c)[0] for c in FL.CODONS}
    rng = _random.Random(seed)
    cache = {}
    hit = 0
    for _ in range(draws):
        sample = [pool[rng.randrange(len(pool))] for _ in range(n)]
        by = defaultdict(list)
        for gene, v in sample:
            for c in v:
                by[c].append(gene)
        for c, genes in by.items():
            k = len(genes)
            if k < 2 or len(set(genes)) < 2:
                continue
            key2 = (k, c)
            if key2 not in cache:
                cache[key2] = FL.binom_sf(k, n, q[c])
            if cache[key2] <= (FL.ALPHA_STOP if c in stops else FL.ALPHA_SENSE):
                hit += 1
                break
    return hit / draws


def lineage_sensitivity(meta, gs, by_acc, aligned, gflags, null) -> tuple[list[dict], list[str], dict]:
    import controls as C  # noqa: E402  (match, covered)
    J = [j for j in csv.DictReader(open(ROOT / "results/audit/known_junctions_refseq238.tsv"), delimiter="\t")
         if j["group"] != "land_plants"]
    jby = defaultdict(list)
    for j in J:
        jby[j["accession"]].append(j)
    CDS = C.cds_parts()
    flagged = defaultdict(set)
    for f in gflags:
        flagged[f["accession"]].add(f"{int(f['shift']):+d}{f['codon']}")
    rows, summ = [], {}
    seq_cache = {}
    line_stats = defaultdict(lambda: [0, 0])  # line -> [documented core sites, counted]
    oos_sites = []

    def seq_of(a):
        code = meta[a]["code"]
        if code not in seq_cache:
            want = {b for b in gs if meta[b]["code"] == code and any(fn(meta[b]) for _, _, _, fn in LINEAGES)}
            seq_cache[code] = S.read_fasta(ROOT / f"data/work/queries/code_{code}.fna", want) if want else {}
            if meta[a]["set"] == "fungal_tpa":
                seq_cache[code].update(S.read_fasta(ROOT / f"data/work/queries/code_{code}.fna", {a}))
        return seq_cache[code].get(a, "")

    groups_all = list(LINEAGES) + [(n_, "+1", "", f_) for n_, f_ in EXTRA_GROUPS]
    for name, dstr, lcodon, fn in groups_all:
        members = [a for a in gs if fn({**meta[a], "accession": a})]
        stat = Counter()
        for a in members:
            evs = eligible_events(by_acc.get(a, []), null)
            al = aligned.get(a, [])
            sites = []
            ann_genes = set()
            for j in jby.get(a, []):
                d = int(j["shift"])
                fold = 1 if d % 3 == 1 else -1
                pos = int(j["genome_pos"]) - 1
                cds = next(((g_, st_, pa) for g_, st_, pa in CDS[a] if min(x for x, _ in pa) <= pos < max(y for _, y in pa)), None)
                gene = norm_site_gene(cds[0]) if cds else "?"
                ann_genes.add(gene)
                m = C.match(evs, a, pos, fold)
                if gene not in FL.CORE:
                    sites.append({"source": "record annotation", "gene": gene, "pos": pos + 1, "shift": f"{fold:+d}",
                                  "outcome": "non-core gene (not tested)"})
                    continue
                if m is None:
                    if any(int(e["shift"]) == fold and e["gene"] == gene for e in evs):
                        out = "misplaced (shift elsewhere in the gene)"
                    else:
                        out = "aligned_no_event" if C.covered(al, pos) else "not_aligned"
                elif m["class"] != "clean":
                    out = m["class"]
                else:
                    v = FL.votes(m, null.cap)
                    own = set()
                    s_ = seq_of(a)
                    if cds and s_:
                        so = s_ if cds[1] == 1 else s_.translate(S.COMP)[::-1]
                        opos = (pos if cds[1] == 1 else len(s_) - 1 - pos) if fold == 1 else (pos + 1 if cds[1] == 1 else len(s_) - pos)
                        tr = F.trigger(F.Event(a, "", cds[1], opos, fold, fold, 0, 0), so, F.codon_table(int(meta[a]["code"])),
                                       frame0=F.old_frame0(cds[2], cds[1], len(s_), opos))
                        own = {tr["codon"], tr["codon_left"]}
                    out = "counted" if (v & own) else ("clean_wide_window" if not v else "clean_other_codon")
                cpos = int(j["cds_pos"]) if str(j.get("cds_pos", "")).isdigit() else -1
                line = ("nad3-174 (birds and turtles)" if name in ("birds", "turtles") and gene == "nad3" and 160 <= cpos <= 190
                        else f"{name} (excl. nad3-174)" if name in ("birds", "turtles") else name)
                sites.append({"source": "record annotation", "gene": gene, "pos": pos + 1, "shift": f"{fold:+d}", "outcome": out,
                              "line": line})
            for st in OOS_SITES.get(a, []):
                gene, fold = st["gene"], st["shift"]
                cdsx = next(((g_, sd_, pa) for g_, sd_, pa in CDS[a] if norm_site_gene(g_) == gene), None)
                if st["cds_pos"] and cdsx:
                    pos = cds_genome_pos(cdsx[2], cdsx[1], st["cds_pos"])
                    m = C.match(evs, a, pos, fold)
                    if m is None:
                        out = ("misplaced (shift elsewhere in the gene)" if any(int(e["shift"]) == fold and e["gene"] == gene for e in evs)
                               else "aligned_no_event" if C.covered(al, pos) else "not_aligned")
                    elif m["class"] != "clean":
                        out = m["class"]
                    else:
                        out = "detected clean at site (codon not checked)"
                    sites.append({"source": st["source"], "gene": gene, "pos": pos + 1, "shift": f"{fold:+d}", "outcome": out,
                                  "line": name, "oos": 1})
                else:
                    ge = [e for e in evs if e["gene"] == gene and int(e["shift"]) == fold]
                    out = ("detected clean in gene (codon not checked)" if any(e["class"] == "clean" for e in ge)
                           else (min(ge, key=lambda e: e["class"])["class"] if ge else "no shift found in the gene"))
                    sites.append({"source": st["source"], "gene": gene, "pos": "", "shift": f"{fold:+d}", "outcome": out,
                                  "line": name, "oos": 1})
            for gene in PAPER_GENE_SITES.get(a, []):
                if gene in ann_genes:
                    continue
                ge = [e for e in evs if e["gene"] == gene and int(e["shift"]) == 1]
                if any(e["class"] == "clean" and lcodon in FL.votes(e, null.cap) for e in ge):
                    out = "counted"
                elif ge:
                    best = min(ge, key=lambda e: e["class"] != "clean")
                    out = best["class"] if best["class"] != "clean" else "clean_other_codon"
                else:
                    out = "no shift found in the gene"
                sites.append({"source": "paper (gene only)", "gene": gene, "pos": "", "shift": "+1", "outcome": out, "line": name})
            oos_sites += [{**s_, "acc": a, "org": meta[a]["organism"]} for s_ in sites if s_.get("oos")]
            sites = [s_ for s_ in sites if not s_.get("oos")] + [s_ for s_ in sites if s_.get("oos")]
            n_cnt = sum(s_["outcome"] == "counted" for s_ in sites)
            sites_t = [s_ for s_ in sites if s_["outcome"] != "non-core gene (not tested)"]
            fl = sorted(flagged.get(a, []))
            rows.append({"lineage": name, "accession": a, "organism": meta[a]["organism"], "documented_sites": len(sites_t),
                         "counted": n_cnt, "flagged": ",".join(fl) or "no",
                         "sites": "; ".join(f"{s_['gene']}{('@' + str(s_['pos'])) if s_['pos'] else ''} [{s_['source']}]: {s_['outcome']}" for s_ in sites)
                         or "no documented sites"})
            stat["genomes"] += 1
            stat["with_sites"] += bool(sites_t)
            stat["with_2plus_sites"] += len(sites_t) >= 2
            stat["sites"] += len(sites_t)
            stat["counted"] += n_cnt
            stat["genomes_2plus_counted"] += n_cnt >= 2
            stat["flagged"] += bool(fl)
            stat["flagged_with_2plus_sites"] += bool(fl) and len(sites_t) >= 2
            for s_ in sites_t:
                stat["out:" + s_["outcome"]] += 1
                if s_.get("oos"):
                    continue
                line_stats[s_["line"]][0] += 1
                line_stats[s_["line"]][1] += s_["outcome"] == "counted"
        summ[name] = stat
    md = ["# Real-data sensitivity in the known lineages (added after the freeze; descriptive)", "",
          "Planted sensitivity for comparison (PREREG 1.12): a planted shift is counted at its codon 0.553 of the time "
          "(random sites) and 0.518 (same-codon trial sites). A genome with 2, 3 or 5 planted same-codon sites flags "
          "0.317, 0.483 or 0.717 of the time.", "",
          "| lineage | genomes scored | with documented sites | with >= 2 | documented sites | counted at its codon | "
          "genomes with >= 2 counted | flagged (any) | flagged, of those with >= 2 sites | outcomes of documented sites |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for name, *_ in groups_all:
        s_ = summ[name]
        outs = ", ".join(f"{k[4:]} {v}" for k, v in sorted(s_.items(), key=lambda x: -x[1]) if k.startswith("out:"))
        md.append(f"| {name} | {s_['genomes']} | {s_['with_sites']} | {s_['with_2plus_sites']} | {s_['sites']} | "
                  f"{s_['counted']} ({s_['counted'] / max(1, s_['sites']):.2f}) | {s_['genomes_2plus_counted']} | {s_['flagged']} | "
                  f"{s_['flagged_with_2plus_sites']} of {s_['with_2plus_sites']} | {outs} |")
    lines_ = sorted(line_stats.items(), key=lambda kv: (kv[0].startswith("nad3-174"), kv[0]))
    md += ["", "### Real per-site rate by lineage (documented core-gene sites counted at their codon)", "",
           "| lineage line | documented sites | counted | rate |", "|---|---|---|---|"]
    md += [f"| {ln} | {n_} | {c_} | {c_ / n_:.2f} |" for ln, (n_, c_) in lines_ if n_]
    head = [(ln, n_, c_) for ln, (n_, c_) in lines_ if n_ and not ln.startswith("nad3-174")]
    rates = [c_ / n_ for _, n_, c_ in head]
    if rates:
        r = sum(rates) / len(rates)
        pooled = sum(c_ for *_, c_ in head) / sum(n_ for _, n_, _ in head)
        if len(rates) >= 5:
            rng = _random.Random(1)
            boots = sorted(sum(rng.choice(rates) for _ in rates) / len(rates) for _ in range(2000))
            lo, hi = boots[49], boots[1949]
            iv = f"bootstrap over {len(rates)} lineages, 95% {lo:.2f}-{hi:.2f}"
        else:
            lo, hi = min(rates), max(rates)
            iv = f"range across {len(rates)} lineages {lo:.2f}-{hi:.2f} (too few lineages for a bootstrap)"
        md += ["", f"**Headline real per-site rate, nad3-174 excluded: {r:.2f}** (mean of per-lineage rates; {iv}; pooled "
               f"{pooled:.2f} over {sum(n_ for _, n_, _ in head)} sites), against 0.553 for planted shifts.",
               "Implied chance that a genome with m documented same-codon sites gets >= 2 counted (and then flags, since B5 "
               "conditional = 1.000): " + "; ".join(f"m = {m}: {p_two_plus(r, m):.2f} ({p_two_plus(lo, m):.2f}-{p_two_plus(hi, m):.2f})"
                                                       for m in (2, 3, 5)) + "."]
    oos = oos_sites
    n_oos = len(oos)
    det = sum(s_["outcome"].startswith("detected clean") or s_["outcome"] == "counted" for s_ in oos)
    md += ["", f"**Out-of-sample sites** (documented, not annotated in the record so not in B2, not in a named positive): "
           f"{n_oos} sites; detected as a clean shift (at the site, or in the gene for gene-level sites) {det} of {n_oos}"
           + (f" ({det / n_oos:.2f})" if n_oos else "") + ". Too few for a stable rate; listed in full:"]
    md += [f"- {s_['acc']} {s_['org']}: {s_['gene']}{('@' + str(s_['pos'])) if s_['pos'] else ''} {s_['shift']} [{s_['source']}]: {s_['outcome']}"
           for s_ in oos]
    md += ["Write-up limits (PREREG log): a negative rules out lineages only at these real-data rates; the pair filter removes real "
           "sites where +1 and +2 shifts sit close together; the chance figures assume errors at random codons, and assembly "
           "errors cluster, so a noisy genome can be strong without being real."]
    md += ["", "Genomes with no documented sites are not misses. Per-genome detail for the small lineages:", ""]
    for r in rows:
        if r["lineage"] in ("CTG yeasts", "glass sponges", "Perkinsus"):
            md.append(f"- {r['lineage']}: {r['organism']} ({r['accession']}): flagged {r['flagged']}; {r['sites']}")
    md.append("")
    md.append("Birds and turtles: full per-genome list in lineage_sensitivity.tsv.")
    return rows, md, summ


def norm_site_gene(g: str) -> str:
    x = g.lower().removeprefix("other:")
    x = re.sub(r"p$", "", x)
    for a_, b_ in (("cytb", "cob"), ("coxi", "cox1"), ("coxii", "cox2"), ("coxiii", "cox3"), ("nd", "nad")):
        if x.startswith(a_) and not x.startswith("nad"):
            x = b_ + x[len(a_):]
    return x


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * ((ph * (1 - ph) / n + z * z / (4 * n * n)) ** 0.5) / den
    return (max(0.0, c - h), min(1.0, c + h))


def p_two_plus(r: float, m: int) -> float:
    """Chance that >= 2 of m documented sites are counted, per-site rate r (independent)."""
    return 1 - (1 - r) ** m - m * r * (1 - r) ** (m - 1)


def mean_identity(accs: set[str], scan_dirs, meta) -> dict[str, float]:
    """Identical / aligned amino-acid columns over each genome's kept hits (same cull as the scan)."""
    if not accs:
        return {}
    _, ref_family, _ = S.load_meta()
    hits = defaultdict(list)
    for d in scan_dirs:
        for hf in Path(d).glob("hits_code*.tsv"):
            with open(hf) as fh:
                for line in fh:
                    q = line.split("\t", 1)[0]
                    if q in accs:
                        f = line.rstrip("\n").split("\t")
                        hits[q].append(F.Hit(f[0], f[1], int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]),
                                             int(f[8]), float(f[9]), f[11]))
    out = {}
    for a in accs:
        kept = F.cull(hits.get(a, []), meta[a]["family_pool"], ref_family)
        ident = tot = 0
        for h in kept:
            _, _, _, cols = F.walk(h, "")
            for c in cols:
                if c[0] == "aa":
                    tot += 1
                    ident += c[1]
        out[a] = ident / tot if tot else float("nan")
    return out


def paired_tally(meta, gs, by_acc, null, scan_dirs=()) -> tuple[list[dict], list[str]]:
    """Exploratory safety net for the pair-filter blind spot (PREREG log; post-freeze, never a flag)."""
    out = []
    for a in gs:
        ps = [e for e in eligible_events(by_acc.get(a, []), null) if e["class"] == "paired" and e["gene"] in FL.CORE]
        if len(ps) < 2:
            continue
        grp = defaultdict(list)
        for e in ps:
            grp[(e["strand"], e["subject"])].append(e)
        pairs = []
        for evs in grp.values():
            evs.sort(key=lambda e: int(e["pos_oriented"]))
            used = set()
            for i, e in enumerate(evs):
                if i in used:
                    continue
                cand = [(abs(int(o["pos_oriented"]) - int(e["pos_oriented"])), j) for j, o in enumerate(evs)
                        if j != i and j not in used and int(o["shift"]) == -int(e["shift"])
                        and abs(int(o["pos_oriented"]) - int(e["pos_oriented"])) <= 150]
                if not cand:
                    continue
                _, j = min(cand)
                used |= {i, j}
                first = e if int(e["pos_oriented"]) <= int(evs[j]["pos_oriented"]) else evs[j]
                pairs.append(first)
        by = defaultdict(list)
        for f in pairs:
            for c in FL.votes(f, null.cap):
                by[(int(f["shift"]), c)].append(f)
        for (d, c), fs in by.items():
            if len(fs) >= 2:
                out.append({"accession": a, "organism": meta[a]["organism"], "family": meta[a]["family_pool_name"],
                            "phylum": meta[a]["phylum_name"] or meta[a]["kingdom_name"], "known_lineage": int(A.is_known(meta[a])),
                            "first_break": f"{d:+d} {c}", "pairs": len(fs), "genes": len({f['gene'] for f in fs}),
                            "total_pairs_in_genome": len(pairs), "total_paired_breaks": len(ps),
                            "sites": "; ".join(f"{f['gene']}@{f['pos_fwd']}" for f in fs)})
    out.sort(key=lambda r: (-r["pairs"], -r["genes"]))
    mid = mean_identity({r["accession"] for r in out}, scan_dirs, meta)
    for r in out:
        r["mean_alignment_identity"] = round(mid.get(r["accession"], float("nan")), 3)
    md = ["# Paired breaks recurring at one first-break codon (exploratory safety net; added after the freeze; never a flag)", "",
          f"{len(out)} genome-codon rows with >= 2 pairs sharing a first-break codon (window rule).", "",
          "Read each recurrence against the genome's total paired breaks and its mean alignment identity: divergent genomes "
          "produce detours by chance.", "",
          "| genome | group | first break | pairs at this codon (genes) | pairs in genome | paired breaks in genome | "
          "mean alignment identity | known lineage | sites |", "|---|---|---|---|---|---|---|---|---|"]
    md += [f"| {r['organism']} ({r['accession']}) | {r['family']}, {r['phylum']} | {r['first_break']} | {r['pairs']} ({r['genes']}) | "
           f"{r['total_pairs_in_genome']} | {r['total_paired_breaks']} | {r['mean_alignment_identity']:.2f} | "
           f"{'yes' if r['known_lineage'] else 'no'} | {r['sites']} |" for r in out]
    return out, md


def eligible_events(evs, null):
    return evs if null.split_votes else [e for e in evs if e["via"] != "split"]


def species_taxid(taxid: str) -> str:
    nodes = getattr(species_taxid, "_nodes", None)
    if nodes is None:
        nodes = {}
        with open(ROOT / "data/raw/taxdump/nodes.dmp") as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                nodes[f[0]] = (f[1], f[2])
        species_taxid._nodes = nodes
    t = taxid
    for _ in range(10):
        if t not in nodes:
            return taxid
        par, rank = nodes[t]
        if rank == "species":
            return t
        t = par
    return taxid


# ---- known groups from results/audit/controls_known_cases.tsv beyond the five lineages (PREREG log, out-of-sample)
EXTRA_GROUPS = [
    ("oysters", lambda m: m["family_pool_name"] == "Ostreidae"),
    ("venerid clams", lambda m: m["family_pool_name"] == "Veneridae"),
    ("Polyrhachis ants", lambda m: m["organism"].startswith("Polyrhachis")),
    ("Arctica", lambda m: m["family_pool_name"] == "Arcticidae"),
    ("Phaeodactylum", lambda m: m["organism"].startswith("Phaeodactylum")),
    ("Dipterophagus", lambda m: m["organism"].startswith("Dipterophagus")),
    ("chytrids (Rhizophydium, Spizellomyces)", lambda m: m["family_pool_name"] in ("Rhizophydiaceae", "Spizellomycetaceae")),
    ("other RefSeq single annotations", lambda m: m["accession"] in ("NC_020309.1", "NC_026914.1")),
]
# out-of-sample documented sites: in controls_known_cases.tsv, not annotated in the record (not in B2), not a named positive
OOS_SITES = {
    "NC_058598.1": [{"gene": "nad5", "cds_pos": 291, "shift": -1, "source": "Towett-Kirui 2022 (K31)"}],
    "NC_001947.1": [{"gene": "nad3", "cds_pos": 135, "shift": 1, "source": "Russell & Beckenbach 2008 (K23; nad3-135)"},
                    {"gene": "nad4l", "cds_pos": 99, "shift": 1, "source": "Russell & Beckenbach 2008 (K23; 'near nt 99')"},
                    {"gene": "nad4l", "cds_pos": 262, "shift": 1, "source": "Russell & Beckenbach 2008 (K23; 'near nt 262')"}],
    "NC_004691.1": [{"gene": "cox2", "cds_pos": None, "shift": 1, "source": "record note, 3' end (K15)"}],
    "NC_016739.1": [{"gene": "cox1", "cds_pos": None, "shift": 1, "source": "Oudot-Le Secq & Green 2011 (K29)"}],
    "NC_003052.1": [{"gene": "cob", "cds_pos": None, "shift": 1, "source": "record note (K33)"}],
}
SPONGE_DOCS = {
    "NC_010769.1": "documented: cox3 aa 140 and nad6 aa 58 (Rosengarten 2008; record annotations)",
    "NC_027419.1": "documented absence: its source paper (Jourda et al. 2015, Genome Announc; full text read) states "
                   "'No occurrence of the +1 frameshifting previously reported in other glass sponge mitochondria was observed.'",
    "NC_028054.1": "unknown: source paper is Haen et al. 2014 (+1 insertions in its new glass-sponge genomes), but only the abstract "
                   "is readable (not open access) and it names no species or positions; the record annotates no shift",
    "NC_072601.1": "none documented: no paper cited by the record; no shift annotated",
}


def cds_genome_pos(parts, strand, cds_pos):
    """0-based genome position of 1-based CDS position cds_pos for a CDS given as forward parts."""
    off = cds_pos - 1
    seq_parts = sorted(parts) if strand == 1 else sorted(parts, reverse=True)
    for x, y in seq_parts:
        n = y - x
        if off < n:
            return x + off if strand == 1 else y - 1 - off
        off -= n
    return None


def removed_breaks(meta, gs, by_acc, null) -> list[str]:
    """Every break DIAMOND wrote in the glass sponges and Polyrhachis, with the filter that removed it (descriptive)."""
    accs = [a for a in gs if meta[a]["class_name"] == "Hexactinellida" or meta[a]["organism"].startswith("Polyrhachis")]
    md = ["# Breaks written in the glass sponges and Polyrhachis (descriptive; added after the freeze)", "",
          "Codons are the old-frame codons holding the leftmost/rightmost equivalent placement; 'next' is the codon after the "
          "rightmost one. **TGGA** marks the A. vastus signature: a TGG codon whose +1 reading (skipping its T) is GGA.", ""]
    verdicts = {}
    for a in accs:
        code = int(meta[a]["code"])
        s_ = S.read_fasta(ROOT / f"data/work/queries/code_{code}.fna", {a}).get(a, "")
        table = F.codon_table(code)
        evs = eligible_events(by_acc.get(a, []), null)
        md.append(f"## {meta[a]['organism']} ({a}; code {code}): {len(evs)} breaks written")
        if meta[a]["class_name"] == "Hexactinellida":
            md.append(f"- documented sites: {SPONGE_DOCS.get(a, 'not in PRIOR_ART')}")
        else:
            md.append("- documented sites: unknown (Beckenbach 2005 reports cytb sites in 12 of 30 Polyrhachis species; only its "
                      "abstract was read, so whether this species is one of them is unknown; the record annotates no shift)")
        sig_removed = sig_clean = 0
        md += ["", "| gene | position | shift | codons L/R | next codon | class (filter) | TGGA signature |", "|---|---|---|---|---|---|---|"]
        for e in sorted(evs, key=lambda e: int(e["pos_fwd"])):
            so = s_ if int(e["strand"]) == 1 else s_.translate(S.COMP)[::-1]
            ev = F.Event(a, e["subject"], int(e["strand"]), int(e["pos_oriented"]), int(e["shift"]), int(e["raw"]), 0, 0)
            tr = F.trigger(ev, so, table)
            f0 = ev.pos % 3
            cs = tr["tpos"] - ((tr["tpos"] - f0) % 3)
            nxt = so[cs + 3:cs + 6]
            # A. vastus signature: an in-frame TGG followed by A (its +1 reading is GGA) within the placement window
            sig = "TGG" in (e["codon"], e["codon_left"]) and any(
                so[t:t + 4] == "TGGA" for t in range(max(0, cs - 9), cs + 10) if (t - f0) % 3 == 0)
            if sig:
                sig_clean += e["class"] == "clean"
                sig_removed += e["class"] != "clean"
            md.append(f"| {e['gene']} | {e['pos_fwd']} | {int(e['shift']):+d} | {e['codon_left']}/{e['codon']} | {nxt} | {e['class']} | "
                      f"{'**TGGA**' if sig else ''} |")
        md.append("")
        if meta[a]["class_name"] == "Hexactinellida" and a != "NC_010769.1":
            doc = SPONGE_DOCS.get(a, "")
            if sig_removed:
                verdicts[a] = f"filters (TGG-signature breaks seen and removed: {sig_removed})"
            elif doc.startswith("documented absence") or doc.startswith("none documented"):
                verdicts[a] = "absence (no sites documented, none seen)" if not sig_clean else "seen, clean"
            else:
                verdicts[a] = "unknown (documentation not readable; no TGG-signature break seen)"
    removed_breaks.verdicts = verdicts
    return md


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", nargs="+", required=True)
    ap.add_argument("--analysis", required=True)
    args = ap.parse_args()
    out = Path(args.analysis)
    meta = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    null = FL.Null(ROOT / "results/controls/plant_null/null.json")
    ev, gs, ins, _ = A.load(args.scan)
    by_acc = defaultdict(list)
    for e in ev:
        by_acc[e["accession"]].append(e)
    clean = {a: FL.clean_core(by_acc.get(a, []), null.split_votes) for a in gs}
    gflags = list(csv.DictReader(open(out / "genome_flags.tsv"), delimiter="\t")) if (out / "genome_flags.tsv").stat().st_size else []
    pflags = [r for r in (csv.DictReader(open(out / "pool_tests.tsv"), delimiter="\t") if (out / "pool_tests.tsv").stat().st_size else [])
              if r["flag"] == "1"]
    srr_of = {}
    for r in csv.DictReader(open(ROOT / "data/raw/fungal_tpa/All_fungal_mitochondrial_genomes_accessions.tsv"), delimiter="\t"):
        srr_of[r["NCBI_accession_number"] + ".1"] = r["seq_id"]
    by_taxon = defaultdict(list)
    for a in gs:
        by_taxon[meta[a]["taxid"]].append(a)
    flagged_keys = {(f["accession"], f["codon"], f["shift"]) for f in gflags}
    unit = U.unit_of()
    plants = list(csv.DictReader(open(ROOT / "results/controls/plant_null/plants.tsv"), delimiter="\t"))
    aligned = {}
    for d_ in args.scan:
        if (Path(d_) / "aligned_fwd.json").exists():
            aligned.update(json.load(open(Path(d_) / "aligned_fwd.json")))
    sp_of = {a: species_taxid(meta[a]["taxid"]) for a in gs}
    by_species = defaultdict(list)
    for a in gs:
        by_species[sp_of[a]].append(a)
    fam_len = defaultdict(list)
    for a in gs:
        fam_len[meta[a]["family_pool"]].append(int(meta[a]["length"]))
    fam_med = {k: sorted(v)[len(v) // 2] for k, v in fam_len.items()}
    comments = gb_comments({f["accession"] for f in gflags})
    junctions = list(csv.DictReader(open(ROOT / "results/audit/known_junctions_refseq238.tsv"), delimiter="\t"))

    rows = []
    for f in gflags:
        a, codon, d = f["accession"], f["codon"], int(f["shift"])
        cl = clean[a]
        sites = [e for e in cl if int(e["shift"]) == d and codon in FL.votes(e, null.cap)]
        rel = []
        for b in gs:
            if b == a or not (meta[b]["family_pool"] == meta[a]["family_pool"] or meta[b]["order_pool"] == meta[a]["order_pool"]):
                continue
            k = sum(1 for e in clean[b] if int(e["shift"]) == d and codon in FL.votes(e, null.cap))
            if k:
                rel.append(f"{b} {meta[b]['organism']} k={k}{' FLAGGED' if (b, codon, str(d)) in flagged_keys else ''}"
                           f"{' (same assembly/strain)' if unit.get(b) == unit.get(a) else ''}")
        pl = platform(a, comments.get(a, ""), srr_of)
        n_dir = sum(int(e["shift"]) == d for e in cl)
        code_a, gc_a = int(meta[a]["code"]), float(gs[a]["gc"])
        stops_a = {c for c, x in F.codon_table(code_a).items() if x == "*"}
        qs_a = {c: null.q(code_a, gc_a, d, c)[0] for c in FL.CODONS}
        al_a = {c: (FL.ALPHA_STOP if c in stops_a else FL.ALPHA_SENSE) for c in FL.CODONS}
        chance_mc = own_chance(n_dir, d, code_a, gc_a, null, plants)
        chance_ub = FL.expected_chance(n_dir, qs_a, al_a)
        same_unit = [b for b in gs if b != a and unit.get(b) == unit.get(a)]
        indep = [b for b in by_species[sp_of[a]] if unit.get(b) != unit.get(a)]
        genus = meta[a]["organism"].split()[0]
        congeners = [b for b in gs if b != a and meta[b]["organism"].split()[0] == genus]

        def shifts_of(b):
            c_ = clean[b]
            return f"{len(c_)} clean (+1 {sum(int(e['shift']) == 1 for e in c_)}, -1 {sum(int(e['shift']) == -1 for e in c_)})"
        thr = int(0.8 * min(int(meta[a]["length"]), fam_med.get(meta[a]["family_pool"], int(meta[a]["length"]))))

        def sp_excl(recs):
            com = gb_comments(set(recs))
            ex = {r_.split(".")[0] for r_ in recs}
            for r_ in recs:
                m_ = re.search(r"(?:identical to|derived from) ([A-Z]{1,2}_?\d+)", com.get(r_, ""))
                if m_:
                    ex.add(m_.group(1))
            return ex
        rows.append({
            "accession": a, "organism": meta[a]["organism"], "family": meta[a]["family_pool_name"], "order": meta[a]["order_pool_name"],
            "flag": f"{d:+d} {codon}", "strong": f["strong"], "known_lineage": f["known_lineage"], "p": f["p"], "E": f["E"],
            "sites_at_codon": "; ".join(f"{e['gene']}@{e['pos_fwd']}" for e in sites),
            "different_genes": int(len({e['gene'] for e in sites}) == len(sites)), "n_genes_at_codon": len({e["gene"] for e in sites}),
            "clean_shifts_total": len(cl), "plus1": sum(int(e["shift"]) == 1 for e in cl), "minus1": sum(int(e["shift"]) == -1 for e in cl),
            **pl,
            "related_same_codon": f"{len(rel)}: " + " | ".join(rel[:6]) if rel else "0",
            "site_sources": " || ".join(f"{e['gene']}@{e['pos_fwd']}: paper: {site_paper(meta[a], e['gene'], codon)}; record: "
                                        f"{record_annotation(a, e['gene'], int(e['pos_fwd']), junctions)}" for e in sites)
            if f["known_lineage"] == "1" else "",
            "n_clean_flag_direction": n_dir, "own_chance_mc": round(chance_mc, 4), "own_chance_upper_bound": round(chance_ub, 4),
            "same_assembly_records": ", ".join(f"{b} ({shifts_of(b)})" for b in same_unit),
            "second_assembly": other_assemblies(a, sp_of[a], thr / 0.8, pl["source_accession"],
                                                [f"{b} {meta[b]['organism']} ({shifts_of(b)})" for b in indep],
                                                exclude=sp_excl(by_species[sp_of[a]])),
            "congeners_in_scan": "; ".join(f"{b} {meta[b]['organism']}{' [same assembly]' if unit.get(b) == unit.get(a) else ''}: "
                                           f"{shifts_of(b)}{' FLAGGED' if any(k_[0] == b for k_ in flagged_keys) else ''}"
                                           for b in congeners[:15]) + (f" (+{len(congeners) - 15} more)" if len(congeners) > 15 else ""),
        })
    A.write(out / "flag_context.tsv", rows)

    # ---- many-shift tally (added after the freeze; descriptive)
    tally = list(csv.DictReader(open(out / "tally_three_gene_genomes.tsv"), delimiter="\t")) \
        if (out / "tally_three_gene_genomes.tsv").exists() and (out / "tally_three_gene_genomes.tsv").stat().st_size else []
    tcom = gb_comments({r["accession"] for r in tally})
    info = {}
    for r in tally:
        a = r["accession"]
        cl = clean[a]
        n_p = sum(int(e["shift"]) == 1 for e in cl)
        n_m = len(cl) - n_p
        info[a] = {"n": len(cl), "maj": "+1" if n_p >= n_m else "-1", "share": max(n_p, n_m) / max(1, len(cl)), "plus1": n_p, "minus1": n_m}
    trows = []
    for r in sorted(tally, key=lambda r: -info[r["accession"]]["n"]):
        a, x = r["accession"], info[r["accession"]]
        rel = [b for b in info if b != a and info[b]["maj"] == x["maj"] and
               (meta[b]["family_pool"] == meta[a]["family_pool"] or meta[b]["order_pool"] == meta[a]["order_pool"])]
        ed = known_editor({**meta[a], "accession": a})
        pl = platform(a, tcom.get(a, ""), srr_of)
        trows.append({"accession": a, "organism": meta[a]["organism"], "family": meta[a]["family_pool_name"],
                      "order": meta[a]["order_pool_name"], "phylum": meta[a]["phylum_name"] or meta[a]["kingdom_name"],
                      "clean_shifts": x["n"], "genes": r["n_genes"], "plus1": x["plus1"], "minus1": x["minus1"],
                      "majority": x["maj"], "majority_share": round(x["share"], 3), **pl,
                      "related_same_pattern": ",".join(f"{b} {meta[b]['organism']}" for b in rel[:6]) + (f" (+{len(rel)-6})" if len(rel) > 6 else ""),
                      "known_editor": ed,
                      "mark": ("possible editing or noisy assembly" if (not ed and x["share"] >= 0.80) else "")})
    A.write(out / "many_shift_tally.tsv", trows)
    tm = ["# Many-shift tally: genomes with clean shifts in >= 3 different core genes (exploratory; added after the freeze)", "",
          f"{len(trows)} genomes, sorted by clean shifts. Marked \"possible editing or noisy assembly\": not a known editor and >= 80% "
          f"of clean shifts in one direction ({sum(bool(t['mark']) for t in trows)} genomes).", "",
          "| genome | group | clean shifts (genes) | +1 / -1 | majority share | platform; assembly | related, same pattern | known editor | mark |",
          "|---|---|---|---|---|---|---|---|---|"]
    for t in trows:
        tm.append(f"| {t['organism']} ({t['accession']}) | {t['family']}, {t['phylum']} | {t['clean_shifts']} ({t['genes']}) | "
                  f"{t['plus1']} / {t['minus1']} | {t['majority_share']:.2f} {t['majority']} | "
                  f"{t['platform'] or 'not given'}; {t['assembly_method'] or 'not given'} | {t['related_same_pattern'] or '-'} | "
                  f"{t['known_editor'] or '-'} | {t['mark']} |")
    (out / "many_shift_tally.md").write_text("\n".join(tm) + "\n")

    # ---- predictions, scored as fixed in the PREREG log before results
    flagged = defaultdict(set)
    for f in gflags:
        flagged[f["accession"]].add(f"{int(f['shift']):+d}{f['codon']}")
    p1 = [a for a in NAMED if flagged.get(a)]
    p2 = [a for a in SPONGES if "+1TGG" in flagged.get(a, set())]
    ctg = [a for a in gs if meta[a]["organism"].split()[0] in CTG_GENERA and a not in NAMED]
    p3 = [a for a in ctg if "+1CTG" in flagged.get(a, set())]
    birdturtle = [f for f in gflags if meta[f["accession"]]["class_name"] == "Aves" or meta[f["accession"]]["order_pool_name"] == "Testudines"]
    bt_nad3 = [f for f in birdturtle if any(e["gene"] == "nad3" and int(e["shift"]) == int(f["shift"]) and f["codon"] in FL.votes(e, null.cap)
                                              for e in clean[f["accession"]])]
    bt_pool = []
    for r in pflags:
        members = [a for a in gs if meta[a][r["rank"]] == r["pool"]]
        if not any(meta[a]["class_name"] == "Aves" or meta[a]["order_pool_name"] == "Testudines" for a in members):
            continue
        sites = A.pool_sites([e for a in members for e in clean[a]], null.cap)
        hit = [s for s in sites if s["shift"] == int(r["shift"]) and r["codon"] in s["codons"]]
        bt_pool.append((r, any(s["gene"] == "nad3" for s in hit)))
    summ = (out / "flags_summary.md").read_text()
    m = re.search(r"own counts: ([\d.]+)", summ)
    e_chance = float(m.group(1)) if m else float("nan")
    new_g = [f for f in gflags if f["known_lineage"] == "0"]
    magn = sum(1 for i in ins if i["accession"] == "NC_024095.1" and i["class"] == "bypass_like")
    dip = Counter(meta[i["accession"]]["family_pool_name"] for i in ins if i["class"] == "bypass_like"
                  and meta[i["accession"]]["family_pool_name"] in ("Dipodascaceae", "Trichomonascaceae"))
    dip_g = {fam: len({i["accession"] for i in ins if i["class"] == "bypass_like" and meta[i["accession"]]["family_pool_name"] == fam})
             for fam in dip}
    P = [
        ("P1 named positives flag", len(p1) == 5, f"{len(p1)}/5 flag ({', '.join(NAMED[a] for a in p1)})"),
        ("P2 another glass sponge flags at +1 TGG", len(p2) >= 1,
         f"{len(p2)}/3 ({', '.join(SPONGES[a] for a in p2) or 'none'}); their flags: " +
         "; ".join(f"{SPONGES[a]} {sorted(flagged.get(a, [])) or 'none'}" for a in SPONGES)),
        ("P3 another CTG-group yeast flags at +1 CTG", len(p3) >= 1,
         f"{len(p3)} of {len(ctg)} scanned ({', '.join(meta[a]['organism'] for a in p3) or 'none'})"),
        ("P4 no bird/turtle flag involves nad3", not bt_nad3 and not any(x for _, x in bt_pool),
         f"bird/turtle genome flags {len(birdturtle)} (with a nad3 site at the codon: {len(bt_nad3)}); bird/turtle pool flags "
         f"{len(bt_pool)} (involving nad3: {sum(x for _, x in bt_pool)})"),
        ("P5 non-known flags <= chance expectation + 5", len(new_g) <= e_chance + 5,
         f"{len(new_g)} non-known single-genome flags; recomputed chance expectation {e_chance}"),
        ("P6 Magnusiomyces in the bypass-like tally", magn >= 1,
         f"{magn} bypass-like inserts in NC_024095; genomes with bypass-like inserts (description): {dip_g or 'none'}"),
    ]
    # completeness and real-data sensitivity (PREREG log, before the analysis ran)
    comp = completeness(meta, gs)
    (out / "completeness.md").write_text("\n".join(comp) + "\n")
    lrows, lmd, lsumm = lineage_sensitivity(meta, gs, by_acc, aligned, gflags, null)
    rbmd = removed_breaks(meta, gs, by_acc, null)
    (out / "removed_breaks_known_groups.md").write_text("\n".join(rbmd) + "\n")
    prows, pmd = paired_tally(meta, gs, by_acc, null, args.scan)
    A.write(out / "paired_break_tally.tsv", prows)
    (out / "paired_break_tally.md").write_text("\n".join(pmd) + "\n")
    A.write(out / "lineage_sensitivity.tsv", lrows)
    (out / "lineage_sensitivity.md").write_text("\n".join(lmd) + "\n")

    def explain(lineage: str, exclude: set) -> str:
        rr = [r for r in lrows if r["lineage"] == lineage and r["accession"] not in exclude]
        if not rr:
            return f"no other {lineage} genomes scored"
        two = [r for r in rr if r["documented_sites"] >= 2]
        cnt2 = [r for r in rr if r["counted"] >= 2]
        if not two:
            return (f"not explained by sensitivity: none of the {len(rr)} other {lineage} genomes scored has >= 2 documented "
                    f"sites, so there was nothing documented to pile up ({sum(r['documented_sites'] for r in rr)} documented sites in all)")
        if len(cnt2) < len(two):
            return (f"explained at least partly by detection: {len(two)} genomes have >= 2 documented sites, but only {len(cnt2)} "
                    f"have >= 2 counted at their codon; see lineage_sensitivity.md for the miss reasons")
        return "not explained by detection: genomes with >= 2 counted sites exist but did not reach alpha"
    expl = {"P1": explain("CTG yeasts", set()) if len(p1) < 5 else "",
            "P2": explain("glass sponges", {"NC_010769.1"}) + "; per sponge: " + "; ".join(
                f"{meta[a_]['organism']}: {v_}" for a_, v_ in getattr(removed_breaks, "verdicts", {}).items()),
            "P3": explain("CTG yeasts", set(NAMED)),
            "P4": "not a sensitivity question: the prediction overlooked documented second AGT sites in some birds and turtles"}
    lines = ["# Predictions (section 1.15), scored by the criteria fixed in the PREREG log before results", "",
             f"Genomes scored: {len(gs)}."]
    for name, ok, detail in P:
        tag = name.split()[0]
        lines.append(f"- **{name}: {'PASS' if ok else 'FAIL'}.** {detail}"
                     + (f" Does sensitivity explain it? {expl[tag]}." if (not ok and expl.get(tag)) else ""))
    (out / "predictions.md").write_text("\n".join(lines) + "\n")

    md = ["# Flag context (added after the freeze; descriptive only, does not change which flags are strong)", ""]
    for r in rows:
        md.append(f"## {r['organism']} ({r['accession']}; {r['family']}, {r['order']}) {r['flag']}"
                  f"{' STRONG' if r['strong'] == '1' else ''}{' [known lineage]' if r['known_lineage'] == '1' else ''}")
        md.append(f"- p = {float(r['p']):.1e}, E = {float(r['E']):.2g}")
        md.append(f"- sites at the codon: {r['sites_at_codon']} ({len(r['sites_at_codon'].split('; '))} sites in {r['n_genes_at_codon']} genes; "
                  f"all in different genes: {'yes' if r['different_genes'] else 'no'})")
        md.append(f"- clean shifts in the genome: {r['clean_shifts_total']} (+1: {r['plus1']}, -1: {r['minus1']})")
        md.append(f"- platform: {r['platform'] or 'not given'}; assembly: {r['assembly_method'] or 'not given'}"
                  f"{'; source ' + r['source_accession'] if r['source_accession'] else ''}{'; SRA ' + r['sra_run'] if r['sra_run'] else ''}")
        md.append(f"- related genomes with a clean shift at the same codon and direction: {r['related_same_codon']}")
        md.append(f"- this genome's own chance of a flag, given its {r['n_clean_flag_direction']} clean shifts in this direction: "
                  f"{r['own_chance_mc']:.4f} (planted-null Monte Carlo, 20,000 draws); analytic upper bound {r['own_chance_upper_bound']:.4f}")
        if r["same_assembly_records"]:
            md.append(f"- other records of the same assembly/strain (one unit; cannot count as related): {r['same_assembly_records']}")
        md.append(f"- independent assemblies: {r['second_assembly']}")
        if r["known_lineage"] != "1" and r["congeners_in_scan"]:
            md.append(f"- same genus in the scan: {r['congeners_in_scan']}")
        if r["site_sources"]:
            md.append("- site sources (paper from PRIOR_ART only; record = RefSeq annotation):")
            md += [f"  - {x}" for x in r["site_sources"].split(" || ")]
        md.append("")
    kt = ["", "# Known-lineage flags: is each site described in a paper or only annotated in a record?", "",
          "| genome | flag | site | paper (PRIOR_ART only; reading depth) | record annotation |", "|---|---|---|---|---|"]
    for r in rows:
        if r["known_lineage"] == "1":
            for x in r["site_sources"].split(" || "):
                site, rest = x.split(": paper: ", 1)
                paper, rec = rest.split("; record: ", 1)
                kt.append(f"| {r['organism']} ({r['accession']}) | {r['flag']} | {site} | {paper} | {rec} |")
    md += kt
    (out / "flag_context.md").write_text("\n".join(md) + "\n")
    print("\n".join(md + lines))


if __name__ == "__main__":
    main()
