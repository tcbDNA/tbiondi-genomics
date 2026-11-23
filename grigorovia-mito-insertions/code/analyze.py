"""Analysis of the frozen scan (rules in PREREG.md section 1). Written before the scan was run.

    uv run python code/analyze.py --scan results/scan/main results/scan/fungal --out results/scan/analysis

Outputs:
  genome_flags.tsv      confirmatory single-genome flags (flags.genome_tests, window rule, fixed alpha),
                        with E (expected chance results at least this extreme across the whole scan) and
                        the strong-flag marks
  pool_tests.tsv        family and order pool tests (distinct orthologous sites), with shuffled-pool null
  flags_summary.md      per the stop rule: known lineages that flagged; every flag not in a known lineage;
                        which are strong; expected vs observed chance flags
  genome_coverage.tsv   per genome: core genes hit, aligned nt, masked share, core-CDS coverage (RefSeq)
  coverage_by_phylum.tsv
  tallies_*.tsv         exploratory tallies: bypass-like inserts, possible code change (stops read in frame),
                        detours around stops, homopolymer-killed shifts (direction x base)

Strong flag (DESIGN.md amendment 4): expected by chance less than once across the whole scan (E < 1), or
two related genomes (same family or order pool) flag the same codon in the same direction. Records of one
assembly or strain (code/units.py) count as one genome and never satisfy the second condition (PREREG log). Only strong
flags outside known lineages are funnel candidates, and the funnel is not started without asking.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fscore as F  # noqa: E402
import flags as FL  # noqa: E402

ROOT = F.ROOT
KNOWN_TAXA = {
    # prior-art lineages (PRIOR_ART.md, results/audit/controls_known_cases.tsv): family or order pool names
    "Aphrocallistidae", "Euplectellidae", "Rossellidae", "Caulophacidae", "Hyalonematidae", "Sceptrulophora",
    "Lyssacinosida", "Amphidiscosida", "Perkinsidae", "Perkinsida", "Saccharomycetaceae", "Dipodascaceae",
    "Testudines", "Veneridae", "Arcticidae", "Ostreidae", "Anomiidae", "Placunidae", "Plicatulidae",
    "Formicidae", "Phaeodactylaceae", "Rhizophydiaceae", "Spizellomycetaceae", "Physaraceae",
}
KNOWN_CLASSES = {"Aves", "Hexactinellida"}  # birds: nad3-174; glass sponges


def poisson_binomial_sf(k: int, qs: list[float]) -> float:
    """P(sum of independent Bernoulli(q_i) >= k)."""
    dist = [1.0]
    for q in qs:
        nd = [0.0] * (len(dist) + 1)
        for i, p in enumerate(dist):
            nd[i] += p * (1 - q)
            nd[i + 1] += p * q
        dist = nd
    return sum(dist[k:]) if k < len(dist) else 0.0


def load(scan_dirs):
    ev, gs, ins, stops = [], {}, [], []
    for d in scan_dirs:
        d = Path(d)
        ev += list(csv.DictReader(open(d / "events.tsv"), delimiter="\t"))
        for g in csv.DictReader(open(d / "genomes.tsv"), delimiter="\t"):
            gs[g["accession"]] = g
        if (d / "inserts.tsv").exists():
            ins += list(csv.DictReader(open(d / "inserts.tsv"), delimiter="\t"))
        if (d / "stop_columns.tsv").exists():
            stops += list(csv.DictReader(open(d / "stop_columns.tsv"), delimiter="\t"))
    return ev, gs, ins, stops


def is_known(meta: dict) -> bool:
    return (meta["family_pool_name"] in KNOWN_TAXA or meta["order_pool_name"] in KNOWN_TAXA
            or meta["class_name"] in KNOWN_CLASSES)


def pool_sites(events, cap):
    """Distinct orthologous sites in a pool: same gene and direction, subject positions within
    POOL_SITE_RESIDUES (50) residues of the site's first shift, as wide as the placement error seen in
    annotated shifts (PREREG log 07:10; was relative position within 2%).
    A site votes for codon c if at least half of its member shifts vote for c (window rule)."""
    sites = []
    for e in sorted(events, key=lambda e: (e["gene"], int(e["shift"]), int(e["spos"]))):
        sp = int(e["spos"])
        for s in sites:
            if s["gene"] == e["gene"] and s["shift"] == int(e["shift"]) and abs(s["spos"] - sp) <= FL.POOL_SITE_RESIDUES:
                break
        else:
            s = {"gene": e["gene"], "shift": int(e["shift"]), "spos": sp, "votes": Counter(), "n": 0,
                 "members": set(), "first": e["accession"]}
            sites.append(s)
        s["n"] += 1
        s["members"].add(e["accession"])
        for c in FL.votes(e, cap):
            s["votes"][c] += 1
    for s in sites:
        s["codons"] = {c for c, v in s["votes"].items() if v >= 0.5 * s["n"]}
    return sites


def pool_tests(clean_by_acc, meta, gc, null, rank_key, alpha_pool):
    """Returns (rows, n_tested_pool_dirs). Flag if k >= 3 sites in >= 2 genes and Poisson-binomial tail
    <= alpha_pool / 4 (kept from section 0 although only one codon rule remains: conservative)."""
    pools = defaultdict(list)
    for a in clean_by_acc:
        pools[meta[a][rank_key]].append(a)
    out, tested = [], 0
    for pool, members in pools.items():
        evs = [e for a in members for e in clean_by_acc[a]]
        if not evs:
            continue
        sites = pool_sites(evs, null.cap)
        for d in (1, -1):
            sd = [s for s in sites if s["shift"] == d]
            if len(sd) < 3:
                continue
            tested += 1
            for codon in set().union(*(s["codons"] for s in sd)):
                hit = [s for s in sd if codon in s["codons"]]
                k = len(hit)
                if k < 3 or len({s["gene"] for s in hit}) < 2:
                    continue
                qs = [null.q(int(meta[s["first"]]["code"]), gc[s["first"]], d, codon)[0] for s in sd]
                p = poisson_binomial_sf(k, qs)
                out.append({"rank": rank_key, "pool": pool, "shift": d, "codon": codon, "k": k,
                            "n_sites": len(sd), "genes": len({s['gene'] for s in hit}), "members": len(members),
                            "p": p, "flag": int(p <= alpha_pool / 4)})
    return out, tested


def scan_E(p_obs: float, groups: Counter, null) -> float:
    """Expected number of (genome, direction, codon) results with p <= p_obs across the whole scan if all
    clean shifts were random, from the scan's own counts. groups: (n, null key) -> genome-directions."""
    tot = 0.0
    for (n, key), cnt in groups.items():
        c_ = null.counts.get(key, Counter())
        N = null.totals.get(key, 0)
        for codon in FL.CODONS:
            q = (c_.get(codon, 0) + 0.5) / (N + 32.0)
            k = FL.kmin(n, q, p_obs)
            if k is not None:
                tot += cnt * FL.binom_sf(k, n, q)
    return tot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha-pool", type=float, default=1e-5)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}
    null = FL.Null(ROOT / os.environ.get("CTRL_OUT", "results/controls") / "plant_null/null.json")
    ev, gs, ins, stops = load(args.scan)
    gc = {a: float(g["gc"]) for a, g in gs.items()}
    by_acc = defaultdict(list)
    for e in ev:
        by_acc[e["accession"]].append(e)
    stops_of = {}

    # ---- confirmatory single-genome test + the scan's own chance expectation
    gflags, groups, exp_flags = [], Counter(), 0.0
    for a in gs:
        code = int(meta[a]["code"])
        stops_c = stops_of.setdefault(code, {c for c, x in F.codon_table(code).items() if x == "*"})
        evs = by_acc.get(a, [])
        for d, (n, qs, al_) in zip((1, -1), FL.genome_dir_inputs(evs, code, gc[a], null, stops_c)):
            if n >= 2:
                groups[(n, null.key(code, gc[a], d))] += 1
                exp_flags += FL.expected_chance(n, qs, al_)
        for r in FL.genome_tests(evs, code, gc[a], null, stops_c):
            if r["flag"]:
                gflags.append({"accession": a, "organism": meta[a]["organism"], "family": meta[a]["family_pool_name"],
                               "order": meta[a]["order_pool_name"], "phylum": meta[a]["phylum_name"],
                               "family_pool": meta[a]["family_pool"], "order_pool": meta[a]["order_pool"],
                               "known_lineage": int(is_known(meta[a])), **r})
    import units as U  # noqa: E402  (PREREG log: records of one assembly/strain are one unit)
    unit = U.unit_of()
    for f in gflags:
        f["unit"] = unit.get(f["accession"], f["accession"])
    for f in gflags:
        f["E"] = scan_E(f["p"], groups, null)
        rel = [g for g in gflags if g is not f and g["unit"] != f["unit"] and g["codon"] == f["codon"] and g["shift"] == f["shift"]
               and (g["family_pool"] == f["family_pool"] or g["order_pool"] == f["order_pool"])]
        f["related_same_codon"] = ",".join(sorted({g["accession"] for g in rel}))
        f["strong"] = int(f["E"] < 1 or bool(rel))
    write(out / "genome_flags.tsv", gflags)

    # ---- lineage (pool) test + shuffled-pool null
    clean = {a: FL.clean_core(by_acc.get(a, []), null.split_votes) for a in gs}
    pflags, n_tested = [], 0
    for rank in ("family_pool", "order_pool"):
        rows, t = pool_tests(clean, meta, gc, null, rank, args.alpha_pool)
        n_tested += t
        for r in rows:
            any_member = next(a for a in gs if meta[a][rank] == r["pool"])
            r["pool_name"] = meta[any_member][rank + "_name"]
            r["known_lineage"] = int(is_known(meta[any_member]))
            pflags.append(r)
    for r in pflags:
        r["E"] = r["p"] * 64 * n_tested
        r["strong"] = int(r["flag"] and r["E"] < 1)
    write(out / "pool_tests.tsv", pflags)
    shuffled = []
    for seed in (1, 2, 3):
        rng = random.Random(seed)
        for rank in ("family_pool", "order_pool"):
            accs = list(gs)
            labels = [meta[a][rank] for a in accs]
            rng.shuffle(labels)
            fake = {a: dict(meta[a], **{rank: lab}) for a, lab in zip(accs, labels)}
            n = sum(r["flag"] for r in pool_tests(clean, fake, gc, null, rank, args.alpha_pool)[0])
            shuffled.append({"seed": seed, "rank": rank, "flags": n})
    write(out / "shuffled_pools.tsv", shuffled)

    # ---- exploratory tallies
    byp = Counter((meta[i["accession"]]["family_pool_name"], meta[i["accession"]]["phylum_name"]) for i in ins if i["class"] == "bypass_like")
    write(out / "tally_bypass_like.tsv", [{"family": f, "phylum": p, "n": n} for (f, p), n in byp.most_common()])
    cc = Counter((meta[s["accession"]]["family_pool_name"], s["subject_aa"]) for s in stops)
    write(out / "tally_stop_read_in_frame.tsv", [{"family": f, "subject_residue": aa, "n": n} for (f, aa), n in cc.most_common()])
    det = Counter(meta[e["accession"]]["family_pool_name"] for e in ev if e["detour_around_stop"] == "1")
    write(out / "tally_detour_around_stop.tsv", [{"family": f, "n_events": n} for f, n in det.most_common()])
    hp = Counter((meta[e["accession"]]["family_pool_name"], e["shift"], e["homopolymer_base"]) for e in ev if e["class"] == "homopolymer")
    write(out / "tally_homopolymer.tsv", [{"family": f, "shift": s, "base": b, "n": n} for (f, s, b), n in hp.most_common()])
    # genomes with >= 3 clean shifts in different core genes, whatever the codons (PREREG log 07:10)
    multi = []
    for a in gs:
        cl = clean[a]
        genes = sorted({e["gene"] for e in cl})
        if len(genes) >= 3:
            multi.append({"accession": a, "organism": meta[a]["organism"], "family": meta[a]["family_pool_name"],
                          "order": meta[a]["order_pool_name"], "phylum": meta[a]["phylum_name"],
                          "known_lineage": int(is_known(meta[a])), "n_clean": len(cl), "n_genes": len(genes),
                          "genes": ",".join(genes),
                          "shifts": ";".join(f"{e['gene']}:{int(e['shift']):+d}:{e['codon_left']}/{e['codon']}" for e in cl)})
    write(out / "tally_three_gene_genomes.tsv", sorted(multi, key=lambda r: -r["n_genes"]))
    # split-hit shifts, when they did not pass probation: tally only
    if not null.split_votes:
        sp = [e for e in ev if e["via"] == "split" and e["class"] == "clean" and e["gene"] in FL.CORE]
        write(out / "tally_split_hit_shifts.tsv", [{"accession": e["accession"], "gene": e["gene"], "shift": e["shift"],
                                                    "codons": f"{e['codon_left']}/{e['codon']}", "window": e["window"]} for e in sp])

    # ---- coverage: per genome and by phylum
    aligned = {}
    for d in args.scan:
        p = Path(d) / "aligned_fwd.json"
        if p.exists():
            aligned.update(json.load(open(p)))
    core_cov = {}
    cov = defaultdict(lambda: [0, 0])
    per = defaultdict(lambda: [0, 0])
    for r in csv.DictReader(open(ROOT / "data/work/cds.tsv"), delimiter="\t"):
        a = r["accession"]
        if a not in gs or r["gene"] not in FL.CORE:
            continue
        parts = [tuple(map(int, x.split("-"))) for x in r["parts"].split(";")]
        L = sum(y - x for x, y in parts)
        ov = min(L, sum(max(0, min(y, b) - max(x, a_)) for x, y in parts for a_, b, _ in aligned.get(a, [])))
        ph = meta[a]["phylum_name"] or meta[a]["kingdom_name"] or "unranked"
        cov[ph][0] += 1
        cov[ph][1] += ov >= 0.5 * L
        per[a][0] += L
        per[a][1] += ov
    for a, (L, ov) in per.items():
        core_cov[a] = round(ov / L, 3) if L else ""
    write(out / "genome_coverage.tsv", [{"accession": a, "set": meta[a]["set"], "phylum": meta[a]["phylum_name"],
                                         "core_genes_hit": len(set(g["genes_hit"].split(",")) & FL.CORE),
                                         "aligned_nt": g["aligned_nt"], "masked_share": g["masked_share"],
                                         "core_cds_coverage": core_cov.get(a, "")} for a, g in gs.items()])
    write(out / "coverage_by_phylum.tsv", sorted(({"phylum": p, "core_cds": n, "aligned_share": round(k / n, 3)} for p, (n, k) in cov.items()), key=lambda r: -r["core_cds"]))

    # ---- summary, by the stop rule
    def gline(f):
        return (f"- {f['organism']} ({f['accession']}; {f['family']}, {f['order']}, {f['phylum']}): {f['shift']:+d} {f['codon']} "
                f"k={f['k']}/{f['n']}, genes {f['genes']}, p={f['p']:.1e}, E={f['E']:.2g}"
                + (f", same codon in related {f['related_same_codon']}" if f["related_same_codon"] else "")
                + (" **STRONG**" if f["strong"] else ""))
    known_g = [f for f in gflags if f["known_lineage"]]
    new_g = [f for f in gflags if not f["known_lineage"]]
    pf = [r for r in pflags if r["flag"]]
    cand_g = [f for f in new_g if f["strong"]]
    cand_p = [r for r in pf if r["strong"] and not r["known_lineage"]]
    cov_vals = sorted(v for v in core_cov.values() if v != "")
    lines = ["# Scan flags", "",
             f"Genomes scanned: {len(gs)}. Window cap {null.cap} nt; alpha_sense = {FL.ALPHA_SENSE:g}, alpha_stop = {FL.ALPHA_STOP:g}; "
             f"alpha_pool = {args.alpha_pool:g} (/4).",
             f"Genome-directions with >= 2 clean core shifts: {sum(groups.values())}. Expected chance flags under the rule, from the "
             f"scan's own counts: {exp_flags:.1f} (upper bound). Observed single-genome flags: {len(gflags)}.",
             f"Shuffled-pool flags (seeds 1-3; family, order): {[s['flags'] for s in shuffled]}. Pool-direction tests: {n_tested}.",
             f"Split-hit shifts {'vote' if null.split_votes else 'do not vote (tally only)'}. Exploratory: genomes with clean shifts "
             f"in >= 3 different core genes: {len(multi)} ({sum(m['known_lineage'] for m in multi)} in known lineages).",
             f"Core-CDS coverage, RefSeq genomes with CDS: median {cov_vals[len(cov_vals)//2] if cov_vals else 'NA'}, "
             f"10th pct {cov_vals[len(cov_vals)//10] if cov_vals else 'NA'}.", "",
             f"## Single-genome flags in known lineages: {len(known_g)}"]
    lines += [gline(f) for f in known_g]
    lines += ["", f"## Single-genome flags NOT in a known lineage: {len(new_g)} ({len(cand_g)} strong)"]
    lines += [gline(f) for f in new_g]
    lines += ["", f"## Pool flags: {len(pf)} ({sum(r['known_lineage'] for r in pf)} in known lineages)"]
    lines += [f"- {r['rank']} {r['pool_name']} ({'known' if r['known_lineage'] else 'NOT known'}): {r['shift']:+d} {r['codon']} "
              f"k={r['k']}/{r['n_sites']} sites, genes {r['genes']}, members {r['members']}, p={r['p']:.1e}, E={r['E']:.2g}"
              + (" **STRONG**" if r["strong"] else "") for r in pf]
    lines += ["", "## Stop rule",
              f"Strong flags outside known lineages (funnel candidates; funnel not started without asking): "
              f"{len(cand_g)} genome, {len(cand_p)} pool."]
    if not cand_g and not cand_p:
        lines.append("None: write the short negative with the measured sensitivity and stop the project.")
    (out / "flags_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
