"""Controls before data: planted null, null genomes with sequencing-error copies, planted same-codon pairs,
recall of RefSeq-annotated shifts, named positives, and the alignment-based intron calibration.

    uv run python code/controls.py <stage> [...]   stages: plant_null null_genomes sites annotated named introns
                                                  intron_summary bars
    uv run python code/controls.py all

Every stage runs genomes through scan.scan / scan.process_genome, the same code as the real scan.

Re-run after the user's decisions (DESIGN.md amendment 4):
  - one window rule (flags.py); the window cap is set here from the planted data: the smallest width that
    includes >= WINDOW_SHARE of clean recovered planted shifts;
  - alpha fixed at flags.ALPHA_SENSE / ALPHA_STOP (the values the pre-written rule gave); `bars` also
    reports what that rule would give on this run, for the record only;
  - planted same-codon trials with 2, 3 and 5 sites (stage `sites`), with each site's position kept so
    per-site detection and the reasons for misses can be reported;
  - B9a and B9c are stated limits, not bars.
`bars` writes results/controls/BARS.md and exits 2 if a scan gate fails (null-genome flags > 0,
B9b > 2%, a named positive not flagging), 1 if only other bars are missed (logged as limits), else 0.

Pre-written threshold rule (first run; kept for the record):
  alpha grid = 1e-6 .. 1e-2 in half-decade steps.
  alpha_sense = smallest grid value at which A. vastus flags AND >= 90% of planted sense-codon pairs flag.
  alpha_stop  = smallest grid value at which >= 90% of planted stop-codon pairs flag.
"""

from __future__ import annotations

import csv
import json
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fscore as F  # noqa: E402
import flags as FL  # noqa: E402
import scan as S  # noqa: E402

ROOT = F.ROOT
import os  # noqa: E402
OUT = ROOT / os.environ.get("CTRL_OUT", "results/controls")
LIMIT = int(os.environ.get("CTRL_LIMIT", "1000000000"))
CTRL = ROOT / "protocol/controls"
GRID = [1e-6, 3.2e-6, 1e-5, 3.2e-5, 1e-4, 3.2e-4, 1e-3, 3.2e-3, 1e-2]
WINDOW_SHARE = 0.95
TRIALS = [("sense", 2), ("sense", 2), ("sense", 3), ("sense", 5), ("stop", 2), ("stop", 3), ("stop", 5)]
NAMED = {"NC_010769.1": "Aphrocallistes vastus", "NC_068742.1": "Perkinsus marinus", "NC_068743.1": "Perkinsus olseni",
         "NC_036748.1": "Huiozyma sinensis (CTGY)", "NC_004918.1": "Monosporozyma servazzii (CTGY)",
         "NC_024095.1": "Magnusiomyces capitatus (bypass)"}


# ---------------------------------------------------------------- helpers
def meta():
    return {r["accession"]: r for r in csv.DictReader(open(ROOT / "data/work/genomes.tsv"), delimiter="\t")}


def seqs_for(accs, G):
    by = defaultdict(set)
    for a in accs:
        by[G[a]["code"]].add(a)
    out = {}
    for code, s in by.items():
        out.update(S.read_fasta(ROOT / f"data/work/queries/code_{code}.fna", s))
    return out


def gc(s: str) -> float:
    return (s.count("G") + s.count("C")) / max(1, len(s) - s.count("N"))


def cds_parts():
    out = defaultdict(list)
    for r in csv.DictReader(open(ROOT / "data/work/cds.tsv"), delimiter="\t"):
        parts = [tuple(map(int, p.split("-"))) for p in r["parts"].split(";")]
        out[r["accession"]].append((r["gene"], int(r["strand"]), parts))
    return out


def hp_near(seq: str, p: int) -> bool:
    return F.homopolymer(seq, p, p)[0]


def apply_edits(seq: str, edits: list[tuple[int, str, str]]) -> tuple[str, list[int]]:
    """edits: (pos, 'ins'|'del', base). Returns mutated seq and each edit's position in the mutated seq."""
    order = sorted(range(len(edits)), key=lambda i: edits[i][0])
    s = list(seq)
    for i in sorted(order, key=lambda i: -edits[i][0]):
        p, kind, b = edits[i]
        if kind == "ins":
            s.insert(p, b)
        else:
            del s[p]
    shift, newpos = 0, [0] * len(edits)
    for i in order:
        p, kind, _ = edits[i]
        newpos[i] = p + shift
        shift += 1 if kind == "ins" else -1
    return "".join(s), newpos


def match(events: list[dict], acc: str, pos: int, d: int, reach: int = 9):
    cands = [e for e in events if e["accession"] == acc and int(e["shift"]) == d and abs(int(e["pos_fwd"]) - pos) <= reach]
    if not cands:
        return None
    return min(cands, key=lambda e: (e["class"] != "clean", abs(int(e["pos_fwd"]) - pos)))


def run(accs, seqs, G, out_dir, code_of=None, family_of=None):
    code_of = code_of or {a: int(G[a]["code"]) for a in accs}
    family_of = family_of or {a: G[a]["family_pool"] for a in accs}
    return S.scan(accs, out_dir, 12, seqs=seqs, code_of=code_of, family_of=family_of)


def save(res, d: Path):
    S.write_tables(res, d)
    with open(d / "aligned_fwd.json", "w") as fh:
        json.dump({r["summary"]["accession"]: r["aligned_fwd"] for r in res}, fh)


def load_events(d: Path):
    return list(csv.DictReader(open(d / "events.tsv"), delimiter="\t")) if (d / "events.tsv").exists() else []


def in_aligned(al: list, pos: int) -> bool:
    return any(a + 45 <= pos <= b - 45 for a, b, _ in al)


def covered(al: list, pos: int) -> bool:
    return any(a <= pos < b for a, b, _ in al)


def group_of(meta: dict) -> str:
    k = meta["kingdom_name"]
    return "animals" if k == "Metazoa" else "fungi" if k == "Fungi" else "protists"


def site_outcome(m, al, pos: int, codon: str, cap: int) -> str:
    """Why a planted site did or did not count: 'votes' = clean and voting for the planted codon."""
    if m is None:
        return "aligned_no_event" if covered(al, pos) else "not_aligned"
    if m["class"] != "clean":
        return m["class"]
    v = FL.votes(m, cap)
    if not v:
        return "clean_wide_window"
    return "votes" if codon in v else "clean_other_codon"


# ---------------------------------------------------------------- stage: planted null + single-shift recall
def stage_plant_null():
    G, CDS = meta(), cds_parts()
    accs = [l.strip() for l in open(CTRL / "plant_genomes.txt") if l.strip()][:LIMIT]
    seqs = seqs_for(accs, G)
    rng = random.Random(101)
    copies, code_of, family_of, plants = {}, {}, {}, []
    for a in accs:
        s = seqs[a]
        comp = [s.count(b) for b in "ACGT"]
        for c in range(3):
            edits, info = [], []
            for gene, strand, parts in CDS[a]:
                if gene not in FL.CORE:
                    continue
                cand = [p for (x, y) in parts for p in range(x + 60, y - 60)]
                rng.shuffle(cand)
                for p in cand[:40]:
                    if hp_near(s, p) or any(abs(p - q) < 90 for q, _, _ in edits):
                        continue
                    kind = rng.choice(["ins", "del"])
                    base = rng.choices("ACGT", weights=comp)[0]
                    edits.append((p, kind, base))
                    info.append({"genome": a, "gene": gene, "strand": strand, "orig_pos": p, "kind": kind,
                                 "_parts": parts})
                    break
            if not edits:
                continue
            m, newpos = apply_edits(s, edits)
            cid = f"{a}__p{c}"
            copies[cid], code_of[cid], family_of[cid] = m, int(G[a]["code"]), G[a]["family_pool"]
            for i, np_ in zip(info, newpos):
                plants.append({**i, "copy": cid, "pos": np_, "shift": 1 if i["kind"] == "ins" else -1, "_edits": edits})
    d = OUT / "plant_null"
    res = run(list(copies), copies, G, d / "scan", code_of, family_of)
    save(res, d / "scan")
    al = {r["summary"]["accession"]: r["aligned_fwd"] for r in res}
    ev = [e for r in res for e in r["events"]]
    by_acc = defaultdict(list)
    for e in ev:
        by_acc[e["accession"]].append(e)
    gc_of = {a: gc(seqs[a]) for a in accs}
    for p in plants:
        evs = by_acc[p["copy"]]
        m = match(evs, p["copy"], p["pos"], p["shift"])
        mh = match([e for e in evs if e["via"] != "split"], p["copy"], p["pos"], p["shift"])  # single-hit path only
        p["aligned"] = int(in_aligned(al.get(p["copy"], []), p["pos"]))
        mseq = copies[p["copy"]]
        # "homopolymers excluded": judged on the mutated sequence with the same kill rule the scan uses
        p["hp_context"] = int(F.homopolymer(mseq, p["pos"] - 2, p["pos"] + 2, p["shift"])[0])
        for suf, mm in (("", m), ("_h", mh)):
            p["found" + suf] = int(mm is not None)
            p["clean" + suf] = int(mm is not None and mm["class"] == "clean")
            p["class" + suf] = mm["class"] if mm else ""
            p["codon_R" + suf] = mm["codon"] if mm else ""
            p["codon_L" + suf] = mm["codon_left"] if mm else ""
            p["window" + suf] = mm["window"] if mm else ""
        p["via"] = m["via"] if m else ""
        p["clean_split_only"] = int(p["clean"] == 1 and m["via"] == "split" and p["clean_h"] == 0)
        p["in_gene"] = int(any(int(e["shift"]) == p["shift"] and e["gene"] == p["gene"] for e in evs))
        p["in_gene_h"] = int(any(int(e["shift"]) == p["shift"] and e["gene"] == p["gene"] and e["via"] != "split" for e in evs))
        # the site's own codons: what trigger() gives for a shift placed exactly at the planted position
        L_ = len(mseq)
        so = mseq if p["strand"] == 1 else mseq.translate(S.COMP)[::-1]
        if p["shift"] == 1:
            opos = p["pos"] if p["strand"] == 1 else L_ - 1 - p["pos"]
        else:
            opos = p["pos"] if p["strand"] == 1 else L_ - p["pos"]
        def mut(x):  # original forward coordinate -> mutated copy coordinate
            return x + sum((1 if k_ == "ins" else -1) for q_, k_, _ in p["_edits"] if q_ < x)
        f0 = F.old_frame0([(mut(x), mut(y)) for x, y in p["_parts"]], p["strand"], L_, opos)
        tr = F.trigger(F.Event(p["copy"], "", p["strand"], opos, p["shift"], p["shift"], 0, 0), so,
                       F.codon_table(int(G[p["genome"]]["code"])), frame0=f0)
        p["true_R"], p["true_L"] = tr["codon"], tr["codon_left"]
        p["covered"] = int(covered(al.get(p["copy"], []), p["pos"]))
        p["group"] = group_of(G[p["genome"]])
        p["code"] = G[p["genome"]]["code"]
        p["gc"] = round(gc_of[p["genome"]], 4)
    for p in plants:
        p.pop("_parts", None)
        p.pop("_edits", None)
    with open(d / "plants.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(plants[0]), delimiter="\t")
        w.writeheader()
        w.writerows(plants)
    # null counts by stratum (code x GC half), from CLEAN recovered planted shifts, for both split regimes
    gc_cut = {}
    for code in {p["code"] for p in plants}:
        vals = sorted(gc_of[a] for a in accs if G[a]["code"] == code)
        gc_cut[code] = vals[len(vals) // 2] if len(vals) >= 20 else None
    variants = {}
    for name, suf in (("split", ""), ("nosplit", "_h")):
        clean = [p for p in plants if p["clean" + suf]]
        ws = sorted(int(p["window" + suf]) for p in clean)
        # window cap: smallest width including >= WINDOW_SHARE of clean recovered planted shifts
        cap = next(w for w in range(ws[-1] + 1) if sum(x <= w for x in ws) >= WINDOW_SHARE * len(ws))
        counts, totals = defaultdict(Counter), Counter()
        for p in clean:
            cut = gc_cut.get(p["code"])
            strata = [f"{p['code']}|all", "ALL"]
            if cut is not None:
                strata.append(f"{p['code']}|{'lo' if p['gc'] < cut else 'hi'}")
            v = FL.votes({"window": p["window" + suf], "codon": p["codon_R" + suf], "codon_left": p["codon_L" + suf]}, cap)
            for st in strata:
                k = f"{st}|{p['shift']}"
                totals[k] += 1
                for c in v:
                    counts[k][c] += 1
        variants[name] = {"counts": counts, "totals": totals, "window_cap": cap, "n_clean": len(clean),
                          "window_dist": dict(sorted(Counter(ws).items()))}
    json.dump({"variants": variants, "gc_cut": {k: v for k, v in gc_cut.items() if v is not None},
               "min_events": 500, "window_share": WINDOW_SHARE}, open(d / "null.json", "w"), indent=0)
    print(f"plant_null: {len(plants)} planted, {sum(p['found'] for p in plants)} found, "
          f"{variants['split']['n_clean']} clean ({variants['nosplit']['n_clean']} single-hit); window cap "
          f"{variants['split']['window_cap']} / {variants['nosplit']['window_cap']}")


# ---------------------------------------------------------------- stage: null genomes + error copies
def stage_null_genomes():
    G = meta()
    accs = [l.strip() for l in open(CTRL / "null_genomes.txt") if l.strip()][:LIMIT]
    seqs = seqs_for(accs, G)
    rng = random.Random(202)
    copies, code_of, family_of = {}, {}, {}
    for a in accs:
        s = seqs[a]
        copies[a] = s
        code_of[a], family_of[a] = int(G[a]["code"]), G[a]["family_pool"]
        runs = [(m.start(), m.end(), m.group(0)[0]) for m in re.finditer(r"(A{5,}|C{5,}|G{5,}|T{5,})", s)]
        n_hp = min(len(runs), max(1, round(0.5 * len(s) / 1000)))
        edits = []
        for st, en, b in rng.sample(runs, n_hp):
            edits.append((st + 1, "del" if rng.random() < 0.7 else "ins", b))
        copies[a + "__hp"], code_of[a + "__hp"], family_of[a + "__hp"] = apply_edits(s, edits)[0], code_of[a], family_of[a]
        n_rnd = max(1, round(0.1 * len(s) / 1000))
        edits = []
        while len(edits) < n_rnd:
            p = rng.randrange(10, len(s) - 10)
            if not hp_near(s, p):
                edits.append((p, rng.choice(["ins", "del"]), rng.choice("ACGT")))
        copies[a + "__rnd"], code_of[a + "__rnd"], family_of[a + "__rnd"] = apply_edits(s, edits)[0], code_of[a], family_of[a]
    d = OUT / "null_genomes"
    res = run(list(copies), copies, G, d / "scan", code_of, family_of)
    save(res, d / "scan")
    json.dump({k: round(gc(v), 4) for k, v in copies.items()}, open(d / "gc.json", "w"))
    print(f"null_genomes: {len(copies)} sequences scanned")


# ---------------------------------------------------------------- stage: 2, 3 or 5 shifts planted at one codon
def stage_sites():
    G, CDS = meta(), cds_parts()
    null = FL.Null(OUT / "plant_null/null.json", split_votes=True)  # only to rank codons by rarity
    accs = [l.strip() for l in open(CTRL / "null_genomes.txt") if l.strip()][:LIMIT]
    seqs = seqs_for(accs, G)
    rng = random.Random(303)
    copies, code_of, family_of, trials = {}, {}, {}, []
    for a in accs:
        s, code = seqs[a], int(G[a]["code"])
        table = F.codon_table(code)
        stops = sorted(c for c, x in table.items() if x == "*")
        g = gc(s)
        qs = {c: null.q(code, g, 1, c)[0] for c in FL.CODONS if c not in stops}
        rare = sorted(qs, key=qs.get)[: len(qs) // 4]
        # in-frame codon starts in single-part core-gene CDS; insertion points in forward coordinates
        sites = defaultdict(list)  # (gene, strand) -> list of (fwd insertion point, codon)
        for gene, strand, parts in CDS[a]:
            if gene not in FL.CORE or len(parts) != 1:
                continue
            x, y = parts[0]
            if strand == 1:
                for p in range(x + 60, y - 60, 3):
                    sites[(gene, 1)].append((p, s[p:p + 3]))
            else:
                rc = s[x:y].translate(S.COMP)[::-1]
                L = y - x
                for o in range(60, L - 60, 3):
                    sites[(gene, -1)].append((y - o, rc[o:o + 3]))
        for ti, (kind, nsites) in enumerate(TRIALS):
            pool = rare if kind == "sense" else stops
            for _ in range(30):
                c = rng.choice(pool)
                groups = [(gs, [p for p, cod in v if cod[:2] == c[1:]]) for gs, v in sites.items()]
                groups = [(gs, ps) for gs, ps in groups if ps]
                if len({gs[0] for gs, _ in groups}) < 2:
                    continue
                rng.shuffle(groups)
                chosen = []  # (fwd insertion point, gene, strand); one per gene first, then round-robin
                for _round in range(nsites):
                    for (gene, strand), ps in groups:
                        if len(chosen) >= nsites:
                            break
                        cand = [p for p in ps if all(abs(p - q) >= 90 for q, _, _ in chosen)]
                        if cand:
                            chosen.append((rng.choice(cand), gene, strand))
                    if len(chosen) >= nsites:
                        break
                if len(chosen) < nsites or len({gn for _, gn, _ in chosen}) < 2:
                    continue
                edits = [(p, "ins", c[0] if strand == 1 else c[0].translate(S.COMP)) for p, _, strand in chosen]
                m, newpos = apply_edits(s, edits)
                cid = f"{a}__{kind}{nsites}_{ti}"
                copies[cid], code_of[cid], family_of[cid] = m, code, G[a]["family_pool"]
                trials.append({"copy": cid, "genome": a, "kind": kind, "n_sites": nsites, "codon": c, "code": code,
                               "group": group_of(G[a]),
                               "sites": [{"gene": gn, "strand": st, "pos": np_} for (_, gn, st), np_ in zip(chosen, newpos)]})
                break
    d = OUT / "sites"
    res = run(list(copies), copies, G, d / "scan", code_of, family_of)
    save(res, d / "scan")
    json.dump(trials, open(d / "trials.json", "w"), indent=0)
    json.dump({k: round(gc(v), 4) for k, v in copies.items()}, open(d / "gc.json", "w"))
    print(f"sites: {len(trials)} trials {dict(Counter((t['kind'], t['n_sites']) for t in trials))}")


# ---------------------------------------------------------------- stage: RefSeq-annotated shifts
def stage_annotated():
    G = meta()
    J = list(csv.DictReader(open(ROOT / "results/audit/known_junctions_refseq238.tsv"), delimiter="\t"))
    accs = sorted({j["accession"] for j in J if j["accession"] in G and G[j["accession"]]["set"] == "refseq_main"})[:LIMIT]
    seqs = seqs_for(accs, G)
    d = OUT / "annotated"
    res = run(accs, seqs, G, d / "scan")
    save(res, d / "scan")
    print(f"annotated: {len(accs)} records scanned")


def stage_named():
    G = meta()
    accs = list(NAMED)
    seqs = seqs_for(accs, G)
    d = OUT / "named"
    res = run(accs, seqs, G, d / "scan")
    save(res, d / "scan")
    json.dump({a: round(gc(seqs[a]), 4) for a in accs}, open(d / "gc.json", "w"))
    print("named: done")


# ---------------------------------------------------------------- stage: intron calibration
def stage_introns():
    G = meta()
    accs = [l.strip() for l in open(CTRL / "intron_genomes.txt") if l.strip()][:LIMIT]
    seqs = seqs_for(accs, G)
    introns = defaultdict(list)
    want = set(accs)
    for r in csv.DictReader(open(ROOT / "data/work/introns.tsv"), delimiter="\t"):
        a, st, en = r["accession"], int(r["start"]), int(r["end"])
        if a in want and r["gene"] in FL.CORE and en - st > 30:
            if not any(abs(st - x) <= 5 and abs(en - y) <= 5 for x, y, _ in introns[a]):
                introns[a].append((st, en, r["gene"]))
    CDS = cds_parts()
    rng = random.Random(404)
    copies, code_of, family_of, plants = {}, {}, {}, []
    for a in accs:
        s = seqs[a]
        copies[a], code_of[a], family_of[a] = s, int(G[a]["code"]), G[a]["family_pool"]
        # plant near-edge shifts: 31..90 nt into the exon next to an annotated intron edge
        exon_parts = [(x, y) for gene, _, parts in CDS[a] if gene in FL.CORE for (x, y) in parts]
        edits, info = [], []
        for st, en, gene in rng.sample(introns[a], min(4, len(introns[a]))):
            for edge, side in ((st, -1), (en, 1)):
                dist = rng.randint(31, 90)
                p = edge + side * dist if side == 1 else edge - dist
                if not any(x + 3 <= p <= y - 3 for x, y in exon_parts) or hp_near(s, p):
                    continue
                if any(abs(p - q) < 120 for q, _, _ in edits):
                    continue
                kind = rng.choice(["ins", "del"])
                edits.append((p, kind, rng.choice("ACGT")))
                info.append({"genome": a, "gene": gene, "edge": edge, "dist": dist, "kind": kind})
                break
        if edits:
            m, newpos = apply_edits(s, edits)
            cid = a + "__edge"
            copies[cid], code_of[cid], family_of[cid] = m, code_of[a], family_of[a]
            for i, np_ in zip(info, newpos):
                plants.append({**i, "copy": cid, "pos": np_, "shift": 1 if i["kind"] == "ins" else -1})
    d = OUT / "introns"
    res = run(list(copies), copies, G, d / "scan", code_of, family_of)
    save(res, d / "scan")
    json.dump({"introns": introns, "plants": plants}, open(d / "truth.json", "w"))
    print(f"introns: {len(accs)} genomes, {sum(len(v) for v in introns.values())} annotated core-gene introns, {len(plants)} near-edge plants")


# ---------------------------------------------------------------- bars
def eligible(events, split_votes: bool):
    return events if split_votes else [e for e in events if e["via"] != "split"]


def flags_for(events, gc_of, codes, null, a_sense=FL.ALPHA_SENSE, a_stop=FL.ALPHA_STOP):
    by = defaultdict(list)
    for e in events:
        by[e["accession"]].append(e)
    out = {}
    for acc, code in codes.items():
        stops = {c for c, x in F.codon_table(code).items() if x == "*"}
        out[acc] = FL.genome_tests(by.get(acc, []), code, gc_of[acc], null, stops, a_sense, a_stop)
    return out


def intron_summary(G=None, verbose=True, split_votes=True) -> dict:
    """B9a, B9b, B9c, near-edge coverage, per-genome core coverage and masked share for the intron set."""
    G = G or meta()
    truth = json.load(open(OUT / "introns/truth.json"))
    iev = defaultdict(list)
    for e in eligible(load_events(OUT / "introns/scan"), split_votes):
        iev[e["accession"]].append(e)
    iins = defaultdict(list)
    if (OUT / "introns/scan/inserts.tsv").exists():
        for i in csv.DictReader(open(OUT / "introns/scan/inserts.tsv"), delimiter="\t"):
            iins[i["accession"]].append(i)
    al = json.load(open(OUT / "introns/scan/aligned_fwd.json"))
    n_int = n_inf = n_edge = n_edge_fs = 0
    edge_fs_split = 0
    for a, ints in truth["introns"].items():
        glen = int(G[a]["length"])
        fwd = []
        for i in iins.get(a, []):
            l, r_ = int(i["left"]), int(i["right"])
            fwd.append((l, r_) if i["strand"] == "1" else (glen - r_, glen - l))
        for st, en, gene in ints:
            n_int += 1
            if any(min(en, y) - max(st, x) >= 0.5 * (en - st) for x, y in fwd):
                n_inf += 1
            for edge in (st, en):
                n_edge += 1
                hit = [e for e in iev[a] if e["class"] == "clean" and abs(int(e["pos_fwd"]) - edge) <= 90]
                n_edge_fs += bool(hit)
                edge_fs_split += any(e["via"] == "split" for e in hit)
    pl = truth["plants"]
    outc = Counter()
    for p in pl:
        m = match(iev[p["copy"]], p["copy"], p["pos"], p["shift"])
        outc["clean" if (m is not None and m["class"] == "clean") else
             (m["class"] if m is not None else ("aligned_no_event" if covered(al.get(p["copy"], []), p["pos"]) else "not_aligned"))] += 1
    cov_near = 1 - outc["not_aligned"] / max(1, len(pl))
    CDS = cds_parts()
    gcov = []
    for a in truth["introns"]:
        tot = ov = 0
        iv = [(x, y) for x, y, _ in al.get(a, [])]
        for gene, _, parts in CDS[a]:
            if gene not in FL.CORE:
                continue
            for x, y in parts:
                tot += y - x
                ov += min(y - x, sum(max(0, min(y, b) - max(x, a_)) for a_, b in iv))
        if tot:
            gcov.append(ov / tot)
    gcov.sort()
    gsum = list(csv.DictReader(open(OUT / "introns/scan/genomes.tsv"), delimiter="\t"))
    ms = sorted(float(g["masked_share"]) for g in gsum if "__" not in g["accession"])
    r = {"B9a": n_inf / max(1, n_int), "n_introns": n_int, "B9b": n_edge_fs / max(1, n_edge), "n_edges": n_edge,
         "B9b_split": edge_fs_split, "B9c": outc["clean"] / max(1, len(pl)), "n_plants": len(pl),
         "near_edge_covered": cov_near, "near_edge_outcomes": dict(outc.most_common()),
         "core_cov_median": gcov[len(gcov) // 2], "core_cov_p10": gcov[len(gcov) // 10], "n_genomes": len(gcov),
         "masked_median": ms[len(ms) // 2], "masked_p90": ms[int(len(ms) * 0.9)], "masked_max": ms[-1]}
    if verbose:
        print(json.dumps(r, indent=1))
    return r


def pct(x: float) -> str:
    return f"{x:.3f}"


def site_clusters(events, tol_rel=None, tol_res=None) -> int:
    """Number of distinct pool sites among clean events (same gene and direction; old rule: relative subject
    position within tol_rel; new rule: subject positions within tol_res residues)."""
    sites = []
    for e in sorted(events, key=lambda e: (e["gene"], int(e["shift"]), int(e["spos"]))):
        key = (e["gene"], int(e["shift"]))
        r = float(e["spos"]) / float(e["slen"])
        for s_ in sites:
            if s_[0] == key and ((tol_rel is not None and abs(s_[1] - r) <= tol_rel) or
                                 (tol_res is not None and abs(s_[2] - int(e["spos"])) <= tol_res)):
                break
        else:
            sites.append((key, r, int(e["spos"])))
    return len(sites)


def stage_bars():
    G = meta()
    path = OUT / "plant_null/null.json"
    null_s, null_h = FL.Null(path, split_votes=True), FL.Null(path, split_votes=False)
    top, gates, bars, limits, notes, extra = [], [], [], [], [], []
    gate_ok = bar_ok = True

    def row(lst, name, value, target, passed=None):
        res = "—" if passed is None else ("PASS" if passed else "**MISS**")
        lst.append(f"| {name} | {value} | {target} | {res} |")

    P = list(csv.DictReader(open(OUT / "plant_null/plants.tsv"), delimiter="\t"))
    trials = json.load(open(OUT / "sites/trials.json"))
    sgc = json.load(open(OUT / "sites/gc.json"))
    sev_all = defaultdict(list)
    for e in load_events(OUT / "sites/scan"):
        sev_all[e["accession"]].append(e)
    sal = json.load(open(OUT / "sites/scan/aligned_fwd.json"))
    nev_all = load_events(OUT / "null_genomes/scan")
    ngc = json.load(open(OUT / "null_genomes/gc.json"))
    codes = {k: int(G[k.split("__")[0]]["code"]) for k in ngc}

    # ---- split-hit probation (PREREG log 07:10)
    split_only_plants = sum(p["clean_split_only"] == "1" for p in P)
    split_only_sites = 0
    for t in trials:
        evs = sev_all[t["copy"]]
        for s_ in t["sites"]:
            ms = match([e for e in evs if e["via"] == "split" and e["class"] == "clean"], t["copy"], s_["pos"], 1)
            mh = match([e for e in evs if e["via"] != "split" and e["class"] == "clean"], t["copy"], s_["pos"], 1)
            split_only_sites += ms is not None and mh is None
    fl_s = flags_for(nev_all, ngc, codes, null_s)
    fl_h = flags_for(nev_all, ngc, codes, null_h)
    added = sorted(k for k in ngc if any(r["flag"] for r in fl_s[k]) and not any(r["flag"] for r in fl_h[k]))
    n_split_events = {"plant_null": sum(1 for _ in open(OUT / "plant_null/scan/events.tsv")) and
                      sum(e["via"] == "split" for e in load_events(OUT / "plant_null/scan")),
                      "null_genomes": sum(e["via"] == "split" for e in nev_all),
                      "sites": sum(e["via"] == "split" for v in sev_all.values() for e in v)}
    split_votes = (split_only_plants + split_only_sites) >= 1 and not added
    null = null_s if split_votes else null_h
    cap = null.cap
    sev = {k: eligible(v, split_votes) for k, v in sev_all.items()}
    nev = eligible(nev_all, split_votes)
    probation = (f"split events: {n_split_events}; planted shifts recovered clean only via a split event: "
                 f"{split_only_plants} (planted null) + {split_only_sites} (same-codon sites); null sequences that flag only "
                 f"with split events voting: {len(added)} {added[:5]}. Decision: split-hit shifts "
                 + ("VOTE." if split_votes else "do NOT vote (exploratory tally only)."))

    # ---- planted single shifts (planted null), under the chosen regime
    sfx = "" if split_votes else "_h"
    elig = [p for p in P if p["aligned"] == "1" and p["hp_context"] == "0"]
    rec = sum(p["clean" + sfx] == "1" for p in elig) / max(1, len(elig))
    found = sum(p["found" + sfx] == "1" for p in elig) / max(1, len(elig))

    def p_votes(p):
        return FL.votes({"window": p["window" + sfx], "codon": p["codon_R" + sfx], "codon_left": p["codon_L" + sfx]}, cap)

    def p_outcome(p):
        if p["found" + sfx] != "1":
            return "aligned_no_event" if p["covered"] == "1" else "not_aligned"
        if p["class" + sfx] != "clean":
            return p["class" + sfx]
        v = p_votes(p)
        if not v:
            return "clean_wide_window"
        return "votes" if v & {p["true_R"], p["true_L"]} else "clean_other_codon"
    single = defaultdict(Counter)
    for p in P:
        single[p["group"]][p_outcome(p)] += 1
        single["all"][p_outcome(p)] += 1
    split4 = defaultdict(lambda: [0, 0, 0, 0, 0])  # planted: n, in gene, at site, clean at site, right codon
    for p in P:
        for gname in (p["group"], "all"):
            x = split4[("random", gname)]
            x[0] += 1
            x[1] += p["in_gene" + sfx] == "1"
            x[2] += p["found" + sfx] == "1"
            x[3] += p["clean" + sfx] == "1"
            x[4] += p_outcome(p) == "votes"

    # ---- same-codon trials
    def trial_eval(nl, a_sense, a_stop):
        out = []
        for t in trials:
            evs = sev[t["copy"]]
            fl = FL.genome_tests(evs, t["code"], sgc[t["copy"]], nl,
                                 {c for c, x in F.codon_table(t["code"]).items() if x == "*"}, a_sense, a_stop)
            flagged = any(r["flag"] and r["codon"] == t["codon"] and r["shift"] == 1 for r in fl)
            outs, in_gene = [], []
            for s_ in t["sites"]:
                outs.append(site_outcome(match(evs, t["copy"], s_["pos"], 1), sal.get(t["copy"], []), s_["pos"], t["codon"], nl.cap))
                in_gene.append(any(int(e["shift"]) == 1 and e["gene"] == s_["gene"] for e in evs))
            out.append({**t, "flagged": flagged, "outcomes": outs, "in_gene": in_gene, "all_vote": all(o == "votes" for o in outs)})
        return out
    TE = trial_eval(null, FL.ALPHA_SENSE, FL.ALPHA_STOP)
    for t in TE:
        for o, ig in zip(t["outcomes"], t["in_gene"]):
            for gname in (t["group"], "all"):
                x = split4[("codon", gname)]
                x[0] += 1
                x[1] += ig
                x[2] += o not in ("not_aligned", "aligned_no_event")
                x[3] += o in ("votes", "clean_other_codon", "clean_wide_window")
                x[4] += o == "votes"

    def rate(kind, n, cond=False, rows=None):
        T = [t for t in (rows or TE) if t["kind"] == kind and t["n_sites"] == n and (t["all_vote"] or not cond)]
        return sum(t["flagged"] for t in T) / max(1, len(T)), len(T)

    top.append("## Headline sensitivity (unconditional)")
    top.append("")
    top.append(f"Window rule, window cap = {cap} nt (smallest width holding >= {WINDOW_SHARE:.0%} of clean planted shifts); "
               f"alpha_sense = {FL.ALPHA_SENSE:g}, alpha_stop = {FL.ALPHA_STOP:g}; split-hit shifts "
               f"{'vote' if split_votes else 'do not vote'}.")
    top.append("")
    top.append("| planted sites at one codon | sense: genome flagged | stop: genome flagged | sense, all sites counted | stop, all sites counted |")
    top.append("|---|---|---|---|---|")
    for n in (2, 3, 5):
        (rs, ns), (rt, nt) = rate("sense", n), rate("stop", n)
        (rsc, nsc), (rtc, ntc) = rate("sense", n, True), rate("stop", n, True)
        top.append(f"| {n} | {pct(rs)} (n={ns}) | {pct(rt)} (n={nt}) | {pct(rsc)} (n={nsc}) | {pct(rtc)} (n={ntc}) |")
    top.append("")
    top.append("\"All sites counted\" = every planted site detected as a clean shift voting for the planted codon (the "
               "conditional B5 definition). The unconditional columns are the headline.")
    top.append("")
    top.append("### Planted shifts: detected in the gene, at the site, and at the right codon (share of planted)")
    top.append("")
    top.append("| planted | group | n | in gene (same direction, anywhere) | at site (±9 nt) | clean at site | counted at right codon |")
    top.append("|---|---|---|---|---|---|---|")
    for kind, label in (("random", "random positions (planted null)"), ("codon", "same-codon trial sites")):
        for gname in ("animals", "fungi", "protists", "all"):
            x = split4[(kind, gname)]
            if x[0]:
                top.append(f"| {label} | {gname} | {x[0]} | " + " | ".join(pct(v / x[0]) for v in x[1:]) + " |")
    top.append("")
    site_out = defaultdict(Counter)
    for t in TE:
        for o in t["outcomes"]:
            site_out[t["group"]][o] += 1
            site_out["all"][o] += 1
    order = ["votes", "clean_other_codon", "clean_wide_window", "not_aligned", "aligned_no_event", "masked_insert_edge",
             "homopolymer", "gene_end", "gene_start", "flank_fail", "paired"]
    for title, tab in (("same-codon trial sites", site_out), ("random planted shifts (planted null)", single)):
        top.append(f"### Reasons, {title} (share of sites)")
        top.append("")
        top.append("| group | sites | counted (clean, votes for its codon) | " + " | ".join(order[1:]) + " |")
        top.append("|---|---|---|" + "---|" * (len(order) - 1))
        for gname in ("animals", "fungi", "protists", "all"):
            c = tab[gname]
            tot = sum(c.values())
            if tot:
                top.append(f"| {gname} | {tot} | " + " | ".join(pct(c[o] / tot) for o in order) + " |")
        top.append("")

    # ---- gates (scan auto-start conditions, amendment 4)
    fl = flags_for(nev, ngc, codes, null)
    bad = {k: [r for r in v if r["flag"]] for k, v in fl.items()}
    bad = {k: v for k, v in bad.items() if v}
    by_kind = Counter(("original" if "__" not in k else k.split("__")[1]) for k in bad)
    ok = len(bad) == 0
    gate_ok &= ok
    row(gates, "B3/B7 null genomes (120) + homopolymer and random-indel copies: genome flags",
        f"{len(bad)} of {len(ngc)} {dict(by_kind)} " + "; ".join(f"{k}:{r['shift']:+d}{r['codon']} k={r['k']}/{r['n']} p={r['p']:.1e}" for k, v in list(bad.items())[:6] for r in v[:1]),
        "0", ok)
    I = intron_summary(G, verbose=False, split_votes=split_votes)
    ok = I["B9b"] <= 0.02
    gate_ok &= ok
    row(gates, "B9b intron edges with a false clean shift within 90 nt", f"{I['B9b']:.4f} (n={I['n_edges']}; via split {I['B9b_split']})", "<= 0.02", ok)
    named_ev = eligible(load_events(OUT / "named/scan"), split_votes)
    named_gc = json.load(open(OUT / "named/gc.json"))
    nfl = flags_for(named_ev, named_gc, {a: int(G[a]["code"]) for a in NAMED}, null)
    ins = list(csv.DictReader(open(OUT / "named/scan/inserts.tsv"), delimiter="\t")) if (OUT / "named/scan/inserts.tsv").exists() else []
    for a, name in NAMED.items():
        if a == "NC_024095.1":
            nb = sum(1 for i in ins if i["accession"] == a and i["class"] == "bypass_like")
            ok = nb > 0
            row(gates, f"B6 {name}: bypass-like inserts", nb, "> 0", ok)
        else:
            hits = [r for r in nfl[a] if r["flag"]]
            ok = bool(hits)
            row(gates, f"B6 {name}: flags", "; ".join(f"{r['shift']:+d}{r['codon']} k={r['k']}/{r['n']} genes={r['genes']} p={r['p']:.1e}" for r in hits) or "no flag", "must flag", ok)
        gate_ok &= ok

    # ---- bars (a miss is logged as a limit)
    ok = rec >= 0.90
    bar_ok &= ok
    row(bars, "B1 planted single shifts, clean recall in aligned genes (homopolymers excluded)", f"{rec:.3f} (n={len(elig)}; detected {found:.3f})", ">= 0.90", ok)
    rsc, nsc = rate("sense", 2, True)
    ok = rsc >= 0.90
    bar_ok &= ok
    row(bars, "B5 two sites at one rare sense codon, flagged, when both counted (conditional)", f"{rsc:.3f} (n={nsc})", ">= 0.90", ok)
    J = [j for j in csv.DictReader(open(ROOT / "results/audit/known_junctions_refseq238.tsv"), delimiter="\t") if j["group"] != "land_plants"]
    aev = defaultdict(list)
    for e in eligible(load_events(OUT / "annotated/scan"), split_votes):
        aev[e["accession"]].append(e)
    scanned = {g["accession"] for g in csv.DictReader(open(OUT / "annotated/scan/genomes.tsv"), delimiter="\t")}
    J = [j for j in J if j["accession"] in scanned]
    aseq = seqs_for(sorted({j["accession"] for j in J}), G)
    CDS = cds_parts()
    res = defaultdict(lambda: [0, 0, 0, 0, 0, 0])  # n, at site (bar), clean at site, via split, in gene, right codon
    offsets = []
    for j in J:
        a = j["accession"]
        d = int(j["shift"])
        fold = 1 if d % 3 == 1 else -1
        grp = "nad3" if re.search(r"ND3|NAD3|SUBUNIT 3$", j["gene"].upper()) and j["group"] == "Metazoa" else "other"
        pos = int(j["genome_pos"]) - 1
        m = match(aev[a], a, pos, fold)
        cds = next(((g_, st_, pa) for g_, st_, pa in CDS[a] if min(x for x, _ in pa) <= pos < max(y for _, y in pa)), None)
        same_dir = [e for e in aev[a] if int(e["shift"]) == fold and cds and e["gene"] == cds[0]]
        in_gene = bool(same_dir)
        right = False
        if m is not None and m["class"] == "clean" and cds:
            seq = aseq[a]
            strand = cds[1]
            so = seq if strand == 1 else seq.translate(S.COMP)[::-1]
            opos = (pos if strand == 1 else len(seq) - 1 - pos) if fold == 1 else (pos + 1 if strand == 1 else len(seq) - pos)
            f0 = F.old_frame0(cds[2], strand, len(seq), opos)
            tr = F.trigger(F.Event(a, "", strand, opos, fold, fold, 0, 0), so, F.codon_table(int(G[a]["code"])), frame0=f0)
            right = bool(FL.votes(m, cap) & {tr["codon"], tr["codon_left"]})
        if in_gene and m is None:
            offsets.append(min(abs(int(e["pos_fwd"]) - pos) for e in same_dir))
        x = res[grp]
        x[0] += 1
        x[1] += m is not None
        x[2] += m is not None and m["class"] == "clean"
        x[3] += m is not None and m["via"] == "split"
        x[4] += in_gene
        x[5] += right
    for grp in ("nad3", "other"):
        n, f_, c, sp, ig, rc = res[grp]
        ok = f_ / max(1, n) >= 0.80
        bar_ok &= ok
        row(bars, f"B2 RefSeq-annotated shifts detected at the site, ±9 nt ({grp})",
            f"{f_/max(1,n):.3f} (n={n}; clean {c/max(1,n):.3f}; via split {sp}) | in gene anywhere {ig/max(1,n):.3f} | "
            f"in the voting window, right codon {rc/max(1,n):.3f}", ">= 0.80 (at site)", ok)
    offsets.sort()

    # ---- pool same-site tolerance check (PREREG log 07:10)
    if offsets:
        within = [o for o in offsets if o <= 150]
        extra.append(f"- Placement offsets of annotated shifts detected in the gene but not at the site (n={len(offsets)}): "
                     f"median {offsets[len(offsets)//2]} nt, 95th pct {offsets[int(len(offsets)*0.95)]} nt, max {offsets[-1]} nt; "
                     f"{len(within)} of {len(offsets)} within 150 nt (= the 50-residue pool tolerance). Offsets beyond 150 nt are "
                     "usually a different shift elsewhere in the gene.")
    birds = defaultdict(list)
    for a, evs in aev.items():
        g_ = G[a]
        if g_["class_name"] == "Aves" or g_["order_pool_name"] == "Testudines":
            for e in evs:
                if e["class"] == "clean" and e["gene"] == "nad3":
                    birds[("family", g_["family_pool_name"])].append(e)
                    birds[("order", g_["order_pool_name"])].append(e)
    old_n = sorted((site_clusters(v, tol_rel=0.02), k) for k, v in birds.items())
    new_n = sorted((site_clusters(v, tol_res=FL.POOL_SITE_RESIDUES), k) for k, v in birds.items())
    if old_n:
        extra.append(f"- Bird and turtle nad3 clean shifts in the annotated set, distinct pool sites per pool ({len(old_n)} family/order pools): "
                     f"old 2% rule max {old_n[-1][0]} ({old_n[-1][1][1]}), pools with >= 3 sites {sum(n_ >= 3 for n_, _ in old_n)}; "
                     f"new 50-residue rule max {new_n[-1][0]} ({new_n[-1][1][1]}), pools with >= 3 sites {sum(n_ >= 3 for n_, _ in new_n)}.")

    # ---- >= 3 genes tally baseline in the null genomes (exploratory)
    by_seq = defaultdict(set)
    for e in FL.clean_core(nev, split_votes):
        by_seq[e["accession"]].add(e["gene"])
    n3 = Counter(("original" if "__" not in k else k.split("__")[1]) for k, gset in by_seq.items() if len(gset) >= 3)
    extra.append(f"- Exploratory baseline: null sequences with clean shifts in >= 3 different core genes: {dict(n3)} (of 120 each).")

    # ---- stated limits (amendment 4, item 5)
    row(limits, "B9a annotated core-gene introns inferred as inserts", f"{I['B9a']:.3f} (n={I['n_introns']})", "limit (was >= 0.80)")
    row(limits, "B9c shifts planted 31-90 nt from intron edges: clean recall", f"{I['B9c']:.3f} (n={I['n_plants']})", "limit (was >= 0.80)")
    row(limits, "near-edge planted shifts inside an aligned hit", f"{I['near_edge_covered']:.3f}", "limit")
    row(limits, "near-edge outcomes", ", ".join(f"{k} {v}" for k, v in I["near_edge_outcomes"].items()), "limit")
    row(limits, "fungal intron genomes: per-genome core-CDS coverage", f"median {I['core_cov_median']:.3f}, 10th pct {I['core_cov_p10']:.3f} (n={I['n_genomes']})", "limit")
    row(limits, "fungal intron genomes: masked coding share", f"median {I['masked_median']:.3f}, 90th pct {I['masked_p90']:.3f}, max {I['masked_max']:.3f}", "limit")

    # ---- what the pre-written alpha rule would give on this run (record only; alpha stays fixed)
    def av(a):
        return any(r["flag"] for r in flags_for(named_ev, named_gc, {"NC_010769.1": 4}, null, a, a)["NC_010769.1"])
    rule_sense = next((a for a in GRID if av(a) and rate("sense", 2, True, trial_eval(null, a, a))[0] >= 0.90), None)
    rule_stop = next((a for a in GRID if rate("stop", 2, True, trial_eval(null, 1e-12, a))[0] >= 0.90), None)
    notes.append(f"- Split-hit probation: {probation}")
    notes.append(f"- Pre-written alpha rule applied to this run (record only; alpha stays {FL.ALPHA_SENSE:g} / {FL.ALPHA_STOP:g}): "
                 f"alpha_sense would be {rule_sense}, alpha_stop {rule_stop}.")

    # ---- expected chance flags in the scan, estimated from the 120 null genomes (originals)
    scan_set = [r for r in G.values() if r["set"] in ("refseq_main", "fungal_tpa") and int(r["code"] or 0) > 0 and int(r["length"]) > 0]
    by_ph = defaultdict(list)
    nmax = Counter()
    for k in ngc:
        if "__" in k:
            continue
        stops = {c for c, x in F.codon_table(codes[k]).items() if x == "*"}
        inputs = FL.genome_dir_inputs([e for e in nev if e["accession"] == k], codes[k], ngc[k], null, stops)
        for n_, _, _ in inputs:
            nmax[n_] += 1
        e_k = sum(FL.expected_chance(n_, qs, al_) for n_, qs, al_ in inputs)
        by_ph[(group_of(G[k]), G[k]["phylum_name"] or G[k]["kingdom_name"] or "unranked")].append(e_k)
    grp_mean = {gname: (lambda v: sum(v) / len(v) if v else 0.0)([x for (gg, _), xs in by_ph.items() if gg == gname for x in xs])
                for gname in ("animals", "fungi", "protists")}
    est = 0.0
    for r in scan_set:
        key = (group_of(r), r["phylum_name"] or r["kingdom_name"] or "unranked")
        v = by_ph.get(key)
        est += (sum(v) / len(v)) if v else grp_mean[group_of(r)]
    notes.append(f"- Expected chance flags across the scan ({len(scan_set)} genomes x 2 directions), estimated from the 120 null "
                 f"genomes' clean-shift counts weighted by phylum (group mean where a phylum is not sampled): {est:.2f}. "
                 f"Clean core shifts per genome-direction in the null originals: {dict(sorted(nmax.items()))}. Upper bound "
                 "(ignores the >= 2 genes condition). Recomputed exactly from the scan's own counts after the scan.")
    json.dump({"window_cap": cap, "alpha_sense": FL.ALPHA_SENSE, "alpha_stop": FL.ALPHA_STOP, "split_votes": split_votes,
               "probation": probation, "rule_sense": rule_sense, "rule_stop": rule_stop,
               "expected_chance_flags_est": est, "gate_ok": gate_ok, "bars_ok": bar_ok, "intron": I,
               "pool_site_residues": FL.POOL_SITE_RESIDUES}, open(OUT / "thresholds.json", "w"), indent=1)

    hdr = ["# Control bars (re-run after DESIGN.md amendment 4 and the PREREG log entry of 07:10)", ""]
    body = top + ["## Scan gates (all must pass for the scan to start)", "", "| Gate | Value | Target | Result |", "|---|---|---|---|"] + gates
    body += ["", "## Bars (a miss is logged in PREREG.md as a limit)", "", "| Bar | Value | Target | Result |", "|---|---|---|---|"] + bars
    body += ["", "## Stated limits (fungal introns)", "", "| Item | Value | Status | |", "|---|---|---|---|"] + limits
    body += ["", "## Misplacement and pool sites", ""] + extra
    body += ["", "## Notes", ""] + notes
    text = "\n".join(hdr + body) + "\n"
    (OUT / "BARS.md").write_text(text)
    print(text)
    print(f"gates {'PASS' if gate_ok else 'FAIL'}; bars {'all pass' if bar_ok else 'some missed (limits)'}")
    sys.exit(2 if not gate_ok else (1 if not bar_ok else 0))


STAGES = {"plant_null": stage_plant_null, "null_genomes": stage_null_genomes, "sites": stage_sites,
          "annotated": stage_annotated, "named": stage_named, "introns": stage_introns,
          "intron_summary": intron_summary, "bars": stage_bars}

if __name__ == "__main__":
    todo = sys.argv[1:] or ["all"]
    if todo == ["all"]:
        todo = ["plant_null", "null_genomes", "named", "annotated", "introns", "sites", "bars"]
    for st in todo:
        STAGES[st]()
