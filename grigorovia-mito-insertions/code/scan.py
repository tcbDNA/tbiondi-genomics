"""Run the frameshift scan on a set of genomes and write per-shift and per-genome tables.

    uv run python code/scan.py --genomes <list.txt|all-main|all-fungal|...> --out <dir> [--threads 12]

`process_genome` is the single per-genome function used by the controls, the planted null and the scan.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fscore as F  # noqa: E402

ROOT = F.ROOT
WORK = ROOT / "data/work"
COMP = str.maketrans("ACGTN", "TGCAN")


def load_meta() -> tuple[dict, dict, dict]:
    genomes = {r["accession"]: r for r in csv.DictReader(open(WORK / "genomes.tsv"), delimiter="\t")}
    ref = {r["ref_id"]: r for r in csv.DictReader(open(WORK / "ref_proteins.tsv"), delimiter="\t")}
    ref_family = {k: v["family_pool"] for k, v in ref.items()}
    ref_gene = {k: v["gene"] for k, v in ref.items()}
    return genomes, ref_family, ref_gene


def oriented(seq: str, strand: int) -> str:
    return seq if strand == 1 else seq.translate(COMP)[::-1]


def process_genome(acc: str, seq: str, code: int, family: str, hits: list[F.Hit], ref_family: dict,
                   ref_gene: dict) -> dict:
    """All per-genome logic: cull, shift events, inserts, masks, filters, trigger codons, tallies."""
    table = F.codon_table(code)
    kept = F.cull(hits, family, ref_family)
    seqs = {1: oriented(seq, 1), -1: oriented(seq, -1)}
    ins = F.inserts(kept, ref_gene)
    events, stops = [], []
    aligned = {1: [], -1: []}
    cols_of = {}
    for i, h in enumerate(kept):
        aligned[h.strand].append((h.ostart, h.oend))
        evs, st, qins, cols_of[i] = F.walk(h, seqs[h.strand])
        gene = ref_gene.get(h.s, "?")
        for a, b in qins:  # insertions inside an HSP (e.g. bypass elements aligned through as gaps)
            n = b - a
            if n >= F.P["bypass_min"]:
                ins.append({"strand": h.strand, "gene": gene, "left": a, "right": b, "len": n, "case": 0,
                            "class": "intron" if n > F.P["bypass_max"] else "bypass_like", "within_hsp": 1})
        for pos, aa in st:
            stops.append({"accession": acc, "gene": gene, "strand": h.strand, "pos": pos, "subject_aa": aa})
        for e in evs:
            events.append((h, gene, e))
    # frame changes between two consecutive hits of the same reference protein (amendment 4, item 4)
    for i, e in F.split_events(kept, cols_of):
        events.append((kept[i], ref_gene.get(e.s, "?"), e))
    edges = defaultdict(list)
    for i in ins:
        i.setdefault("within_hsp", 0)
        if i["len"] > F.P["insert_min"]:
            edges[i["strand"]] += [i["left"], i["right"]]
    rows = []
    for h, gene, e in events:
        so = seqs[e.strand]
        trig = F.trigger(e, so, table)
        a, b = min(e.pos, trig["tpos"]) - 1, max(e.pos, trig["tpos"]) + 1
        hp, hp_base = F.homopolymer(so, a, b, e.shift)
        # single-hit shifts are classified without reference to split events (PREREG log 07:10)
        partners = [o for _, _, o in events if o is not e and o.strand == e.strand and o.shift == -e.shift and o.s == e.s
                    and (e.via == "split" or o.via == "hsp")]
        paired = any(abs(o.pos - e.pos) <= F.P["pair_window"] for o in partners)
        masked = any(abs(e.pos - x) <= F.P["mask_nt"] for x in edges[e.strand])
        lid = sum(e.left_cols) / max(1, len(e.left_cols))
        rid = sum(e.right_cols) / max(1, len(e.right_cols))
        flank_ok = (len(e.left_cols) >= F.P["flank_cols"] and len(e.right_cols) >= F.P["flank_cols"]
                    and lid >= F.P["flank_identity"] and rid >= F.P["flank_identity"])
        gene_start = e.spos <= F.P["gene_start_aa"]
        gene_end = e.spos > e.slen - F.P["gene_end_aa"]
        cls = ("paired" if paired else "masked_insert_edge" if masked else "gene_start" if gene_start else
               "gene_end" if gene_end else "flank_fail" if not flank_ok else "homopolymer" if hp else "clean")
        detour_stop = False
        if paired:
            other = min(partners, key=lambda o: abs(o.pos - e.pos))
            first = e if e.pos <= other.pos else other
            lo, hi = first.pos, max(e.pos, other.pos)
            # old frame (the frame before the first shift) read straight through the detour
            old = so[lo:hi + 3]
            detour_stop = "*" in "".join(table.get(old[i:i + 3], "X") for i in range(0, len(old) - 2, 3))
        fwd_pos = e.pos if e.strand == 1 else len(seq) - 1 - e.pos
        rows.append({
            "accession": acc, "gene": gene, "subject": e.s, "strand": e.strand, "pos_oriented": e.pos,
            "pos_fwd": fwd_pos, "shift": e.shift, "raw": e.raw, "spos": e.spos, "slen": e.slen,
            "codon": trig["codon"], "codon_left": trig["codon_left"], "window": trig["window"],
            "is_stop_trigger": int(table.get(trig["codon"], "X") == "*"),
            "class": cls, "homopolymer_base": hp_base, "left_id": round(lid, 2), "right_id": round(rid, 2),
            "detour_around_stop": int(detour_stop), "bits": h.bits, "via": e.via,
        })
    # coverage: aligned oriented bases (union) and masked share
    cov = 0
    masked_nt = 0
    for s in (1, -1):
        iv = sorted(aligned[s])
        merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        cov += sum(b - a for a, b in merged)
        for x in edges[s]:
            lo, hi = x - F.P["mask_nt"], x + F.P["mask_nt"]
            masked_nt += sum(max(0, min(b, hi) - max(a, lo)) for a, b in merged)
    genes_hit = sorted({ref_gene.get(h.s, "?") for h in kept})
    aligned_fwd = [((h.qstart - 1, h.qend) if h.qframe > 0 else (h.qend - 1, h.qstart)) + (ref_gene.get(h.s, "?"),)
                   for h in kept]
    acgt = sum(seq.count(b) for b in "ACGT")
    summary = {"accession": acc, "gc": round((seq.count("G") + seq.count("C")) / max(1, acgt), 4), "n_hits": len(kept), "aligned_nt": cov, "masked_nt": masked_nt,
               "masked_share": round(masked_nt / cov, 4) if cov else 0.0, "genes_hit": ",".join(genes_hit),
               "n_genes_hit": len(genes_hit), "n_events": len(rows),
               "n_clean": sum(r["class"] == "clean" for r in rows),
               "n_split": sum(r["via"] == "split" for r in rows),
               "n_inserts_intron": sum(i["class"] == "intron" for i in ins),
               "n_inserts_bypass_like": sum(i["class"] == "bypass_like" for i in ins),
               "n_stop_columns": len(stops)}
    return {"events": rows, "summary": summary, "inserts": [{"accession": acc, **i} for i in ins],
            "stops": stops, "aligned_fwd": aligned_fwd}


def write_tables(results: list[dict], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for key, name in (("events", "events.tsv"), ("inserts", "inserts.tsv"), ("stops", "stop_columns.tsv")):
        rows = [r for res in results for r in res[key]]
        if not rows:
            continue
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    with open(out / "genomes.tsv", "w", newline="") as fh:
        rows = [res["summary"] for res in results]
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def read_fasta(path: Path, wanted: set[str] | None = None) -> dict[str, str]:
    out, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name and (wanted is None or name in wanted):
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name and (wanted is None or name in wanted):
        out[name] = "".join(buf)
    return out


def scan(accs: list[str], out: Path, threads: int = 12, seqs: dict[str, str] | None = None,
         code_of: dict[str, int] | None = None, family_of: dict[str, str] | None = None) -> list[dict]:
    """Scan genomes (by accession; or custom sequences with given code/family) and return results."""
    genomes, ref_family, ref_gene = load_meta()
    out.mkdir(parents=True, exist_ok=True)
    by_code = defaultdict(list)
    for a in accs:
        c = code_of[a] if code_of else int(genomes[a]["code"])
        by_code[c].append(a)
    results = []
    for code, group in sorted(by_code.items()):
        t0 = time.time()
        if seqs is None:
            s = read_fasta(WORK / f"queries/code_{code}.fna", set(group))
        else:
            s = {a: seqs[a] for a in group}
        fa = out / f"q_code{code}.fna"
        with open(fa, "w") as fh:
            for a in group:
                fh.write(f">{a}\n{s[a]}\n")
        hits_path = out / f"hits_code{code}.tsv"
        if not (os.environ.get("REUSE_HITS") == "1" and hits_path.exists()):  # rescore saved hits, same DB
            F.run_diamond(fa, code, hits_path, threads)
        hits = defaultdict(list)
        for h in F.read_hits(hits_path):
            hits[h.q].append(h)
        for a in group:
            fam = family_of[a] if family_of else genomes[a]["family_pool"]
            results.append(process_genome(a, s[a], code, fam, hits.get(a, []), ref_family, ref_gene))
        fa.unlink()
        print(f"  code {code}: {len(group)} genomes, {sum(len(v) for v in s.values())/1e6:.1f} Mb, "
              f"{time.time()-t0:.0f}s", flush=True)
    return results


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True, help="file with accessions, or all-main / all-fungal / all-landplant")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=12)
    args = ap.parse_args()
    genomes, _, _ = load_meta()
    sets = {"all-main": "refseq_main", "all-fungal": "fungal_tpa", "all-landplant": "refseq_landplant"}
    if args.genomes in sets:
        accs = [a for a, r in genomes.items() if r["set"] == sets[args.genomes] and int(r["code"]) > 0
                and int(r["length"]) > 0]
    else:
        accs = [l.strip() for l in open(args.genomes) if l.strip()]
    res = scan(accs, Path(args.out), args.threads)
    write_tables(res, Path(args.out))
    import json
    with open(Path(args.out) / "aligned_fwd.json", "w") as fh:
        json.dump({r["summary"]["accession"]: r["aligned_fwd"] for r in res}, fh)


if __name__ == "__main__":
    main()
