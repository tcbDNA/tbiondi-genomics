"""Read support at the four G. jiainica sites with PCR duplicates removed.

Same exact-match counts as step f (`grigorovia.py step_f`) and `grigorovia_h1h2.py reads_strand`: 20 nt + extra base +
20 nt (or the 40 nt without it), either strand; and for the two close cob sites, the window 26462-26511 with both, one
or neither base. Counted twice: over all reads, and over reads not flagged as duplicates by samtools markdup in
`data/work/trna_gly_tcc/dedup.bam` (made by `trna_gly_tcc.py`; a read pair is a duplicate when another pair has the same
mapped start positions). The raw counts must equal the published ones (checked below).

Writes results/explore/grigorovia/step_f_reads_dedup.md and adds the duplicate-removed numbers to step_f_reads.md,
sites_for_lab.md and h1h2/reads_strand.md.

    uv run python src/explore/site_reads_dedup.py      (run trna_gly_tcc.py first)
"""
import gzip
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / "results/explore/grigorovia"
BAM = ROOT / "data/work/trna_gly_tcc/dedup.bam"
SITES = {"cox1@24710": 24710, "cob@26495": 26495, "cob@26476": 26476, "cob@26265": 26265}   # 0-based forward
PUBLISHED = {"cox1@24710": 1440, "cob@26495": 1718, "cob@26476": 1733, "cob@26265": 1785, "both": 1516}
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(x):
    return x.translate(COMP)[::-1]


def main():
    seq = "".join(l.strip() for l in open(G / "BK073075.1.fna") if not l.startswith(">")).upper()
    pats = {}
    for k, p in SITES.items():
        w, wo = seq[p - 20:p + 21], seq[p - 20:p] + seq[p + 1:p + 21]
        pats[k] = {"with": (w, rc(w)), "without": (wo, rc(wo))}
    a, b = 26476, 26495
    joint = seq[a - 15:b + 16]
    ia, ib = 15, 15 + (b - a)
    jv = {"both bases": joint, "lacks 26477 only": joint[:ia] + joint[ia + 1:],
          "lacks 26496 only": joint[:ib] + joint[ib + 1:], "lacks both": joint[:ia] + joint[ia + 1:ib] + joint[ib + 1:]}
    jv = {k: (v, rc(v)) for k, v in jv.items()}
    dups = set(subprocess.run(f"samtools view -f 1024 {BAM} | cut -f1", shell=True, check=True, capture_output=True,
                              text=True).stdout.split())
    cnt = {t: {k: {"with": 0, "without": 0} for k in SITES} for t in ("all", "dedup")}
    jc = {t: {k: 0 for k in jv} for t in ("all", "dedup")}
    for fq in sorted((ROOT / "data/reads/SRR16974510").glob("*.fastq.gz")):
        with gzip.open(fq, "rt") as fh:
            name = None
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    name = line[1:].split()[0]
                    continue
                if i % 4 != 1:
                    continue
                tags = ("all",) if name in dups else ("all", "dedup")
                for k, d in pats.items():
                    for kind, (x, y) in d.items():
                        if x in line or y in line:
                            for t in tags:
                                cnt[t][k][kind] += 1
                for k, (x, y) in jv.items():
                    if x in line or y in line:
                        for t in tags:
                            jc[t][k] += 1
    for k in SITES:
        assert cnt["all"][k]["with"] == PUBLISHED[k], (k, cnt["all"][k]["with"])
    assert jc["all"]["both bases"] == PUBLISHED["both"], jc["all"]["both bases"]
    L = ["# Read support at the four sites with PCR duplicates removed (SRR16974510; exact 20 + base + 20 nt)", "",
         f"Duplicates: read pairs flagged by samtools markdup after competitive mapping to BK073075.1 plus the nuclear "
         f"assembly ({len(dups)} read names flagged, each a pair or a single read). Raw counts equal the published step f counts.", "",
         "| site | with the base | without | with, duplicates removed | without, duplicates removed |",
         "|---|---|---|---|---|"]
    for k in SITES:
        c, d = cnt["all"][k], cnt["dedup"][k]
        L.append(f"| {k} | {c['with']} | {c['without']} | {d['with']} | {d['without']} |")
    L += ["", "Reads spanning both cob sites (26477 and 26496; window 26462-26511):", "",
          "| version | reads | duplicates removed |", "|---|---|---|"]
    for k in jv:
        L.append(f"| {k} | {jc['all'][k]} | {jc['dedup'][k]} |")
    (G / "step_f_reads_dedup.md").write_text("\n".join(L) + "\n")
    # add the duplicate-removed numbers to the tables that quote the raw counts
    sf = G / "step_f_reads.md"
    s = sf.read_text()
    if "duplicates removed" not in s:
        s = s.replace("| shift | reads with the extra base | reads without it | share with |\n|---|---|---|---|",
                      "| shift | reads with the extra base | reads without it | share with | with / without, duplicates removed |\n|---|---|---|---|---|")
        for k in SITES:
            d = cnt["dedup"][k]
            s = re.sub(rf"(\| {re.escape(k)} \| \d+ \| \d+ \| [0-9.]+ \|)\n", rf"\1 {d['with']} / {d['without']} |\n", s)
        sf.write_text(s)
    sl = G / "sites_for_lab.md"
    s = sl.read_text()
    if "duplicates removed" not in s:
        s = s.replace("| reads with / without |", "| reads with / without (duplicates removed) |")
        for k, p in SITES.items():
            c, d = cnt["all"][k], cnt["dedup"][k]
            s = s.replace(f"| {c['with']} / {c['without']} |", f"| {c['with']} / {c['without']} ({d['with']} / {d['without']}) |")
        sl.write_text(s)
    rs = G / "h1h2/reads_strand.md"
    s = rs.read_text()
    if "duplicates removed" not in s:
        s = s.replace("| version | reads |\n|---|---|", "| version | reads | duplicates removed |\n|---|---|---|")
        for k in jv:
            s = re.sub(rf"(\| {re.escape(k)} \| \d+ \|)\n", rf"\1 {jc['dedup'][k]} |\n", s)
        rs.write_text(s)
    print("\n".join(L))


if __name__ == "__main__":
    main()
