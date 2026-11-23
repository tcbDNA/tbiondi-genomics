"""Is there a TCC (UCC) tRNA-Gly in G. jiainica's mitochondrion, and what are the reads with TCC at the anticodon?

Steps (reads SRR16974510; mitogenome BK073075.1; nuclear assembly GCA_030571715.1):
1. Map all reads competitively to the mitogenome plus the nuclear assembly (minimap2 -x sr), sort, and mark
   duplicates (samtools fixmate -m / markdup).
2. Pile up bases (base quality >= 20, mapping quality >= 20) over 2 kb either side of the tRNA-Gly gene, with and
   without duplicates. Position 16416 (+ strand) is anticodon position 34 of the minus-strand tRNA-Gly.
3. Background: the C>A fraction (the change a TCC read represents on the + strand) at every other reference-C position
   with depth >= 200, and the any-mismatch fraction at every other position.
4. For each read with A at 16416: its alignment, mismatches, mate, and whether other reads of the same molecule (its
   mate, or PCR copies with the same start positions) read C there.
5. ARAGORN 1.2.41 at a relaxed threshold (-ps) on the mitogenome; minimap2 -x asm5 of the mitogenome against the nuclear
   assembly, to look for a nuclear copy of the mitochondrial tRNA-Gly.

    uv run python code/explore/trna_gly_tcc.py    -> results/explore/grigorovia/trna_gly_tcc.md
"""
import collections as C
import re
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / "results/explore/grigorovia"
W = ROOT / "data/work/trna_gly_tcc"
W.mkdir(parents=True, exist_ok=True)
MITO = G / "BK073075.1.fna"
READS = sorted((ROOT / "data/reads/SRR16974510").glob("*.fastq.gz"))
POS, TRNA = 16416, (16378, 16449)
LO, HI = TRNA[0] - 2000, TRNA[1] + 2000


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw).stdout


def nuclear_fasta():
    out = W / "nuclear.fna"
    if not out.exists():
        z = zipfile.ZipFile(ROOT / "data/nuclear/GCA_030571715.1/GCA_030571715.1.zip")
        name = next(n for n in z.namelist() if n.endswith("_genomic.fna"))
        out.write_bytes(z.read(name))
    return out


def parse_pile(path):
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        pos, ref, s = int(f[1]), f[2].upper(), re.sub(r"\^.", "", f[4]).replace("$", "")
        names = f[6].split(",") if len(f) > 6 and f[6] else []
        bases, i = [], 0
        while i < len(s):
            if s[i] in "+-":
                m = re.match(r"[+-](\d+)", s[i:])
                i += len(m.group(0)) + int(m.group(1))
                continue
            ch = s[i]
            bases.append(ref if ch in ".," else ch.upper())
            i += 1
        out[pos] = (ref, C.Counter(b for b in bases if b in "ACGT"), list(zip(names, bases)))
    return out


def base_at(read_fields, pos):
    start, seq, b = int(read_fields[3]), read_fields[9], None
    q, r = 0, start
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", read_fields[5]):
        n = int(n)
        if op in "M=X":
            if r <= pos < r + n:
                b = seq[q + pos - r]
            q += n
            r += n
        elif op in "IS":
            q += n
        elif op in "DN":
            r += n
    return b


def main():
    nuc = nuclear_fasta()
    comb = W / "combined.fna"
    comb.write_text(MITO.read_text() + nuc.read_text())
    sh(f"samtools faidx {comb}")
    bam, dedup = W / "reads.bam", W / "dedup.bam"
    if not dedup.exists():
        sh(f"minimap2 -ax sr -t 10 {comb} {' '.join(map(str, READS))} 2>/dev/null | samtools sort -@4 -o {bam} -")
        sh(f"samtools collate -@4 -O -u {bam} {W}/col | samtools fixmate -@4 -m -u - - | samtools sort -@4 -u - | "
           f"samtools markdup -@4 -s - {dedup} 2> {W}/markdup.txt")
        sh(f"samtools index {dedup}")
    dup = {l.split(":")[0].strip(): l.split(":")[1].strip() for l in open(W / "markdup.txt") if ":" in l}
    res = {}
    for tag, flt in (("all reads", "UNMAP,SECONDARY,QCFAIL"), ("duplicates removed", "UNMAP,SECONDARY,QCFAIL,DUP")):
        pile = W / f"pile_{tag.split()[0]}.txt"
        sh(f"samtools mpileup -a -B -Q 20 -q 20 --ff {flt} --output-QNAME -f {comb} -r BK073075.1:{LO}-{HI} {dedup} "
           f"2>/dev/null > {pile}")
        res[tag] = parse_pile(pile)
    anti = {"C": "GCC", "A": "TCC", "G": "CCC", "T": "ACC"}
    L = ["# Is there a TCC tRNA-Gly in G. jiainica's mitochondrion? (post hoc, exploratory)", "",
         f"Reads SRR16974510 ({dup.get('READ', '?')} reads; {dup.get('DUPLICATE TOTAL', '?')} marked as duplicates) mapped to "
         "BK073075.1 plus the nuclear assembly GCA_030571715.1 (minimap2 -x sr). Base quality >= 20, mapping quality "
         ">= 20. Position 16416 (+ strand) is anticodon position 34 of the minus-strand tRNA-Gly (c[16378,16449]).", "",
         "| reads | depth at 16416 | GCC | TCC | CCC | ACC | TCC share | other C positions with a C>A share at least as high | median C>A share elsewhere | 95th percentile | any-mismatch share at 16416 | median any-mismatch share elsewhere |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for tag, P in res.items():
        ref, cnt, per = P[POS]
        d = sum(cnt.values())
        fr = sorted(sum(v for b, v in c.items() if b == "A") / sum(c.values())
                    for p, (r, c, _) in P.items() if r == "C" and p != POS and sum(c.values()) >= 200)
        al = sorted(sum(v for b, v in c.items() if b != r) / sum(c.values())
                    for p, (r, c, _) in P.items() if p != POS and sum(c.values()) >= 200)
        x = cnt["A"] / d
        L.append(f"| {tag} | {d} | {cnt['C']} | {cnt['A']} | {cnt['G']} | {cnt['T']} | {x:.2%} | "
                 f"{sum(f >= x for f in fr)} of {len(fr)} | {fr[len(fr) // 2]:.2%} | {fr[int(len(fr) * 0.95)]:.2%} | "
                 f"{(d - cnt['C']) / d:.2%} | {al[len(al) // 2]:.2%} |")
    # the TCC reads and the other reads of the same molecules
    names = [n for n, b in res["all reads"][POS][2] if b == "A"]
    (W / "tcc_names.txt").write_text("\n".join(names) + "\n")
    rows = [l.split("\t") for l in sh(f"samtools view -N {W}/tcc_names.txt {dedup}").strip().split("\n") if l]
    near = [l.split("\t") for l in sh(f"samtools view {dedup} BK073075.1:{POS - 300}-{POS + 300}").strip().split("\n") if l]
    L += ["", "## The reads with TCC at the anticodon", "",
          "| read | mate covering 16416 reads | PCR copies of the same molecule (same start positions) read | mismatches in the TCC read | verdict |",
          "|---|---|---|---|---|"]
    for n in names:
        mine = [r for r in rows if r[0] == n]
        tcc = next(r for r in mine if base_at(r, POS) == "A")
        nm = next((x[5:] for x in tcc[11:] if x.startswith("NM:i:")), "?")
        mate = [r for r in mine if r is not tcc]
        mate_b = base_at(mate[0], POS) if mate else None
        key = (tcc[3], tcc[7])
        copies = [base_at(r, POS) for r in near if (r[3], r[7]) == key and r[0] != n]
        copies = [b for b in copies if b]
        mb = f"{ {'C': 'GCC', 'A': 'TCC'}.get(mate_b, mate_b)}" if mate_b else "does not cover it"
        cb = ", ".join(f"{anti[b]} {k}" for b, k in C.Counter(copies).items()) if copies else "none"
        contra = (mate_b == "C") or ("C" in copies)
        v = "sequencing error (same molecule reads GCC)" if contra else f"no second read of the molecule; {nm} mismatches in the read"
        L.append(f"| {n} | {mb} | {cb} | {nm} | {v} |")
    # ARAGORN relaxed, and nuclear copies
    ar = sh(f"{ROOT}/tools/bin/aragorn -t -ps -gc3 -c -l {MITO}")
    gly = re.findall(r"tRNA-Gly\((\w+)\)", ar)
    ntr = len(re.findall(r"tRNA-\S+\(", ar))
    paf = sh(f"minimap2 -x asm5 {nuc} {MITO} 2>/dev/null").strip().split("\n")
    hits = [p.split("\t") for p in paf if p]
    over = [h for h in hits if int(h[2]) < TRNA[1] and int(h[3]) > TRNA[0]]
    L += ["", "## Other copies", "",
          f"- ARAGORN 1.2.41, relaxed threshold (-ps), yeast mitochondrial code: {ntr} tRNAs in BK073075.1, of them "
          f"{len(gly)} tRNA-Gly (anticodon {', '.join(g.upper() for g in gly)}). tRNAscan-SE 2.0.12 (earlier run, "
          "`trnascan_BK073075.1.txt`) also finds one, GCC.",
          f"- Mitogenome against the nuclear assembly (minimap2 -x asm5): {len(hits)} match(es) in all: "
          + "; ".join(f"mito {h[2]}-{h[3]} to {h[5]}:{h[7]}-{h[8]} ({h[9]}/{h[10]} identical)" for h in hits)
          + f". Matches overlapping the tRNA-Gly gene: {len(over)}.", ""]
    (G / "trna_gly_tcc.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
