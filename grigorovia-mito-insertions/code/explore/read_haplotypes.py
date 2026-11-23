"""Group reads by the exact sequence they carry across a reference window, using each read's alignment (CIGAR),
not sequence anchors (anchors miss the other allele when it differs inside them).
Only primary alignments with MAPQ >= 20 that cover the window plus 5 nt each side, with no clipping inside, count.

    uv run python code/explore/read_haplotypes.py BAM CONTIG START END   (1-based, inclusive)
"""
import re
import subprocess
import sys
from collections import Counter


def haplotypes(bam, contig, a, b, pad=5, minq=20):
    out = subprocess.run(["samtools", "view", "-F", "0x904", "-q", str(minq), bam, f"{contig}:{a - pad}-{b + pad}"],
                         capture_output=True, text=True, check=True).stdout
    hap = Counter()
    for line in out.splitlines():
        f = line.split("\t")
        pos, cigar, seq = int(f[3]), f[5], f[9]
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        r, q = pos, 0
        qa = qb = None
        ok = True
        for n, op in ops:
            n = int(n)
            if op in "M=X":
                for i in range(n):
                    if r + i == a:
                        qa = q + i
                    if r + i == b:
                        qb = q + i
                r += n
                q += n
            elif op == "I":
                q += n
            elif op in "DN":
                if r <= a < r + n or r <= b < r + n:
                    ok = False
                r += n
            elif op == "S":
                if q > 0 and r <= b + pad:
                    ok = ok and r > b + pad
                q += n
        start_ok = pos <= a - pad
        end_ok = r - 1 >= b + pad
        if ok and start_ok and end_ok and qa is not None and qb is not None:
            hap[seq[qa:qb + 1]] += 1
    return hap


if __name__ == "__main__":
    bam, c, a, b = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    h = haplotypes(bam, c, a, b)
    tot = sum(h.values())
    print(f"{c}:{a}-{b}: {tot} reads")
    for s, n in h.most_common(8):
        print(f"  {n:4d} ({n / tot:.0%})  len {len(s)}  {s}")
