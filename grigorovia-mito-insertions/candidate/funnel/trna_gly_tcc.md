# Is there a TCC tRNA-Gly in G. jiainica's mitochondrion? (post hoc, exploratory)

Reads SRR16974510 (5235450 reads; 1470627 marked as duplicates) mapped to BK073075.1 plus the nuclear assembly GCA_030571715.1 (minimap2 -x sr). Base quality >= 20, mapping quality >= 20. Position 16416 (+ strand) is anticodon position 34 of the minus-strand tRNA-Gly (c[16378,16449]).

| reads | depth at 16416 | GCC | TCC | CCC | ACC | TCC share | other C positions with a C>A share at least as high | median C>A share elsewhere | 95th percentile | any-mismatch share at 16416 | median any-mismatch share elsewhere |
|---|---|---|---|---|---|---|---|---|---|---|---|
| all reads | 1249 | 1243 | 4 | 2 | 0 | 0.32% | 74 of 564 | 0.15% | 0.42% | 0.48% | 0.21% |
| duplicates removed | 911 | 908 | 2 | 1 | 0 | 0.22% | 216 of 564 | 0.14% | 0.50% | 0.33% | 0.21% |

## The reads with TCC at the anticodon

| read | mate covering 16416 reads | PCR copies of the same molecule (same start positions) read | mismatches in the TCC read | verdict |
|---|---|---|---|---|
| SRR16974510.322023 | GCC | none | 1 | sequencing error (same molecule reads GCC) |
| SRR16974510.1889620 | does not cover it | GCC 2 | 2 | sequencing error (same molecule reads GCC) |
| SRR16974510.2313667 | does not cover it | GCC 6 | 1 | sequencing error (same molecule reads GCC) |
| SRR16974510.833658 | does not cover it | none | 8 | no second read of the molecule; 8 mismatches in the read |

## Other copies

- ARAGORN 1.2.41, relaxed threshold (-ps), yeast mitochondrial code: 24 tRNAs in BK073075.1, of them 1 tRNA-Gly (anticodon GCC). tRNAscan-SE 2.0.12 (earlier run, `trnascan_BK073075.1.txt`) also finds one, GCC.
- Mitogenome against the nuclear assembly (minimap2 -x asm5): 1 match(es) in all: mito 18903-18974 to JAJLYE010000938.1:6135-6208 (43/73 identical). Matches overlapping the tRNA-Gly gene: 0.

