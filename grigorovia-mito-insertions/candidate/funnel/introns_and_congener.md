# Intron positions and the congener at the four sites (post hoc; prompted by the referee review)

## Intron positions in *S. cerevisiae* (NC_001224, CDS joins from RefSeq)

Extra-base positions are given in *S. cerevisiae* CDS coordinates: just before the codon of the residue read after the
skip.

| site (BK073075.1) | residue after the extra base (S. cerevisiae numbering) | extra base between CDS nt | nearest S. cerevisiae intron (inserted after CDS nt) | distance |
|---|---|---|---|---|
| cob 26496 (AGG) | Gly137 | 408 / 409 | bI1, after nt 415 (within codon 139) | 7 nt |
| cob 26477 (TGG) | Gly143 | 426 / 427 | bI2, after nt 429 (end of codon 143) | 3 nt |
| cox1 24711 (AGG) | residue 78 | 231 / 232 | aI3, after nt 243 (end of codon 81); aI2 after nt 205 | 12 nt; 26 nt |
| cob 26266 (ACT/AGG) | residue 214 | 639 / 640 | after nt 506 and after nt 756 | 133 nt; 117 nt |

- S. cerevisiae cob introns are inserted after CDS nt 415, 429, 506, 756 and 807.
- Its cox1 introns are inserted after nt 169, 205, 243, 720, 972, 1107 and 1132.
- In *G. yakushimaensis* (OR703946), cob is interrupted by an inferred 1,026-nt intron at the TGG site (candidate/heldout/).
- G. jiainica has no intron at any of these positions in cox1 or cob.

**Three of the four extra bases sit within 12 nt of an intron insertion position in S. cerevisiae.** One of them is
also at an intron in the congener.

## The congener *G. yakushimaensis* at the sites (local alignment of the coding strand)

- **cox1 24711.**
  - *G. jiainica*: `ATT GGA AGG CTT TGG TAA`. With the A of AGG skipped: `ATT GGA GGC TTT GGT AAT` = Ile Gly **Gly** Phe Gly
    Asn.
  - *G. yakushimaensis*: `ATT GGT AGT TTT GGT AAT` = Ile Gly **Ser** Phe Gly Asn.
  - **The congener has Ser where the skip-A reading gives Gly.**
  - Skipping instead the G of AGG (`ATT GGA AGC TTT GGT`) gives Ser and fits the congener equally well or better. The
    position of this extra base is ambiguous, and the most-related sequence favours a reading that is not AGG→GGC.
  - The distant reference proteins mostly have Gly here (8 of 10).
- **cob 26266.**
  - *G. jiainica*: `GGA ATT [A]CT AGG CAA TTT`.
  - *G. yakushimaensis*: `GGT ATC CTT GGT AAT TTA`.
  - The congener has Thr (CTT, code 3) then Gly (GGT), the same two residues as either placement in G. jiainica.
  - Its codon structure (`CTT GGT`) resembles the ACT placement (`CTA GGC`).
- **cob 26496 and 26477.** See candidate/heldout/summary.md: the extra A is absent in the congener, and the TGG site is
  interrupted by its intron.
- **Glycine codons at the GGC positions in the six closest relatives** (blinded check; genome_1 = *K. bromeliacearum*,
  genome_3 = *M. solicola*): GGA, GGT, GGT, GGT and GGT, GGA, GGA, GGA. In *G. jiainica* the GGC is derived.
