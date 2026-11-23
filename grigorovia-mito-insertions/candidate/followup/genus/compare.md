# Test 6: the four Grigorovia mitogenomes compared

## Frozen scan pipeline (unchanged) on each mitogenome

| species | record | length | core genes aligned | clean core shifts | +1 at AGG (clean) | flags | stop columns in core genes |
|---|---|---|---|---|---|---|---|
| G. jiainica | BK073075.1 | 29,214 | atp6, atp9, cob, cox1, cox2, cox3 | 4 | 3 | 1 | 0 |
| G. yakushimaensis | OR703946.1 | 49,874 | atp6, atp9, cob, cox1, cox2, cox3 | 0 | 0 | 0 | 0 |
| G. humatica | humatica_NODE_74_length_53413_cov_109.893787, humatica_NODE_ | 62,112 | atp6, atp8, atp9, cob, cox1, cox2, cox3 | 0 | 0 | 0 | 0 |
| G. transvaalensis | transvaalensis_bait_NODE_1_length_15361_cov_93.661415, trans | 37,049 | atp6, atp8, atp9, cob, cox1, cox2, cox3 | 0 | 0 | 0 | 0 |

## Do the other species carry G. jiainica's extra bases?

G. jiainica window of +-60 nt aligned locally (both strands) to each species' spliced CDS of the same gene (miniprot exon structure, introns removed). Across 15 nt either side of the extra base, 'lacks' = G. jiainica has one base more than the other species (net difference 1, mod 3); 'same length' = no net difference, i.e. the other species has a base there too; 'not spanned' = the alignment does not cover that stretch.

| site | G. yakushimaensis | G. humatica | G. transvaalensis |
|---|---|---|---|
| 24711 | lacks (120 nt, 0.90) | lacks (120 nt, 0.93) | lacks (120 nt, 0.88) |
| 26496 | lacks (107 nt, 0.87) | lacks (105 nt, 0.83) | lacks (119 nt, 0.82) |
| 26477 | lacks (106 nt, 0.84) | lacks (104 nt, 0.81) | lacks (118 nt, 0.80) |
| 26266 | lacks (120 nt, 0.90) | lacks (120 nt, 0.89) | lacks (119 nt, 0.91) |

## cox1 and cob introns (miniprot exon structures; junctions mapped onto G. jiainica, extra bases removed)

| species | gene | exons (oriented) | introns: position after G. jiainica coding nt (nearest site, distance) |
|---|---|---|---|
| G. jiainica | cox1 | BK073075.1: 4273-5863 (protein 1-534) | no intron |
| G. jiainica | cob | BK073075.1: 2311-3465 (protein 1-384) | no intron |
| G. yakushimaensis | cox1 | OR703946.1: 10700-11714, 13327-13943 (protein 1-534) | after 1029 (cox1 24711: 742 nt) |
| G. yakushimaensis | cob | OR703946.1: 20853-21279, 22771-23489 (protein 1-385) | cob 26477: 2 nt (refined; equal-score splits 2-3 nt) |
| G. humatica | cox1 | humatica_NODE_74_length_53413_cov_109.893787: 32358-33283, 35787-35858, 38042-38660 (protein 1-534) | after 982 (cox1 24711: 695 nt); unmapped |
| G. humatica | cob | humatica_NODE_74_length_53413_cov_109.893787: 39470-39896, 41393-42111 (protein 1-385) | cob 26477: 2 nt (refined; equal-score splits 2-3 nt) |
| G. transvaalensis | cox1 | transvaalensis_bait_NODE_2_length_8313_cov_84.941730: 766-1444, 3143-3570, 6738-7220 (protein 1-534) | after 735 (cox1 24711: 448 nt); unmapped |
| G. transvaalensis | cob | transvaalensis_bait_NODE_1_length_15361_cov_93.661415: 2937-4091 (protein 1-385) | no intron |

## tRNA-Gly (tRNAscan-SE -O and ARAGORN)

| species | tRNAscan-SE tRNAs | tRNA-Gly anticodons (tRNAscan-SE) | ARAGORN tRNA-Gly |
|---|---|---|---|
| G. jiainica | 23 | GCC | gcc |
| G. yakushimaensis | 23 | TCC | tcc |
| G. humatica | 23 | TCC | tcc |
| G. transvaalensis | 16 | TCC | tcc |
