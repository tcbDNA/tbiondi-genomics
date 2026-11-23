# Test 5: the rest of the mitogenome (G. jiainica BK073075.1 and the congener OR703946.1)

## Protein genes (miniprot 0.18, -T 3; S. cerevisiae proteins; best alignment per gene; 'relaxed' = --outs=0.3 -p 0.1 -n 1 -l 3; VAR1 is assessed by its open frame below, as it is too divergent for miniprot)

| genome | gene | settings | query length | aligned | target | frameshifts | stops |
|---|---|---|---|---|---|---|---|
| G. jiainica | atp6 | default | 259 | 1-259 | 936-1709 (-) | 0 | 0 |
| G. jiainica | atp8 | default | 48 | 1-48 | 1758-1901 (-) | 0 | 0 |
| G. jiainica | atp9 | default | 76 | 1-76 | 27697-27924 (+) | 0 | 0 |
| G. jiainica | cob | default | 385 | 1-384 | 25750-26904 (-) | 3 | 0 |
| G. jiainica | cox1 | default | 534 | 1-534 | 23355-24942 (-) | 1 | 0 |
| G. jiainica | cox2 | default | 251 | 1-251 | 2624-3376 (+) | 0 | 0 |
| G. jiainica | cox3 | default | 269 | 1-269 | 22079-22885 (-) | 0 | 0 |
| congener | atp6 | default | 259 | 1-259 | 13342-14118 (+) | 0 | 0 |
| congener | atp8 | relaxed | 48 | 1-48 | 13001-13144 (+) | 0 | 0 |
| congener | atp9 | default | 76 | 1-76 | 29144-29371 (-) | 0 | 0 |
| congener | cob | default | 385 | 1-385 | 20853-23489 (+) | 0 | 0 |
| congener | cox1 | default | 534 | 1-534 | 35935-39175 (-) | 0 | 0 |
| congener | cox2 | default | 251 | 1-251 | 17914-18663 (-) | 0 | 0 |
| congener | cox3 | default | 269 | 1-269 | 11478-12284 (+) | 0 | 0 |

## VAR1 (rps3)

- G. jiainica: other:rps3 hit 19206-19868 (bits 86); longest open frame over it: 259 codons (strand +1)
- congener: other:VAR1 hit 25329-26006 (bits 84); longest open frame over it: 272 codons (strand +1)

## RNA genes (local alignment of S. cerevisiae sequences; match 2, mismatch -3, gap -5/-2)

| genome | gene (S. cerevisiae length) | score | strand | coverage of the S. cerevisiae sequence | identity |
|---|---|---|---|---|---|
| G. jiainica | 15S_RRNA (1649) | 763 | -1 | 0.79 | 0.82 |
| G. jiainica | 21S_RRNA (3296) | 1000 | -1 | 0.79 | 0.83 |
| G. jiainica | RPM1 (483) | 107 | -1 | 0.53 | 0.82 |
| congener | 15S_RRNA (1649) | 681 | +1 | 0.83 | 0.84 |
| congener | 21S_RRNA (3296) | 1051 | -1 | 0.72 | 0.82 |
| congener | RPM1 (483) | 106 | +1 | 0.53 | 0.79 |

## Every scan event and stop column in G. jiainica

| type | gene | forward position (0-based) | shift / residue | codon | class |
|---|---|---|---|---|---|
| event | cox1 | 24710 | 1 | AGG | clean |
| event | cob | 26495 | 1 | AGG | clean |
| event | cob | 26476 | 1 | TGG | clean |
| event | cob | 26265 | 1 | AGG | clean |
| event | other:orf146 | 11167 | 1 | CGC | flank_fail |
| stop column | other:orf146 | 18167 | I | | |

## tRNAs (tRNAscan-SE 2.0.12 -O)

G. jiainica 23 tRNAs, congener 23; notes (pseudogene etc.): 0.

| type | anticodon | G. jiainica (scores) | congener (scores) |
|---|---|---|---|
| Ala | TGC | 1 (55.9) | 1 (50.9) |
| Arg | TCT | 1 (55.0) | 1 (59.6) |
| Asn | GTT | 1 (63.0) | 1 (65.6) |
| Asp | GTC | 1 (47.2) | 1 (54.2) |
| Cys | GCA | 1 (26.7) | 1 (35.1) |
| Gln | TTG | 1 (41.5) | 1 (43.6) |
| Glu | TTC | 1 (55.5) | 1 (61.6) |
| Gly | GCC | 1 (52.0) | 0 () |
| Gly | TCC | 0 () | 1 (56.1) |
| His | GTG | 1 (46.3) | 1 (50.2) |
| Ile | GAT | 1 (55.6) | 1 (56.0) |
| Leu | TAA | 1 (47.0) | 1 (48.1) |
| Lys | TTT | 1 (60.1) | 1 (66.1) |
| Met | CAT | 2 (58.3, 44.5) | 2 (53.6, 55.6) |
| Phe | GAA | 1 (52.9) | 1 (59.5) |
| Pro | TGG | 1 (48.2) | 1 (35.3) |
| Ser | GCT | 1 (38.0) | 1 (46.0) |
| Ser | TGA | 1 (64.9) | 1 (63.3) |
| Sup | TCA | 1 (57.9) | 1 (53.0) |
| Sup | TTA | 1 (37.9) | 1 (47.1) |
| Thr | TGT | 1 (61.2) | 1 (59.0) |
| Tyr | GTA | 1 (45.9) | 1 (39.8) |
| Val | TAC | 1 (59.9) | 1 (61.3) |

Score ranges: G. jiainica 26.7-64.9; congener 35.1-66.1.
