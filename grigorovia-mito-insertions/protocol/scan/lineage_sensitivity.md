# Real-data sensitivity in the known lineages (added after the freeze; descriptive)

Planted sensitivity for comparison (PREREG 1.12): a planted shift is counted at its codon 0.553 of the time (random sites) and 0.518 (same-codon trial sites). A genome with 2, 3 or 5 planted same-codon sites flags 0.317, 0.483 or 0.717 of the time.

| lineage | genomes scored | with documented sites | with >= 2 | documented sites | counted at its codon | genomes with >= 2 counted | flagged (any) | flagged, of those with >= 2 sites | outcomes of documented sites |
|---|---|---|---|---|---|---|---|---|---|
| CTG yeasts | 18 | 4 | 3 | 11 | 4 (0.36) | 1 | 6 | 2 of 3 | counted 4, clean_other_codon 3, flank_fail 2, misplaced (shift elsewhere in the gene) 1, detected clean in gene (codon not checked) 1 |
| glass sponges | 4 | 1 | 1 | 2 | 2 (1.00) | 1 | 1 | 1 of 1 | counted 2 |
| Perkinsus | 2 | 2 | 2 | 38 | 12 (0.32) | 2 | 2 | 2 of 2 | paired 14, counted 12, misplaced (shift elsewhere in the gene) 8, not_aligned 2, clean_other_codon 2 |
| birds | 1107 | 526 | 1 | 527 | 316 (0.60) | 1 | 2 | 1 of 1 | counted 316, misplaced (shift elsewhere in the gene) 130, homopolymer 47, flank_fail 13, clean_other_codon 8, clean_wide_window 6, aligned_no_event 5, paired 2 |
| turtles | 129 | 95 | 7 | 104 | 57 (0.55) | 1 | 2 | 2 of 7 | counted 57, homopolymer 16, misplaced (shift elsewhere in the gene) 16, clean_wide_window 7, flank_fail 3, clean_other_codon 2, detected clean at site (codon not checked) 2, paired 1 |
| oysters | 23 | 0 | 0 | 0 | 0 (0.00) | 0 | 0 | 0 of 0 |  |
| venerid clams | 28 | 0 | 0 | 0 | 0 (0.00) | 0 | 0 | 0 of 0 |  |
| Polyrhachis ants | 1 | 0 | 0 | 0 | 0 (0.00) | 0 | 0 | 0 of 0 |  |
| Arctica | 1 | 1 | 0 | 1 | 0 (0.00) | 0 | 0 | 0 of 0 | flank_fail 1 |
| Phaeodactylum | 1 | 1 | 0 | 1 | 0 (0.00) | 0 | 0 | 0 of 0 | homopolymer 1 |
| Dipterophagus | 1 | 1 | 0 | 1 | 0 (0.00) | 0 | 0 | 0 of 0 | misplaced (shift elsewhere in the gene) 1 |
| chytrids (Rhizophydium, Spizellomyces) | 4 | 2 | 0 | 2 | 0 (0.00) | 0 | 0 | 0 of 0 | no shift found in the gene 1, not_aligned 1 |
| other RefSeq single annotations | 2 | 2 | 0 | 2 | 0 (0.00) | 0 | 0 | 0 of 0 | clean_other_codon 2 |

### Real per-site rate by lineage (documented core-gene sites counted at their codon)

| lineage line | documented sites | counted | rate |
|---|---|---|---|
| Arctica | 1 | 0 | 0.00 |
| CTG yeasts | 10 | 4 | 0.40 |
| Perkinsus | 38 | 12 | 0.32 |
| birds (excl. nad3-174) | 14 | 2 | 0.14 |
| chytrids (Rhizophydium, Spizellomyces) | 1 | 0 | 0.00 |
| glass sponges | 2 | 2 | 1.00 |
| other RefSeq single annotations | 2 | 0 | 0.00 |
| turtles (excl. nad3-174) | 16 | 7 | 0.44 |
| nad3-174 (birds and turtles) | 598 | 364 | 0.61 |

**Headline real per-site rate, nad3-174 excluded: 0.29** (mean of per-lineage rates; bootstrap over 8 lineages, 95% 0.09-0.54; pooled 0.32 over 84 sites), against 0.553 for planted shifts.
Implied chance that a genome with m documented same-codon sites gets >= 2 counted (and then flags, since B5 conditional = 1.000): m = 2: 0.08 (0.01-0.29); m = 3: 0.20 (0.02-0.56); m = 5: 0.44 (0.07-0.86).

**Out-of-sample sites** (documented, not annotated in the record so not in B2, not in a named positive): 7 sites; detected as a clean shift (at the site, or in the gene for gene-level sites) 3 of 7 (0.43). Too few for a stable rate; listed in full:
- NC_004691.1 Nakaseomyces glabratus CBS 138: cox2 +1 [record note, 3' end (K15)]: detected clean in gene (codon not checked)
- NC_001947.1 Pelomedusa subrufa: nad3@9630 +1 [Russell & Beckenbach 2008 (K23; nad3-135)]: detected clean at site (codon not checked)
- NC_001947.1 Pelomedusa subrufa: nad4l@10015 +1 [Russell & Beckenbach 2008 (K23; 'near nt 99')]: detected clean at site (codon not checked)
- NC_001947.1 Pelomedusa subrufa: nad4l@10178 +1 [Russell & Beckenbach 2008 (K23; 'near nt 262')]: flank_fail
- NC_016739.1 Phaeodactylum tricornutum: cox1 +1 [Oudot-Le Secq & Green 2011 (K29)]: homopolymer
- NC_058598.1 Dipterophagus daci: nad5@7189 -1 [Towett-Kirui 2022 (K31)]: misplaced (shift elsewhere in the gene)
- NC_003052.1 Spizellomyces punctatus: cob +1 [record note (K33)]: no shift found in the gene
Write-up limits (PREREG log): a negative rules out lineages only at these real-data rates; the pair filter removes real sites where +1 and +2 shifts sit close together; the chance figures assume errors at random codons, and assembly errors cluster, so a noisy genome can be strong without being real.

Genomes with no documented sites are not misses. Per-genome detail for the small lineages:

- CTG yeasts: Huiozyma sinensis (NC_036748.1): flagged +1CTG; cox3@16878 [record annotation]: counted; atp6@29007 [record annotation]: counted; cob@34612 [record annotation]: misplaced (shift elsewhere in the gene)
- CTG yeasts: Monosporozyma unispora (NC_036378.1): flagged no; atp6@7780 [record annotation]: clean_other_codon; cox2@14542 [record annotation]: counted; rps3@26764 [record annotation]: non-core gene (not tested)
- CTG yeasts: Nakaseomyces nivariensis (NC_036379.1): flagged no; no documented sites
- CTG yeasts: Nakaseomyces bracarensis (NC_036380.1): flagged no; no documented sites
- CTG yeasts: Nakaseomyces glabratus CBS 138 (NC_004691.1): flagged no; cox2 [record note, 3' end (K15)]: detected clean in gene (codon not checked)
- CTG yeasts: Nakaseomyces delphensis (NC_012619.1): flagged no; no documented sites
- CTG yeasts: Monosporozyma servazzii (NC_004918.1): flagged +1CTG; cox2@4697 [record annotation]: clean_other_codon; cob@8588 [record annotation]: flank_fail; cob@8610 [record annotation]: flank_fail; cox3@18298 [record annotation]: counted; atp6@28460 [record annotation]: clean_other_codon
- CTG yeasts: Huiozyma naganishii (BK073131.1): flagged +1CTG; no documented sites
- CTG yeasts: Kazachstania bromeliacearum (BK073182.1): flagged no; no documented sites
- CTG yeasts: Kazachstania kunashirensis (BK073183.1): flagged no; no documented sites
- CTG yeasts: Kazachstania taianensis (BK073184.1): flagged no; no documented sites
- CTG yeasts: Kazachstania viticola (BK073185.1): flagged +1CTG; no documented sites
- CTG yeasts: Monosporozyma aerobia (BK073484.1): flagged +1CTG; no documented sites
- CTG yeasts: Monosporozyma siamensis (BK073485.1): flagged no; no documented sites
- CTG yeasts: Monosporozyma solicola (BK073486.1): flagged +1CTG; no documented sites
- CTG yeasts: Nakaseomyces cf. glabratus (BK073615.1): flagged no; no documented sites
- CTG yeasts: Nakaseomyces kungkrabaensis (BK073616.1): flagged no; no documented sites
- CTG yeasts: Nakaseomyces uthaithaninus (BK073617.1): flagged no; no documented sites
- glass sponges: Aphrocallistes vastus (NC_010769.1): flagged +1TGG; cox3@3513 [record annotation]: counted; nad6@11256 [record annotation]: counted
- glass sponges: Caulophacus iocasicus (NC_072601.1): flagged no; no documented sites
- glass sponges: Oopsacas minuta (NC_027419.1): flagged no; no documented sites
- glass sponges: Vazella pourtalesii (NC_028054.1): flagged no; no documented sites
- Perkinsus: Perkinsus marinus (NC_068742.1): flagged +1AGG; cox1@1589 [record annotation]: paired; cox1@1608 [record annotation]: paired; cox1@1748 [record annotation]: paired; cox1@1827 [record annotation]: counted; cox1@2404 [record annotation]: counted; cox1@2447 [record annotation]: counted; cox1@2526 [record annotation]: misplaced (shift elsewhere in the gene); cox1@2593 [record annotation]: misplaced (shift elsewhere in the gene); cox1@2648 [record annotation]: misplaced (shift elsewhere in the gene); cox1@2664 [record annotation]: misplaced (shift elsewhere in the gene); cox3@15034 [record annotation]: not_aligned; cob@27089 [record annotation]: paired; cob@27132 [record annotation]: paired; cob@27212 [record annotation]: paired; cob@27414 [record annotation]: counted; cob@27445 [record annotation]: counted; cob@27476 [record annotation]: counted; cob@27669 [record annotation]: clean_other_codon; cob@27925 [record annotation]: misplaced (shift elsewhere in the gene)
- Perkinsus: Perkinsus olseni (NC_068743.1): flagged +1AGG; cox1@8483 [record annotation]: paired; cox1@8502 [record annotation]: paired; cox1@8642 [record annotation]: paired; cox1@8721 [record annotation]: counted; cox1@9298 [record annotation]: counted; cox1@9341 [record annotation]: counted; cox1@9420 [record annotation]: misplaced (shift elsewhere in the gene); cox1@9487 [record annotation]: misplaced (shift elsewhere in the gene); cox1@9542 [record annotation]: paired; cox1@9558 [record annotation]: paired; cox3@19548 [record annotation]: not_aligned; cob@61789 [record annotation]: paired; cob@61746 [record annotation]: paired; cob@61666 [record annotation]: paired; cob@61464 [record annotation]: counted; cob@61433 [record annotation]: counted; cob@61402 [record annotation]: counted; cob@61209 [record annotation]: clean_other_codon; cob@60953 [record annotation]: misplaced (shift elsewhere in the gene)

Birds and turtles: full per-genome list in lineage_sensitivity.tsv.
