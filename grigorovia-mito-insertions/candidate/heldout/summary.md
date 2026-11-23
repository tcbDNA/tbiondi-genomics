# Held-out relative: Grigorovia yakushimaensis OR703946.1 (frozen pipeline, unchanged)

- record: UNVERIFIED: Grigorovia yakushimaensis mitochondrion sequence; 49874 bp; CDS features: 0; comment: GenBank staff is unable to verify sequence and/or annotation provided by the submitter.
- core genes aligned: atp6,atp9,cob,cox1,cox2,cox3,other:orf146,other:rps5; breaks: {'flank_fail': 3, 'gene_end': 1}; clean core shifts: 0; flags: 0; clean +1 at AGG: 0
- site 24711: 120 nt aligned at 0.90; automatic call: lacks
- site 26496: 68 nt aligned at 0.88; automatic call: carries
- site 26477: 56 nt aligned at 0.91; automatic call: carries
- site 26266: 120 nt aligned at 0.90; automatic call: lacks
- tRNA-Gly (ARAGORN): tcc

## Checked by eye (coding strand; the automatic calls above are wrong at 26496 and 26477)

- **26496 (AGG site): lacks.** *G. jiainica* `...TTGTGTATAT AGGCCAAAT` aligns to *G. yakushimaensis*
  `...TTGCGTATA- CGGGCAAAT` (OR703946 21200-21268). The gap sits 1 nt from the extra base in a TA repeat.
- **26477 (TGG site): not comparable.** In *G. yakushimaensis*, cob is interrupted here by a 1,026-nt insert that the
  pipeline infers as an intron (OR703946 21737-22763). Its downstream exon (from 22763) matches *G. jiainica* starting
  4 nt after the TGG site.
- **24711 and 26266: lacks** (120 nt aligned at 90% identity).

**Outcome against the prediction written before the record was opened** (protocol/frozen_rules.md):
- The predicted outcome was not met. No homologous site carries the extra base, and there are no clean +1 shifts at
  AGG. By the logged reading this means "sites probably recent in *G. jiainica*, not a refutation".
- tRNA-Gly is UCC, as predicted for a genome that does not skip. tRNAscan-SE -O gives score 56.1.
- Record: UNVERIFIED, linear, 49,874 bp, no CDS features. Submitted by Wolters, LaBella, Opulente, Rokas and Hittinger
  (UW-Madison, 2023). Assembled with iWGS 1.1, NOVOPlasty 4.2 and plasmidSPAdes 3.9.0 from Illumina reads; BioSample
  SAMN08343276.
