# candidate/

The exploratory follow-up of *Grigorovia jiainica*, the one candidate. Every number about it in `../README.md` and
`../notebook.md` comes from one of these files.

- `grigorovia_jiainica.md`: the six follow-up tests, results and reading. Predictions for each are in
  `../protocol/frozen_rules.md` (log), written before it ran.

## funnel/ (the steps after the result)
- `step_a_sites.md` … `step_f_reads.md`: the funnel steps in their logged order.
- `step_f_reads_dedup.md`: the step f read counts with PCR duplicates removed, also for the two close *cob*
  sites (added after an outside review).
- `agg_usage_relatives.md`: in-frame AGG and glycine codon use in *G. jiainica* and its six closest relatives.
- `closest_relatives.md`, `closest_relatives.json`: relatives ranked by protein identity.
- `trna_gly_structure.md`, `fam_trna_gly.json`: tRNA-Gly loops, base 37 and anticodon read support; tRNA-Gly
  anticodons in all 101 Saccharomycetaceae.
- `trna_gly_tcc.md`: whether a TCC tRNA-Gly exists: read counts at the anticodon with and without duplicates, the
  background error rate, each TCC read against other reads of the same molecule, and searches for other copies
  (added after an outside review).
- `completeness_and_data.md`: completeness of the mitogenome, and the other public data on the species and genus.
- `introns_and_congener.md`: intron positions near the sites (*S. cerevisiae*), and the congener's codons at them
  (post hoc, prompted by the referee).
- `sites_for_lab.md`: strain, culture collection, coordinates and sequence context of the four sites.
- `blind/`: the blinded re-check: its task (`TASK.md`) and verdict (`VERDICT.md`); `blind_key.json` decodes the
  labels.
- `figure_site.json`: the data drawn in `../plots/`: codons, code-3 residues, *S. cerevisiae* residues, read counts,
  and the anticodons of the tRNAs the *G. jiainica* mitogenome encodes (tRNAscan-SE).

## followup/ (six follow-up tests: bypassed or broken)
- `test1_*`: the congener's codons at the four sites. `test1_sister_codons.tsv` has one row per site;
  `test1_placements.md` has every candidate placement with its score and residues.
- `test2_*`: introns.
  - `test2_introns.md`: the null (share of coding positions near intron positions) and the observed distances.
  - `test2b_annotated_introns_on_gj.tsv`: annotated introns mapped onto *G. jiainica*.
  - `test2b_relatives_miniprot.md`: exon structures in the congener and the six closest relatives.
  - `test2b_intron_orfs.md`: intron-ORF remnants, the congener's intron and the Rfam check.
  - `test2b_servazzii_cob.md`: *M. servazzii*'s cob skips.
  - `test2_summary.md`: the reading.
- `test3_nuclear_genes.md` and `nuclear/`: nuclear genes.
  - `pass1_calls.*`: *S. cerevisiae* queries.
  - `close_queries.*`: relatives' queries.
  - `congener_query_calls.*`: congener queries.
  - `reads_check.md`: read consensus, read hits and the PET122 local assembly.
  - `heterozygosity.md`: genome-wide heterozygosity, heterozygous variants and read versions at each flagged spot.
  - `rates.md`: identity to the congener.
  - `rbh.md`: reciprocal best hits against the *S. cerevisiae* proteome (added after review).
  - `yak_queries_sources.tsv`: where each congener query protein came from.
  - The query proteins: `sc_nuclear_panel.faa`, `close_queries.faa`, `yak_queries.faa`.
- `completeness/`: compleasm results per assembly, with one reference protein per BUSCO gene.
- `test4_phenotypes.md`: growth on non-fermentable carbon in the congeners' descriptions and trait tables.
- `test5_*`: the rest of the *G. jiainica* and congener mitogenomes (genes, VAR1, rRNAs, tRNAs, scan events).
- `reads_strand.md`: reads at the four sites by strand, and reads spanning both close *cob* sites (added after review).
- `test6_genus.md` and `genus/`: the four *Grigorovia* mitogenomes.
  - `mtreads.md`: the share of mtDNA reads per run.
  - `pick.md`, `bait.md`, `confirm_*.md`, `genes_in_contigs.md`: assembly and confirmation.
  - `humatica_mito.fna`, `transvaalensis_mito.fna`: the assembled mitochondrial contigs.
  - `compare.md`: the frozen scan, the four sites, introns and tRNA-Gly.
  - `refine_cob_junction.md`: refined cob intron junctions.
  - `codon_alignments/`, `dnds.md`/`.json`, `dnds_profile.md`/`_jiainica.json`: dN/dS.
  - `dnds_sensitivity.md`/`.json`: the same test under other codon models and without outgroups (added after review).
  - `ggc_sites.md`: GGC codons in *G. jiainica*'s core genes and every taxon's codons at the four sites (added after
    review).
  - `codon_positions.md`: codons per gene, glycine codons by type and third-position bases (added after an outside
    review).
  - `trnascan_*.txt`: tRNAscan-SE output for the assembled mitogenomes.
  - `downloads.md`: what was downloaded and its size.

## heldout/
- `summary.md`, `genomes.tsv`, `events.tsv`, `inserts.tsv`: *Grigorovia yakushimaensis* (OR703946), run through the
  frozen pipeline after a written prediction.

## search_log/
- `grigorovia_literature.md`, `grigorovia_epmc_query_grid.json`: the *Grigorovia* literature and data search, and
  its 140-query grid.
