# Post-hoc outputs for cold review 2 (2026-11-25) [P]

Everything here is post-hoc. No output here changes a locked output, PREREG.md or a gate report. REPORT.md, the figures and the dossiers cite these files as [P]. Scripts: `src/cr2_*`, `src/tyscan.c`. Reads and alignments are in `results/work/r2/` (not tracked); downloads are ledger rows with gate `6x`.

| Check | Outputs |
|---|---|
| 1 locus gate | `locus_gate_corroboration.tsv`, `locus_gate_retry_b.tsv`, `locus_gate_retry_ab.tsv`, `locus_gate_summary.txt`; `bufotes_framing.txt` |
| 2 *Craniella* | `craniella_depth.txt`, `craniella_depth2.txt`, `craniella_depth3.txt` (mitogenome written twice: use these counts), `craniella_rna.txt` |
| 3 starfish | `starfish_li2019.txt`, `starfish_asteroids.tsv` |
| 4 second copies | `secondcopy/summary.tsv`, one text file per run, `sensitivity.tsv`; high-copy contigs: `highcopy_summary.txt` (use this), `highcopy_*.txt` per assembly, `highcopy_craniella.txt` (windows of the deposited assemblies), `beanii_contig.txt` |
| 5 read power | `read_power.tsv`, `read_power_check4.tsv`, `read_power_tiers.txt` |
| 6 *Rotylenchulus* | `roty_criteria.md` (written 2026-11-25 07:26, before mapping), `roty_eval.txt` |
| 7 titles | `title_audit.tsv`, `title_audit_summary.txt` |
| 8 code and counts | `pvulnus_models.txt`, `cox1_reconcile.txt`, `numt_truecov.txt`, `repro_{b1,cols,stem,elig}.txt`, `fig3_bootstrap.tsv`, `weak_audit_strict.txt` |
| Pre-push fixes (2026-11-25) | `conspecific_specimen_audit.tsv`, `conspecific_specimen_notes.txt` (Tier 1 specimens and panel *A. vulpes*); `mtrf1a_stop_denominator.tsv` (the 47 UAA=Stop species) |
| Preprint, six-fish stops (2026-11-25) | `sixfish_polyA_criteria.md` (rules, written before any RNA read was looked at, with two dated notes); `sixfish_terminal_uaa.tsv` (per species, Table 3 of the preprint), `sixfish_terminal_genes.tsv` (per gene and annotation), `sixfish_terminal_refs.tsv` (14 reference fish), `sixfish_polya_junctions.tsv` and `sixfish_polya_positions.tsv` (poly(A) ends in *A. glossodonta* and *L. verrillii* RNA). Scripts `src/preprint_sixfish_{pats,polya,stops}.py` |
| Preprint, number trace (2026-11-25) | `preprint_number_trace.tsv`: every number in the manuscript text, tables, captions, references and figure labels, with its location, source file, match (yes; no; label, for part of a name such as T216; untraced), how it was checked, the surrounding words and a note; `preprint_refs_check.tsv`: each reference against its Europe PMC, Crossref or Zenodo record. Scripts: `src/preprint_trace/` (`extract.py`, `refs.py`, `run.py`, in that order). The scripts take the path of the LaTeX file in `manuscript/`, which was added after the trace was made. |
| Preprint, context (2026-11-25) | `elopomorpha_c33.tsv` (C33 over resolved UAA=Stop vertebrate genomes and minimum gains by the T3 counter, Elopomorpha against the rest; `src/preprint_elopomorpha.py`); `tylench_trny_identity.tsv` (tRNA-Tyr identities by boundary set; `src/preprint_trny_identity.py`); `prior_art_c33_fish.md` (literature search) |
