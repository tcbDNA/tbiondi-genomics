# Design decisions (agreed 2026-11-21 19:11 EST, before any data were downloaded or any code written)

**Question.** Does any mitochondrial lineage shift reading frame at a recurring codon, the way
glass-sponge mitochondria do at UGG? Single known sites (bird and turtle *nad3*, oyster *cytb*,
*Polyrhachis* ants) are controls. A lineage whose shifts pile up at one codon is the prize. A clean
"no, and here is what the scan could have seen" is also a result.

The method follows ry-xu/genome-research `cgg-skip-tanb77` (README, notes, data/prereg.md); its code is
not reused. Order of work, each step written up before the next: 1 prior art (`PRIOR_ART.md`),
2 controls before data, 3 frozen rules (`PREREG.md`), 4 scan and kill funnel, 5 blinded independent
check, 6 write-up. Python with uv. Ask before downloading > 20 GB or starting anything that runs > 1 h.

## Data
- Every RefSeq mitochondrial genome (release 238 at the time of writing) and the 2,695 fungal
  mitogenomes of Sarhan et al. 2025 (*Sci Data*, doi:10.1038/s41597-025-06447-x; TPA BK072095–BK074789),
  each translated with its own mitochondrial genetic code.
- Raw reads and RNA-seq from SRA, for candidates only.
- GenBank mitogenomes of species not in the scan are a **held-out replication set** (item 7).

## Decisions (numbers follow the interview)
1. **Single-genome flag:** ≥ 2 clean shifts in different genes at the same codon, where that codon is
   rare in frame (catches *Aphrocallistes vastus* alone).
2. **Lineage pools:** NCBI family **and** order, corrected for the two tests. A shift at the same
   orthologous site counts once per pool. Where NCBI gives no family (or order), use the nearest
   named clade above genus.
3. **Shift types:** every net frame change, reported by offset (+1 or −1; +2 ≡ −1, −2 ≡ +1).
4. **Confirmatory test: the codon route** (old-frame codon at the shift). The 4-base motif route runs
   on pools only, needs ≥ 8 distinct sites, and is **exploratory**.
5. **Main scan: everything except land plants** (Embryophyta); protists stay in (*Perkinsus*). Land
   plants are scanned afterwards as a separate stratum with its own null, **exploratory**.
6. **RNA-editing lineages** (*Physarum* and other myxomycetes, kinetoplastids) stay in as
   **known-answer tests for the RNA step at single sites**, not as lineages that must flag. If they do
   not flag, log why; do not tune the scan until they do.
7. **Held-out set:** species in GenBank but not in the scan. Second assemblies of scanned species stay
   the step-4 source. Before opening any held-out relative of a flagged lineage, write the prediction
   in `PREREG.md`. After the main result is frozen, scan the rest of the held-out set as a second,
   labelled pass.
8. **Aligners:** scan with DIAMOND blastx in frameshift mode under each genome's own code; the blinded
   check uses miniprot.
9. **References:** RefSeq mitochondrial proteins plus the fungal set, **excluding the genome's own
   family**.
10. **Introns:** a break within 30 nt of an annotated intron edge is dropped, as is any break that
    spans > 30 nt of unaligned sequence.
11. **Control bars** (one script, exits non-zero if any bar is missed):
    - ≥ 90% recall of planted single shifts in genes that align (homopolymers excluded);
    - ≥ 80% recall of frameshifts RefSeq already annotates;
    - 0 flags among genomes with no known shift;
    - ≤ 1 lineage flag among shuffled pools (empirical null);
    - two shifts planted at one codon flag in ≥ 90% of trials;
    - **real positives:** *A. vastus* must flag; also run a CTGY yeast (*Kazachstania sinensis*) and
      *Perkinsus* if they can be fetched, and report hit or miss;
    - **harder null:** ≥ 100 named genomes fixed in advance, plus copies with homopolymer and random
      indels planted at realistic sequencing-error rates;
    - **coverage:** per phylum, the share of genes that align at all (the 90% bar is conditional on
      alignment, so it cannot show where the scan is blind).
12. **Planting:** 1, 2 and 5 shifts per genome, at random positions or at one chosen codon, +1 and −1,
    in genomes from 5 phyla.
13. **Fungal reads** built their own assemblies, so read support is necessary but not sufficient;
    independence needs a second run or a second assembly.
14. **RNA-seq:** ≥ 10 reads across a site to call "kept" or "edited". First show each library is RNA,
    not DNA (coverage shaped like transcripts, not flat across the genome). Then the RNA step must get
    two known answers right: a site known to stay in the mRNA and a site known to be edited.
15. **Blinded check:** 2 decoys per candidate (≥ 1 genome with no known shift, ≥ 1 known single-site
    genome such as a bird), coded labels, miniprot, told only "find and describe frameshifts". It runs
    in an empty folder **outside the repo** holding only the coded genomes and reference proteins, with
    no CLAUDE.md and no key. It writes its verdict before unblinding.
- **Step 4 funnel addition:** drop a shifted gene when the genome has an intact second copy of it.
16. Local git; each frozen step is committed before the next runs. Nothing is pushed.
17. This folder. `../tils-tada-scan` stays paused; its `data/proteomes_pilot` is kept.
18. Write-up as `README.md`, `notes.md` and `data/`, laid out like the CGG folder.
19. Tools from Homebrew; exact versions go in `PREREG.md`.

## Amendment 1 (2026-11-21 20:11 EST, after the prior-art audit; before any scan code)

Decided by the user after reading `PRIOR_ART_AUDIT.md`:
- **Introns in the fungal set: alignment-based rule, no MFannot on all 2,695.**
  - **Calibrate first** on RefSeq fungal mitogenomes that have annotated introns. Report how often an
    annotated intron is inferred as an intron, how often an intron edge produces a false shift, and the
    recall of shifts planted near intron edges.
  - **Fix the bar** for the rule in `PREREG.md` before the 2,695 are scanned.
- **Masked coding:** report, per genome, the share of coding bases hidden by the 30-nt mask (intron-rich
  genes such as *cox1* may be mostly masked). This is part of what the scan could not see.
- **MFannot as a second detector,** run only on (a) every fungal genome that produces a candidate and (b) a
  stratified sample of about 100 others. Its "potential frameshift" flag counts as an independent
  detector; agreement is reported both ways.
- **Authors' MFannot output first.** Before any MFannot run, check the Sarhan et al. deposit for it.
  Checked 20:05–20:11:
  - Figshare 10.6084/m9.figshare.28750034 (version 1) holds only `Assembled_fungal_mitochondrial_genomes.fasta`,
    `fungal_mitochondrial_genomes_phylogeny.tre` and `All_fungal_mitochondrial_genomes_accessions.tsv`.
  - The Europe PMC copy of the paper has no supplementary material.
  - github.com/msabrysarhan/fungal_mtDNA has only `README.md` and `assembly_workflow.py`.
  - The TPA records carry no features.

  **Not public.** MFannot is run locally instead (amendment 2).
- **Named positives:**
  - ***Perkinsus*** (in the scan set) must flag, or be logged as a blind spot.
  - ***Magnusiomyces*** must appear in the bypass-insert report, or be logged as a blind spot.
- **One shift-codon convention** for real, planted and shuffled-null shifts alike.
- **Wolters 2023:** `PRIOR_ART.md` must state exactly what it did, with the passage quoted.

## Amendment 2 (2026-11-21 20:28 EST)

- **Standing rule** (also in `CLAUDE.md`): no contact with anyone outside the working session, by email
  or any other channel, and no drafting of messages for the user to send.
- **MFannot locally.** Wherever the Sarhan et al. authors' MFannot output would have been used, MFannot
  is run locally instead: on every fungal genome that produces a candidate, and on the stratified sample
  of about 100 others (amendment 1).

## Amendment 3 (2026-11-21 20:36 EST, user's conditions on the red-team additions; before any scan code)

1. **Stop codons as triggers: yes.** The expected rate for every trigger codon, stops above all, comes
   from **random shifts planted and run through the same pipeline and shift-codon convention**, never from
   in-frame codon counts. A random extra base sometimes creates a stop at the shift site itself, while the
   in-frame stop count is zero, so an in-frame null would call any two such sites significant. Stop
   triggers get **their own single-genome threshold** from the planted null. Both thresholds go in
   `PREREG.md`.
2. **"Possible code change" is defined by what the alignment does:**
   - frame unchanged after the stop → possible code change (record the residue relatives have at that
     column);
   - frame changed → frameshift candidate, even when the trigger is a stop;
   - paired detours around a stop → a third tally.

   These are reported, not pursued this round.
3. **Homopolymer tallies: exploratory.** "Possible poly(U) editing" is flagged only when all three hold:
   - the excess runs one way (DNA shorter than the protein needs);
   - it sits in coding-strand T runs, not A runs;
   - it survives a check of the sequencing platform, since nanopore-style homopolymer deletions would
     otherwise give the same picture.
4. **Held-out sets: yes.** They are labelled **unseen known positives, a test of the pipeline**, not
   replication of our own flags. Only already-deposited assembled sequences are used; nothing is
   assembled from raw reads for this.
5. **New funnel step.** Before any candidate is called new, a fresh agent searches the literature for that
   lineage and gene with every synonym now known: frameshift, bypass, jumping, hopping, editing,
   insertion, extra nucleotide, pseudogene, unusual code.
6. **Scope rule.** One confirmatory test: **the codon pile-up of clean shifts.** Bypass inserts, code
   changes and poly(U) editing are tallies in the report, each labelled exploratory. No new pipelines for
   them this round.

- **RNA known-answer test:** chicken *nad3*-174 is a documented "kept" site (Russell & Beckenbach 2008
  sequenced the mRNA).
- **Build order:** only what reaches the frozen scan, namely the planted null, the named positives, the
  intron calibration, `PREREG.md` and the scan itself. MFannot, the held-out sets, reads, RNA and the blind
  check wait until a flag needs them.
- **Stop rule:**
  - After the scan, stop and report which known lineages flagged and every flag that is not a known
    lineage.
  - If there are none, write a short negative with the measured sensitivity and **stop the project
    there**.

Extra leads for prior art: Gornik 2022 (*Perkinsus*), Andreu-Sánchez 2021 (bird and turtle *nad3*),
Russell & Beckenbach 2008 (turtles), Rosengarten 2008 and Haen 2014 (glass sponges), Szabóová 2018
(yeasts), starting from Al Arab et al. ("Accurate annotation of protein-coding genes in mitochondrial
genomes").

## Amendment 4 (2026-11-22 06:50 EST, user's decisions after the first control run missed bars; before the re-run)

The first full control run (commit 972ae02) missed B5, B2, B9a and B9c, and no alpha met the rule. The user decided:

1. **One codon rule.** A shift counts at either codon of its placement window (the old-frame codons holding
   the leftmost and the rightmost protein-equivalent placement). The planted null uses the identical rule.
   The width of window that may vote is capped, with the cap set from the planted data and written in
   PREREG.md. Shifts with wider windows count as detected (they are in n) but do not vote.
2. **B5 is scored on trials where both shifts were detected**, but the **unconditional number is the
   headline sensitivity** and goes at the top of the report. Also reported:
   - the unconditional single-site detection rate by group (animals, fungi, protists);
   - the reasons for misses (gene not aligned, masked, homopolymer, gene end, ...);
   - the chance of flagging a genome with 2, 3 and 5 planted sites.
3. **Alpha stays as the pre-written rule gave it: 1e-3 for sense codons, 3.2e-3 for stop codons.** The
   expected number of chance flags under exactly that rule is stated in PREREG.md. (Stop pairs at 0.86 vs
   1.00 is 36 vs 42 trials and does not drive the choice.)
4. **Frame changes between two hits of the same gene are detected**, with limits:
   - same reference protein;
   - contiguous protein coordinates;
   - the two hits abut, or overlap by no more than about 10 bases on the genome.
   Wider gaps stay introns or bypass-like inserts. Planted shifts and nulls go through the same path. B9b
   is re-checked, because this is where intron edges could start producing false shifts.
5. **The 2,695 fungal genomes stay in the main scan.** The false-positive controls pass, so the blind spot
   costs sensitivity only. B9a and B9c are recorded in PREREG.md as stated limits, with per-genome
   coverage and masked share, not as bars to chase. One fix is tried in the same re-run: drop reference
   proteins much longer than their gene's typical length (they carry intron-encoded sequence). If
   coverage does not improve, it is left.

Also:
- **NC_079576** ("Pholiota spumosa", likely a mislabelled ascomycete) is removed from the references and the
  controls, and logged.
- **One re-run only.** After it, the rules are frozen with whatever the numbers are, and any bar still
  missed is logged as a limit. No second round of tuning on the same controls.
- **Strong flags.** Several flags are expected by chance. In the scan report a flag is **strong** if it would
  be expected by chance less than once across the whole scan, or if two related genomes flag the same
  codon. Only strong flags that are not known lineages trigger the funnel. If there are none, the short
  negative is written and the project stops.
- **Scan start (user's condition).** The full scan (about 90 min) starts automatically after the freeze only
  if the re-run shows:
  - 0 flags in the null genomes;
  - B9b ≤ 2%;
  - every named positive still flagging.
  If any of these fails: stop before freezing and show the user. The unconditional sensitivity numbers go
  into PREREG.md before the scan starts. After the scan: report by the stop rule, and do not begin the
  funnel (reads, RNA, MFannot) without asking.
