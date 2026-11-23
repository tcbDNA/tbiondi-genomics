# Pre-registered rules, predictions and changes

**Provenance.**
- This file is a copy of `PREREG.md` from the project's working repository. It is reproduced verbatim except for this
  header and the commit tables at the end.
- The working repository is a **local git repository**. The commit hashes and times below come from it and were **not
  externally timestamped**: there is no public remote and no third-party timestamping service. They show the order in
  which entries were committed, as recorded by local git; they are not independent proof of timing.
- **History rewrite.** On 2026-11-22 the local history after commit d6e06a9 was rewritten to remove a 611 MB
  intermediate file committed by mistake. Messages, authors and dates were preserved, and no hash cited in this file
  changed. The old-to-new hash table is in the log entry "Local history rewritten in this session".
- **Paths.** The log refers to files by their paths in the working repository. This release renamed several of
  them; `protocol/corrections.md` has the old-to-new table. The log itself is not edited.
- Times are US Eastern (EST), 2026-11-21 to 2026-11-22. Where an entry's typed time differed from its commit time, the
  log says so (entry "correction of timestamps").

---

# Pre-registered rules, predictions and changes

Each section is committed to local git before the data it governs are scanned. Later changes are added
as dated log entries with the reason, never edited in place. Times are US Eastern.

## 0. Scope and fixed decisions (written 2026-11-21 20:36, before any scan code; committed with DESIGN.md amendment 3)

**Question.** Does any mitochondrial lineage shift reading frame at a recurring codon, as glass-sponge
mitochondria do at UGG?

**Data (main scan).**
- RefSeq mitochondrion release 238: 18,220 records, of which 17,648 are not land plants (Embryophyta).
- The 2,695 fungal TPA records of Sarhan et al. 2025 (BK072095–BK074789).
- Each genome is translated with its own mitochondrial genetic code: the record's CDS `transl_table`, or
  NCBI taxonomy's mitochondrial code where the record has no CDS.
- **Land plants:** 572 records, scanned afterwards as a separate stratum with its own null, exploratory.
- **Held-out set:** GenBank mitogenomes of species not in the scan. These are unseen known positives used
  as a pipeline test, not replication. Deposited assemblies only; no assembly from reads.

**The single confirmatory test: the codon pile-up of clean shifts** (amendment 3, item 6).
- **Clean shift:** a single frame change inside one DIAMOND frameshift alignment, passing the clean-shift
  filters fixed in section 1.
- **Trigger codon:** the old-frame codon assigned to the shift by one shift-codon convention (section 1),
  applied identically to real, planted and null shifts.
- **Single-genome flag:** ≥ 2 clean shifts of the same direction at the same trigger codon, in different
  genes. Their count must exceed the expectation **from planted random shifts run through the same
  pipeline and convention** (never from in-frame codon counts). The threshold is in section 1, with a
  separate threshold for stop-codon triggers.
- **Lineage flag:** NCBI family and order pools (nearest named clade above genus where a rank is missing),
  corrected for the two ranks. A shift at the same orthologous site counts once per pool. Expectations come
  from the members' planted nulls. Thresholds are in section 1.

**Exploratory tallies only** (reported, not tested, no new pipelines):
- bypass-like inserts;
- possible code change, defined by the alignment: frame unchanged after a stop = possible code change,
  with relatives' residue recorded; frame changed = frameshift candidate; paired detours around a stop =
  third tally;
- homopolymer-killed shifts and "possible poly(U) editing". This is flagged only when the excess is
  one-way (DNA shorter than the protein needs), in coding-strand T runs rather than A runs, and survives a
  check of the sequencing platform;
- the 4-base motif route, on pools with ≥ 8 distinct sites.

**Named positives.**
- *Aphrocallistes vastus* (NC_010769) must flag.
- *Perkinsus marinus* (NC_068742) and *P. olseni* (NC_068743) must flag, or be logged as a blind spot.
- The CTGY yeasts (*Huiozyma sinensis* NC_036748, *Monosporozyma servazzii* NC_004918) are run and
  reported, hit or miss.
- *Magnusiomyces capitatus* (NC_024095) must appear in the bypass-insert tally, or be logged as a blind
  spot.

**Funnel before anything is called new** (only once a flag exists):
1. not a homopolymer, intron edge or incomplete stop;
2. no intact second copy of the gene in the genome;
3. raw DNA reads carry it;
4. a second assembly carries it;
5. RNA keeps the base (shift) or lacks it (editing);
6. **a fresh agent searches the literature** for that lineage and gene with every known synonym
   (frameshift, bypass, jumping, hopping, editing, insertion, extra nucleotide, pseudogene, unusual code);
7. the blinded check (miniprot, coded decoys, isolated folder).

RNA known answers: kept = chicken *nad3*-174 (Russell & Beckenbach 2008), *Perkinsus*; edited =
*Physarum* (*Badhamia*) insertion sites.

**Build order and stop rule.**
- Built before the scan: the planted null, the named positives, the intron calibration, this file and the
  scan. MFannot, the held-out sets, reads, RNA and the blind check wait until a flag needs them.
- After the scan: stop and report which known lineages flagged, and every flag that is not a known lineage.
- If there are none: a short negative with the measured sensitivity, and the project stops.

**Tools.** DIAMOND 2.2.9 (Homebrew), miniprot 0.18-r281 (Homebrew; blinded check only), Python 3.12 via
uv (pyproject.toml / uv.lock), Biopython 1.86.

## 1. Calibrated rules and thresholds (frozen 2026-11-22 07:50, after the control re-run, before the scan)

These rules govern the scan of the 20,343 genomes: 17,648 RefSeq non-land-plant plus 2,695 fungal TPA.
The code is in `src/` at the freeze commit. The control results are in `results/controls/BARS.md` and
`thresholds.json`.

**1.1 Alignment.**
- DIAMOND 2.2.9: `blastx --sensitive -F 15 --range-culling -k 10 -e 1e-5 --query-gencode <genome's code>`.
- References: `data/work/ref_final`, 53,641 proteins. These are prep.py's family representatives (one RefSeq
  genome per family pool, the one with the most proteins, then the longest), minus NC_079576.
- The length filter was tried and not kept (PREREG log, length-filter pre-check result).
- Hits to the query's own family are dropped. Hits are then culled greedily by bitscore: a hit is kept if it
  overlaps kept hits on the same strand by less than 50% of its length.

**1.2 Shift events.**
- (a) **Inside one alignment.** Frameshift operations in the BTOP. Adjacent operations are merged, and the
  net offset is folded: +1/−2 → +1, −1/+2 → −1, multiples of 3 dropped.
- (b) **Between two hits** (split events), when all of the following hold:
  - consecutive hits of the same reference protein on one strand;
  - protein gap d_s = b.sstart − a.send − 1 in [−3, +1];
  - genome gap d_q ≥ −10;
  - residual r = d_q − 3·d_s in {±1, ±2}.
  Split events vote only if they passed probation (1.9).

**1.3 Classes**, assigned in this order (first match wins):
1. paired: an opposite-direction shift within 150 nt on the same reference protein. A single-hit shift is
   paired only with single-hit shifts.
2. masked_insert_edge: within 30 nt of an inferred insert longer than 30 nt.
3. gene_start: subject residue ≤ 3.
4. gene_end: within the last 10 residues.
5. flank_fail: fewer than 15 aligned columns on either side, or identity below 0.20 on either side.
6. homopolymer: a run of ≥ 5 identical bases for +1 shifts, or ≥ 4 for −1 shifts, reaching within 6 nt of
   the placement window.
7. clean: everything else.

**1.4 Inserts** (between hits of one gene, and query insertions ≥ 20 nt inside a hit):
- bypass-like: 20–100 nt;
- intron: > 100 nt.
These are tallies and masks only.

**1.5 Window rule.**
- Placement window: every placement within ±9 nt that gives the same protein over the local window
  (−24/+27 nt).
- R and L are its rightmost and leftmost placements; window = R − L.
- A clean shift votes for the old-frame codon holding R and the one holding L (one vote if they are the
  same codon), if window ≤ **cap = 4 nt**.
- The cap is the smallest width holding ≥ 95% of clean recovered planted shifts in the planted null.
- Wider windows count in n but do not vote.

**1.6 Planted null.**
- 620 planted genomes, stratified by code (seed 11). Three copies each, with one random insertion or
  deletion per core gene per copy (seed 101):
  - ≥ 60 nt inside the CDS, ≥ 90 nt apart, not next to a homopolymer;
  - inserted base drawn from the genome's composition.
- The copies go through the identical pipeline. Clean recovered shifts give vote counts by stratum:
  - strata: code × GC half (cut = median GC of the code's planted genomes, when ≥ 20 genomes);
  - backoff to code, then all, requiring ≥ 500 clean plants;
  - q(c) = (votes for c + 0.5) / (N + 32).
- Stops get their own q from the same counts, never from in-frame codon counts.

**1.7 Single-genome test (confirmatory).**
- Core genes only: cox1–3, cob, atp6/8/9, nad1–6, nad4L.
- Per genome and direction: n = clean shifts, k = those voting for codon c, g = distinct genes among them.
- Flag if k ≥ 2, g ≥ 2 and P(Binomial(n, q(c)) ≥ k) ≤ α.
- **α_sense = 1e-3** for sense codons and **α_stop = 3.2e-3** for the genome code's stop codons. These are
  the values the pre-written rule gave on the first run under the window convention (amendment 4,
  item 3). One test; no correction across conventions.

**1.8 Pool test.**
- Family pools and order pools.
- Distinct sites: same gene, same direction, subject positions within **50 residues** (PREREG log 07:10).
  A site votes for c if at least half of its member shifts vote for c.
- Flag if k ≥ 3 sites in ≥ 2 genes and the Poisson-binomial tail over the direction's sites ≤ α_pool/4 =
  2.5e-6, with each site's q taken from its first member's stratum.
- The /4 is kept from section 0 although only one codon rule remains. This is conservative.
- Shuffled-pool null: labels permuted, seeds 1–3, family and order.

**1.9 Split-hit probation result: split-hit shifts do NOT vote.**
- The re-run found 0 split events in the planted null, the null genomes and the same-codon trials.
- So 0 planted shifts were recovered only through a split event, and 0 null sequences flag only with
  split events voting.
- Split-hit shifts are therefore an exploratory tally. They are excluded from n, from the votes and from
  the pools, and the null is built the same way.

**1.10 Strong flags and the stop rule.**
- A genome flag is **strong** if either holds:
  - E < 1, where E = the expected number of (genome, direction, codon) results with p at most the flag's p
    across the whole scan, computed from the scan's own clean-shift counts under the planted null;
  - another genome in the same family or order pool flags the same codon in the same direction.
- A pool flag is strong if p × 64 × (number of pool-direction tests) < 1.
- Only strong flags outside the known lineages trigger the funnel. Known lineages are the family and order
  names in `src/analyze.py` KNOWN_TAXA, plus classes Aves and Hexactinellida.
- The funnel (reads, RNA, MFannot) is not begun without asking the user.
- If there are no strong flags outside known lineages: a short negative with the measured sensitivity
  below, and the project stops.

**1.11 Exploratory tallies** (reported, never tested):
- bypass-like inserts;
- stops read in frame (possible code change), with the relatives' residue;
- detours around stops;
- homopolymer-killed shifts by direction and base;
- genomes with clean shifts in ≥ 3 different core genes, whatever the codons;
- split-hit shifts, if they did not pass probation.

**1.12 Measured sensitivity (control re-run, finished 07:50; `results/controls/BARS.md`).** Headline,
unconditional:

| planted sites at one codon | sense codon: genome flagged | stop codon: genome flagged |
|---|---|---|
| 2 | **0.317** (n = 240) | 0.342 (n = 120) |
| 3 | **0.483** (n = 120) | 0.575 (n = 120) |
| 5 | **0.717** (n = 120) | 0.767 (n = 120) |

When every planted site is counted (clean, voting for the planted codon), the genome flags in 1.000 of
sense trials (n = 76, 22, 12 for 2, 3, 5 sites) and in 0.976–1.000 of stop trials. This is B5, conditional.
So the losses are almost all in single-site detection:

| planted shift | in gene (same direction, anywhere) | at the site (±9 nt) | clean at the site | **counted at the right codon** |
|---|---|---|---|---|
| random positions, animals (n = 14,355) | 0.906 | 0.886 | 0.680 | **0.564** |
| random positions, fungi (n = 6,054) | 0.832 | 0.816 | 0.650 | **0.528** |
| random positions, protists (n = 4,251) | 0.903 | 0.876 | 0.655 | **0.553** |
| same-codon trial sites, all (n = 2,640) | 0.900 | 0.853 | 0.618 | **0.518** |

Reasons a random planted shift is not counted (share of all planted):
- homopolymer kill: 0.149;
- clean but voting for another codon (misplaced): 0.100;
- not inside any aligned hit: 0.092 (fungi 0.151);
- flank filter: 0.045;
- aligned but no shift written: 0.041;
- window wider than the cap: 0.015;
- other: < 0.01.

**1.13 Bars, gates and limits.**
- **Gates** (all pass; the scan starts on them, amendment 4):
  - null genomes plus error copies: 0 flags of 360;
  - B9b: 0.0018 (n = 17,742);
  - named positives all flag:
    - *A. vastus* +1 TGG, p = 4.0e-4;
    - *P. marinus* and *P. olseni* +1 AGG, p = 2.7e-6;
    - *H. sinensis* and *M. servazzii* +1 CTG, p = 8.5e-6;
    - *Magnusiomyces* 65 bypass-like inserts.
- **Bars passed:**
  - B1 0.901 (n = 17,796);
  - B5 conditional 1.000 (n = 76).
  - The pre-written alpha rule, applied to this run, gives the same alpha (1e-3 / 3.2e-3).
- **Bar missed, logged as a limit: B2.**
  - Annotated shifts detected at the site: nad3 0.756 (n = 619), other 0.788 (n = 66).
  - Detected anywhere in the gene: 0.989 and 0.894.
  - Counted at the right codon: 0.590 and 0.409.
  - Cause: placement error in low-similarity stretches (PREREG log, B2 correction). For the 155 shifts
    detected in the gene but not at the site, the offset has median 75 nt, 95th percentile 120 nt, and
    153 of 155 are within 150 nt.
- **Stated limits (fungal introns):**
  - B9a 0.619 (n = 8,871);
  - B9c 0.279 (n = 2,101);
  - near-edge planted shifts inside an aligned hit 0.438;
  - per-genome core-CDS coverage, median 0.795 and 10th percentile 0.647 (n = 663);
  - masked coding share, median 0.030, 90th percentile 0.061, max 0.288.
- **Pool tolerance check.** Bird and turtle nad3 clean shifts in the annotated set:
  - old 2% rule: 13 of 108 family or order pools have ≥ 3 distinct sites (max 6, Charadriiformes), so they
    could flag by themselves;
  - 50-residue rule: 0 pools have ≥ 3 (max 2).

**1.14 Expected chance flags under exactly this rule: about 0.**
- Estimated from the 120 null genomes, weighted by phylum: 233 of their 240 genome-directions have 0 clean
  core shifts and 7 have 1. None can reach k ≥ 2, so the estimate is 0.00 across 20,343 genomes.
- Caveat: these null genomes are well-annotated records (≥ 8 CDS, > 5 kb). Lower-quality assemblies in the
  scan will have more clean shifts.
- The number is recomputed exactly from the scan's own counts and reported beside the observed flags.
- (My earlier figure of "about 6 at alpha = 1e-3" to the user was not computed this way and should not be
  used.)

**1.15 Predictions, written before the scan.**
- P1: the five named positives flag as in the controls; same code and same genomes.
- P2: at least one of the other three glass sponges in the scan (*Oopsacas minuta*, *Vazella pourtalesii*,
  *Caulophacus iocasicus*) flags at +1 TGG.
- P3: some other CTG-group Saccharomycetaceae (*Kazachstania*, *Nakaseomyces* and relatives, RefSeq or
  fungal TPA) flag at +1 CTG.
- P4: no bird or turtle genome flags from nad3, and no bird or turtle pool flags.
- P5: the number of single-genome flags outside known lineages is small: of the order of the recomputed
  chance expectation plus a few assembly artifacts.
- P6: *Magnusiomyces* appears in the bypass-like tally. Other Dipodascaceae and Trichomonascaceae in the
  fungal set may too.

**1.16 Tools.** DIAMOND 2.2.9 (Homebrew); Python 3.12.2 via uv (pyproject.toml / uv.lock); Biopython 1.88.
miniprot 0.18-r281 is for the blinded check only.

## Log
- 2026-11-21 20:36: section 0 written.
- 2026-11-21 21:29: first full control run (commit 972ae02) exited 1. Missed: B5 (0.28), B2 (nad3 0.756,
  other 0.788), B9a (0.62), B9c (0.28); no alpha met the rule.
- 2026-11-22 06:50: user decisions recorded in DESIGN.md amendment 4. Changes to section 0:
  - the two placement conventions (R, L) are replaced by one window rule, so the single-genome test is
    no longer split across two conventions;
  - a clean shift may also be a frame change between two hits of the same reference protein (limits in
    section 1), not only a change inside one DIAMOND alignment;
  - B9a and B9c become stated limits, not bars;
  - "strong flag" defined; only strong flags outside known lineages enter the funnel.
- 2026-11-22 06:50: correction to section 0 "Tools": the Biopython installed in the uv environment is
  1.88, not 1.86.
- 2026-11-22 06:50: NC_079576 ("Pholiota spumosa", Strophariaceae) is removed from the reference proteins
  and from the intron control set. All 13 of its proteins are 95–98% identical to *Colletotrichum*
  (Glomerellaceae, Ascomycota), so the record is almost certainly a mislabelled ascomycete. It stays in
  the scan set under its database taxonomy. Its Strophariaceae pool assignment is suspect, and any flag
  in it will be read with that in mind.
- 2026-11-22 06:58: correction to the NC_079576 entry above. Re-checked with DIAMOND blastp against the reference set
  (itself excluded): NC_079576 has 15 proteins.
  - 14 have *Colletotrichum graminicola* (NW_007361658) as their best hit, at 71–100% identity.
  - 10 of the 14 core-gene proteins are ≥ 95% identical.
  - The best hit for nad3 is *Ophiocordyceps sinensis*, at 84%.
  - The conclusion stands; the "13" and "95–98%" in the entry above were wrong.
- 2026-11-22 07:40: **correction to my B2 diagnosis** (the report that led to amendment 4, item 4). I said most missed
  annotated nad3 shifts were cases where DIAMOND splits the gene into two hits in different frames. That was wrong.
  Re-checked on the first run's own hits, before the re-run:
  - Missed nad3 sites: 146 of 151 have a same-direction shift written 10–150 nt from the annotated site, mostly
    downstream in the variable stretch after nad3-174, often failing the flank filter. The other 5 have only an
    opposite-direction shift nearby.
  - Other missed annotated sites: 7 of 14 are the same kind; 6 have no event within 150 nt; 1 has only an
    opposite-direction shift.
  - No pair of hits meets the split limits anywhere in the 632 annotated genomes (0 split events).
  - B2's misses are therefore misplacement in low-similarity stretches, not split hits.
  - Split detection is kept as decided: it is cheap, it goes through the same path for real and planted shifts,
    and B9b is re-checked. It is not expected to move B2.
  - Not added: a fix for misplacement. That would be a new design change; amendment 4 allows one re-run with no
    further tuning, so B2 will be logged as a limit if it misses again.
- 2026-11-22 07:40: rule fixed before the intron pre-check (amendment 4, item 5: "if coverage does not improve,
  leave it"):
  - The length filter (core-gene reference proteins > 1.5 x their gene's median length within the same kingdom
    group: 189 proteins) is kept only if the intron stage, run with it, raises near-edge coverage (planted
    near-edge shifts inside an aligned hit) or B9a by >= 0.05 absolute over the first run (0.438 and 0.619).
  - Otherwise the references revert to the unfiltered set (minus NC_079576), and every stage, including the
    intron stage, is run on that.
- 2026-11-22 08:02: **user changes after the B2 correction**, logged before any re-run numbers were seen. The
  intron pre-check was running; its output, which concerns B9 only, had not been read.
  1. **Split-hit shifts on probation.**
     - They vote in the codon test only if the re-run shows both of the following:
       (a) they recover planted shifts that the single-hit path misses: at least one planted shift (planted
           null or same-codon trials) recovered as a clean shift only through a split event, with no clean
           single-hit event within 9 nt;
       (b) they add no null flags: no null-genome sequence (original or error copy) flags with split events
           voting that does not flag without them.
     - Otherwise they are an exploratory tally that neither votes nor flags. They are then excluded from n,
       from the votes and from the pools, and the planted null is built the same way.
     - Either way, single-hit shifts are classified without reference to split events: the opposite-direction
       pairing check for a single-hit shift looks only at single-hit shifts. The single-hit test is therefore
       identical whether or not split events exist.
     - The decision is computed by `bars` and recorded in thresholds.json.
  2. **Misplacement.** A shift placed 10–150 nt off votes for the wrong codon. So:
     - **B2 is reported two ways:**
       - annotated shifts detected anywhere in the gene: a same-direction shift event, any class, inside the
         annotated CDS holding the site;
       - detected within the voting window: a clean same-direction shift within 9 nt whose votes include the
         site's own codons. "Own codons" are what trigger() gives for a shift placed exactly at the annotated
         position.
       The bar itself is unchanged (detection within 9 nt, >= 0.80). The two extra rates are reported beside it.
     - **The same split for planted shifts** (planted null and same-codon trials): in gene / at site (±9 nt) /
       clean at site / counted at the right codon. The right codon is the trial codon, or for random plants
       the codons trigger() gives at the true planted position.
     - **Pool "same site" tolerance.** It was relative subject position within 2% (≈ 2 residues in nad3), far
       narrower than a 150-nt placement error. It becomes: same gene, same direction, subject positions within
       50 residues (150 nt / 3, the widest placement error examined). The re-run reports the observed
       placement offsets of annotated shifts (detected in gene but off site), to show 50 residues is at least
       as wide. It also reports how many distinct sites bird and turtle nad3 form per family and order pool
       under the old and the new tolerance.
     - **New exploratory tally:** genomes with >= 3 clean shifts in different core genes, whatever the
       codons. Reported, not tested.
  Nothing else changes. The conditional go-ahead for the scan stands.
- 2026-11-22 07:12: **correction of timestamps in this log.** I typed the times on several entries above instead of
  reading the clock, and they ran ahead. The commit times are the true ones:
  - "06:50" entries (amendment 4, Biopython, NC_079576): committed 06:48 (57c8b8a).
  - "06:58" NC_079576 correction: committed 06:48 (942d74d).
  - "07:40" B2 correction and the length-filter rule: committed 06:56 (e878133).
  - "08:02" split-probation and misplacement entry: committed 07:10 (43f9003).
  The order relative to the data is unchanged. All of these were committed before any re-run output was read.
  The intron pre-check started at about 06:58, after e878133. From here on, times are read from the clock.
- 2026-11-22 07:13: **length-filter pre-check result: not kept.** The intron stage run with the length-filtered references
  (finished 07:12) gave:
  - near-edge coverage 0.439 (first run 0.438);
  - B9a 0.637 (first run 0.619).
  Neither rose by 0.05, so by the rule logged before the pre-check, the references revert to prep.py's set
  minus NC_079576 (53,641 proteins), and every stage is re-run on that, the intron stage included. The
  pre-check's B9b (0.0019) and B9c (0.284) are not used. Output kept in results/controls/precheck/rerun_intron_precheck.log.
- 2026-11-22 07:50: **section 1 frozen.**
  - The control re-run finished 07:50 with exit status 1: B2 missed, logged as a limit. All three scan
    gates passed.
  - Split-hit shifts failed probation (0 split events) and are a tally.
  - Window cap 4 nt. Alpha 1e-3 / 3.2e-3.
  - Per amendment 4, the scan (main + fungal) starts automatically after this commit.
- 2026-11-22 07:56: **added after the freeze (user request), while the scan was running and before any scan output was
  read. These are descriptive context for the report only; they do not change which flags are strong.**
  - For every flag:
    - the sites at the flagged codon (gene, position) and whether they are in different genes;
    - the genome's total clean shifts and its +1/−1 mix;
    - sequencing platform and assembly method where the record gives them:
      - RefSeq: COMMENT / Assembly-Data fields and the GenBank source accession;
      - fungal TPA: the SRA run named in the record, with platform from SRA run info, plus the Sarhan et al.
        2025 assembly pipeline;
    - whether any related genome (same family or order pool) has a clean shift voting for the same codon in
      the same direction, flagged or not;
    - whether a second assembly of the species exists: another record of the same taxon in the scan set, or
      a GenBank mitochondrial record of that taxon at least 80% of the flagged genome's length, other than
      the RefSeq record's own source entry.
  - **Prediction scoring**, fixed now before results. P2, P3 and P5 were worded loosely in section 1.15:
    - P1 pass: all five named positives flag (count x/5).
    - P2 pass: at least 1 of *Oopsacas minuta*, *Vazella pourtalesii* and *Caulophacus iocasicus* flags at
      +1 TGG (count x/3).
    - P3 pass: at least 1 genome of *Kazachstania*, *Nakaseomyces*, *Huiozyma* or *Monosporozyma*, other
      than the two named positives, flags at +1 CTG (count flagged / scanned).
    - P4 pass: no Aves or Testudines genome flag and no Aves or Testudines pool flag involves a nad3 shift.
      All Aves/Testudines flags are counted and listed either way.
    - P5 pass: single-genome flags outside known lineages ≤ (recomputed chance expectation + 5). "A few" is
      fixed here as at most 5.
    - P6 pass: *Magnusiomyces capitatus* has ≥ 1 bypass-like insert in the scan. Other Dipodascaceae and
      Trichomonascaceae with bypass-like inserts are counted, as description only.
- 2026-11-22 08:12: **added after the freeze (user request), descriptive only, written before any full-scan output was
  read.** The only scan data seen at this point is the code-1 peek logged in PRIOR_ART (08:09).
  - The exploratory many-shift tally (genomes with clean shifts in ≥ 3 different core genes) is listed in full,
    sorted by number of clean shifts. Each genome gets:
    - the share of its clean shifts in the majority direction;
    - platform and assembly method, as for flags;
    - related genomes showing the same pattern: genomes in the same family or order pool that are also in the
      tally with the same majority direction;
    - whether it is a known editing lineage, with the citation.
  - Any genome that is not a known editor and has ≥ 80% of its clean shifts in one direction is marked
    **"possible editing or noisy assembly"**.
  - Known insertion/deletion editing lineages, fixed now from PRIOR_ART section 5 and the 08:09 entry:
    - order Physariida (myxomycetes; *Badhamia*/*Physarum*, *Didymium*): Mahendran et al. 1991; Gott et al.
      2005; Bundschuh et al. 2011;
    - Heterolobosea records sequenced in Yang et al. 2017 (records citing PMID 28453770: NC_034794,
      NC_034798): "a novel RNA editing system in deep-branching heteroloboseids". Other Heterolobosea are
      not labelled;
    - class Kinetoplastea: Benne et al. 1986;
    - class Diplonemea: Moreira et al. 2016;
    - Calcaronean sponges (families Sycettidae, Leucosoleniidae, Grantiidae, Heteropiidae, Amphoriscidae,
      Jenkinidae, Lelapiidae): Lavrov et al. 2016.
      The scan's Calcarea are Clathrinidae (Calcinea), which are not labelled; only tRNA editing is reported
      for them.
- 2026-11-22 08:54: **added after the freeze (user request), descriptive only.** The report's table of known-lineage flags gets a
  per-site column: is the site described in a paper (cited, with our reading depth) or only annotated in a
  database record?
  - Sources are only what is already in PRIOR_ART.md (paper, by lineage and gene) and the local RefSeq
    annotation table (results/audit/known_junctions_refseq238.tsv: an annotated shift in the same record, same
    gene, within 150 nt).
  - No new literature search is run for this column. Sites not covered say "not in PRIOR_ART (not searched)".
  - Context at the time of writing: the vertebrate-batch peek (08:47) showed four known-lineage flags (two grebes,
    two emydid turtles, +1 AGT).
  - **Not started, on the user's instruction:** a possible later cheap exploratory look at the grebe nad1 site (a
    fresh-agent literature search, and a comparison across all deposited grebe and flamingo mitogenomes, no read
    downloads). It will only run if the user asks after the scan report.
- 2026-11-22 09:14: **four additions (user request), descriptive or conservative, logged before the analysis runs.** The scan
  was still in the code-5 batch. Peeks seen so far: codes 1–4.
  1. **Completeness.** Per batch (set × code), the report gives:
     - genomes expected (scan-set rows in data/work/genomes.tsv);
     - scored (rows in the scan's genomes.tsv);
     - failed (expected but not scored);
     - scored with no core gene aligned.
     A genome with no result row is "not scored", never 0 shifts. Every negative statement gives the number of
     genomes scored. (Cause: my code-4 peek read the still-empty code-5 hits file as 0 shifts for two glass
     sponges. That peek line was wrong and was withdrawn.)
  2. **Unit for the strong rule (conservative).** Records that share a BioSample or an Assembly accession
     (RefSeq DBLINK; for the fungal TPA records, the BioSample in the Sarhan accession table), or that have an
     identical organism name, are one unit. Example: the three *V. dahliae* VdLs.17 scaffolds NW_009276968,
     NW_009276969 and NW_009276970, BioSample SAMN02953716, GCF_000150675.1. Records of one unit can never
     satisfy "two related genomes flag the same codon". E is left as computed over records, which is
     conservative.
  3. **Real-data sensitivity for the known lineages**:
     - lineages: CTG yeasts (*Kazachstania*, *Nakaseomyces*, *Huiozyma*, *Monosporozyma*), glass sponges
       (Hexactinellida), *Perkinsus* (Perkinsidae), birds (Aves) and turtles (Testudines);
     - for every scanned genome of these lineages:
       - its documented sites: RefSeq record annotations, plus paper sites in PRIOR_ART (gene-level where no
         position is given: *K. sinensis* cox3/atp6/cob; *K. servazzii* atp6/cox2/cox3/cob; *A. vastus* cox3
         and nad6);
       - per site: detected at the site (±9 nt), clean, counted at its own codon, or only elsewhere in the gene,
         with the reason for each miss;
       - flagged or not.
     - A genome with no documented sites is "no documented sites", not a miss.
     - Shown next to the planted sensitivity. For each failed prediction, the report says whether this
       explains the failure.
  4. **Each flag's own chance of flagging.** Monte Carlo from the planted null: 20,000 draws of the flag's n
     clean shifts (its direction), with replacement, from the clean planted shifts of its stratum (with their
     genes and votes), run through the identical test. The share that flag at any codon is reported, beside
     the analytic upper bound.
  - *Verticillium* context: the independent assemblies are named instead of the same-assembly scaffolds.
    *V. dahliae* NC_008248 (Pantou et al. 2006, 27,184 bp) and *V. nonalfalfae* NC_029238 are both in the scan
    set, and their clean-shift counts are reported, with the other *Verticillium* records in the scan
    (NW_003315039 *V. alfalfae*; BK074514–BK074519). Comparing bases at the VdLs.17 shift positions waits
    until after the report, with the grebe and tilapia checks, and only if asked.
- 2026-11-22 09:21: **two refinements to the descriptive context, before the analysis runs.**
  - "Second assembly":
    - the length threshold is 80% of the smaller of the flagged record's length and the median length of its
      family in the scan set, because the VdLs.17 scaffold (42 kb) is longer than the whole *V. dahliae*
      mitogenome (27 kb);
    - the GenBank search is at species level (strains included);
    - records already in the scan, and their own GenBank sources, are listed separately, not counted again.
  - Documented sites in non-core genes (e.g. *rps3*) are shown as "non-core gene (not tested)" and are left out of
    the sensitivity counts. Yeast gene names such as "Cox3p" are normalised to the core names.
  - Tested on the code-3 and code-4 batches rebuilt from their saved hits (not the final analysis). In that partial
    set the VdLs.17 flag had E = 0.21. E for the report comes from the full scan.
- 2026-11-22 09:24: **user instructions before the analysis runs** (the scan was in the code-5 batch).
  1. **Pre-approved, conditional: the second-assembly step for *V. dahliae* VdLs.17, and nothing else from the
     funnel.** It runs only if VdLs.17 (NW_009276968/69/70, one unit) comes out strong.
     - Positions: every clean shift on all three VdLs.17 scaffolds.
     - Comparison: for each position, a 121-nt window (±60 nt) is aligned locally (Biopython PairwiseAligner,
       local, match 2 / mismatch −3 / gap open −5 / gap extend −2) to *V. dahliae* NC_008248 (Pantou et al.
       2006) and to *V. nonalfalfae* NC_029238. Both are local RefSeq records; no reads, no downloads.
     - A window is comparable if the local alignment covers ≥ 100 VdLs.17 bases at ≥ 85% identity.
     - At a −1 shift (VdLs.17 lacks a base), the other assembly "carries the base" if the alignment has a gap in
       VdLs.17 within 3 nt of the shift. At a +1 shift, it "carries the correct form" if the alignment has a gap
       in the other assembly there.
     - A position counts as carried when NC_008248 carries it, or, where NC_008248 is not comparable, when
       NC_029238 does.
     - Kill rule: carried at > 50% of comparable positions → reported as "strong, killed at the second-assembly
       step: assembly error". Otherwise the flag is reported with the counts and the funnel waits for the user.
  2. **Limits for the write-up**, from the real-data table:
     - Documented real sites are counted at their codon less often than planted ones (16 of 48 in the code-3/4
       test, against 0.553 planted). What a negative rules out is stated with the real per-site rate and its
       sample size, and with the genome-level flag chance that rate implies for 2, 3 and 5 sites.
     - The pair filter removes real sites in lineages that mix +1 and +2 shifts close together (*Perkinsus*: 14 of
       38 in the test).
     - The strong rule's chance figures assume errors land at random codons. Assembly errors cluster, so a noisy
       genome can be strong without being real.
  3. **Exploratory safety net (post-freeze, never a flag):** paired breaks per genome, by the codon(s) of the
     first break of each pair (window rule). Genomes where ≥ 2 pairs share a first-break codon are listed.
  4. **Units corrected:** grouping by organism name merges records only when the name carries the same strain,
     i.e. a token after the binomial that contains a digit (e.g. "VdLs.17", "CBS 138", "sp. BB2").
     - Subspecies trinomials and bare binomials are not merged by name.
     - BioSample and Assembly links still merge.
     - So NC_008248 and VdLs.17 are separate units.
- 2026-11-22 09:34: **two adjustments to descriptive outputs (user), before the analysis runs.**
  1. **Real per-site rate.**
     - Reported per lineage. Lines: CTG yeasts, glass sponges, *Perkinsus*, birds excluding nad3-174, turtles
       excluding nad3-174.
     - **nad3-174 gets its own line**: bird or turtle nad3 sites at annotated CDS position 160–190. It is one
       inherited site seen many times.
     - **The headline rate excludes nad3-174.** It is the mean of the per-lineage rates over lineages with ≥ 1
       documented core-gene site, each lineage weighted equally, because sites within a lineage are not
       independent. The pooled rate is shown beside it.
     - Interval: a bootstrap over lineages (2,000 resamples, seed 1) if ≥ 5 lineages have sites; otherwise the
       range of the per-lineage rates.
     - The 2/3/5-site implications use the headline rate and its interval bounds.
  2. **Paired-break tally:** each listed genome also gets:
     - its total paired breaks (core genes) and pairs;
     - the number of pairs at the recurring codon;
     - the mean identity of its kept alignments to the references: identical over aligned amino-acid columns,
       from the saved DIAMOND hits.
     A recurrence is read only against these, since divergent genomes produce detours by chance.
- 2026-11-22 09:41: **two descriptive additions (user), before the analysis runs.**
  1. **Removed breaks in known groups.**
     - For the three glass sponges besides *A. vastus* (*Oopsacas minuta* NC_027419, *Vazella pourtalesii*
       NC_028054, *Caulophacus iocasicus* NC_072601), for *A. vastus* itself, and for *Polyrhachis*, every break
       DIAMOND wrote is listed: gene, position, codons (L/R), the codon after it, and the class (filter) that
       removed it, or "clean".
     - TGG followed by GGA (the *A. vastus* signature) is marked.
     - Documented sites are checked in PRIOR_ART and in Haen et al. 2014 (as far as its text can be read; the
       reading depth is stated).
     - The P2 line says whether the failure is due to filters (sites seen and removed), absence (no sites
       documented and none seen), or unknown.
  2. **Out-of-sample sensitivity.**
     - Every known group in results/audit/controls_known_cases.tsv joins the real-data table, single-site groups
       included: oysters, venerid clams, *Polyrhachis*, *Arctica*, turtle nad4/nad4l, *Phaeodactylum*,
       *Dipterophagus*, the chytrids, and the other RefSeq single annotations.
     - The per-site rate is reported separately for **out-of-sample sites**: documented sites that played no part
       in setting the rules, i.e. not in a named-positive genome and not in the B2 set (known_junctions_refseq238).
       This is the cleanest estimate of real sensitivity; its sample size will be small and is stated.
  - **Not started, on the user's instruction:** running the frozen pipeline on held-out known positives with
    deposited assemblies (*Heteranomia* first). It runs only if asked after the report.
- 2026-11-22 10:03: **user decisions after the scan report** (results/scan/REPORT.md, commit 21ba959). **The confirmatory
  result stands as frozen; everything below is exploratory.**
  1. ***Grigorovia jiainica* (BK073075, fungal TPA), +1 AGG flag.** The user's statement, written before any
     check: the family rule was meant to cover the CTG mechanism; this flag is at AGG, a codon not described in
     yeasts, so the user is sending it to the funnel **after seeing the result**.
     Steps, in order, stopping at the first failure:
     a. **Sites.** For each shift: ±30 nt, the codon after it, whether a CTG sits in the placement window, and
        whether the shifted frame runs to the gene's end. Stop if fewer than 2 shifts are really at AGG.
     b. **Residue.** What do relatives have at the position facing AGG+N (skipping the A leaves GGN, glycine)?
        Continue only if relatives have glycine at ≥ 2 sites.
     c. **In-frame AGG.** How many AGG codons are read in frame elsewhere in its aligned genes?
     d. **tRNAs.** Does its mitogenome keep a tRNA that can read AGG (anticodon CCU or UCU)? Compared with its
        closest relatives in the scan.
     e. **Relatives.** Do its closest relatives carry the extra base at the same positions?
     f. **Reads.** Only if a and b pass: the SRA run behind the assembly (SRR16974510), asking first if it is over
        20 GB. Stop if the reads lack the extra base.
     g. **Literature.** A fresh agent searches with all synonyms, for this species under any former name.
     The report comes after g, or at the first stop. No RNA without asking.
  2. **Then the held-out run**, *Heteranomia* first, with the prediction written before running; deposited
     assemblies only.
  3. **Then the grebe nad1 look**, as specified earlier (fresh-agent literature search; comparison across all
     deposited grebe and flamingo mitogenomes; no reads).
  4. **Tilapia last**, skipped if time runs short.
- 2026-11-22 10:31: ***Grigorovia* step f read rule (user), logged before any read was counted** (download still running).
  - Per site, from the assembly's own run SRR16974510: exact matches of 20 nt + extra base + 20 nt versus
    20 + 20, both strands.
  - **Pass:** ≥ 10 spanning reads and ≥ 90% with the extra base.
  - **Mix:** 10–90%, reported as a mix, not a pass.
  - **Fail:** < 10% with the extra base.
  - A site with fewer than 10 spanning reads is reported as "too few reads".
  - Also owed and added (descriptive):
    - in-frame AGG per 1,000 aligned codons in the six closest relatives;
    - completeness of the *Grigorovia* mitogenome (circular, core genes, size against relatives) and a second
      tRNA finder (tRNAscan-SE, organellar mode) for the missing tRNA-Pro and tRNA-Cys;
    - a search for any other data on this species or genus.
- 2026-11-22 10:43: ***Grigorovia* step f result: PASS, 4 of 4 sites** (SRR16974510, 5,193,602 reads; results/explore/grigorovia/
  step_f_reads.md). Reads with / without the extra base:
  - cox1@24710: 1,440 / 0;
  - cob@26495: 1,718 / 0;
  - cob@26476: 1,733 / 0;
  - cob@26265: 1,785 / 0.
  The extra bases are in the sequenced DNA, not assembly errors. Whether they are translated by frameshifting, or
  removed by RNA editing, is not tested (no RNA exists).
- 2026-11-22 10:43: **blinded check (user), set up before running.**
  - Three genomes with coded labels, assigned in random order (seed 20261122): *Grigorovia jiainica* BK073075.1;
    decoy 1, the closest relative with no clean shifts, *Kazachstania bromeliacearum* BK073182.1; decoy 2, a CTG
    yeast, *Monosporozyma solicola* BK073486.1.
  - Reference proteins: the core mitochondrial proteins (cox1, cox2, cox3, cob, atp6, atp8, atp9) of
    *Saccharomyces cerevisiae* S288C (NC_001224) and *Saccharomycodes ludwigii* (NC_037703), with coded headers.
  - Aligner: miniprot 0.18 (Homebrew). Genetic code 3 is stated for all three.
  - The folder is outside the repo and holds only the coded genomes, the proteins and a task file. No CLAUDE.md, no
    key; the key is in results/explore/grigorovia/blind_key.json.
  - A fresh agent is told only to find and describe frameshifts. For each one it gives the position, direction,
    nearby codons and the residue the reference protein has at the aligned position after the shift: the blind
    residue call.
  - It writes VERDICT.md before unblinding.
  - The literature search (step g, both names, all synonyms) runs as a separate fresh agent in parallel. The two
    agents share nothing.
- 2026-11-22 11:04: **blinded check result (unblinded after VERDICT.md was written; key committed beforehand in 452342e).**
  - genome_1 = *K. bromeliacearum* (no-shift decoy): 0 frameshifts in 7 core genes.
  - genome_2 = ***Grigorovia jiainica***: 4 frameshifts, all +1, the same four sites as the scan (cox1 24711; cob
    26496, 26477, 26263, where the scan has 26266 for the same AGG/ACT site). Motif: AGG (×3) or TGG (×1) followed by C;
    skipping the first base leaves GGC.
    - **Blind residue call: glycine at all four** (cox1 G78; cob G137, G143, G214), conserved in both references and
      in genomes 1 and 3, which use GGA/GGT there.
    - GGC is used in frame in genome_2 only at these four sites, and AGG and TGG never.
    - The agent called it systematic. Counts were stable at miniprot -F 12 to 35.
  - genome_3 = *M. solicola* (CTG decoy): 5 frameshifts, all +1 at the C of CTG, read as TGT/TGC (Cys). This is the
    known CTG mechanism, recovered blind; the agent noted it as systematic but less tight.
  - Different aligner (miniprot) and different references; no access to our calls. The blinded check agrees with the
    scan on every site and direction.
  - Saved to results/explore/grigorovia/blind/.
- 2026-11-22 11:33: **held-out relative of a flagged lineage (original design, decision 7): *Grigorovia yakushimaensis* OR703946.**
  Written before the record was downloaded or opened. Its existence and the claim that its annotation reads cox1 and
  cob straight through came from the literature agent's report (Wolters 2023 data); that report has been read, the
  record has not.
  1. **Prediction:** if the AGG skip is a lineage trait, OR703946 has extra bases at AGG, i.e. either:
     - at least one at a position homologous to the four *G. jiainica* sites; or
     - at least 2 clean +1 shifts at AGG anywhere.
  2. **Exploratory prediction:** its mitochondrial tRNA-Gly is GCC if it skips at AGG, and UCC if it does not.
  3. **What each outcome means:**
     - shared sites = a lineage trait;
     - AGG shifts elsewhere = a trait with turnover;
     - no shifts = the sites are probably recent in *G. jiainica*. That is not a refutation.
  - **Method:** the frozen pipeline (src/scan.py process_genome, ref_final, code 3, own family Saccharomycetaceae
    excluded) run unchanged on the raw sequence, not the annotation. Homology to the *G. jiainica* sites by local
    alignment of ±60-nt windows, as in step e. tRNAs by ARAGORN and tRNAscan-SE -O.
  - **The record is read for:** any CDS notes, exceptions or pseudo flags; the submitter; the assembly method.
- 2026-11-22 11:33: **respiration (literature agent, follow-up).** No source reports growth of *G. jiainica* on non-fermentable
  carbon:
  - glycerol, ethanol, DL-lactate, succinate and citrate are all negative in Lee et al. 2008 (species description)
    and in Opulente et al. 2018 Table S10; acetate was never tested;
  - no petite test found;
  - *G. transvaalensis*: the collection strain does not grow on non-fermentable substrate (Fekete et al. 2007,
    petite-positive);
  - recorded as a major caveat: if *G. jiainica* does not respire, cox1/cob may be under relaxed selection and the
    extra bases could be disabling mutations rather than programmed shifts.
- 2026-11-22 11:36: **OR703946 (*G. yakushimaensis*) result**; the record had not been opened before the 11:31 entry.
  - **Record:** UNVERIFIED, linear, 49,874 bp, no CDS features (so no notes, exceptions or pseudo flags). Submitted by
    Wolters, LaBella, Opulente, Rokas and Hittinger (UW-Madison) in 2023. Assembly iWGS 1.1, NOVOPlasty 4.2 and
    plasmidSPAdes 3.9.0, from Illumina reads; BioSample SAMN08343276.
  - **Frozen pipeline, unchanged** (code 3, Saccharomycetaceae excluded):
    - core genes atp6, atp9, cob, cox1, cox2, cox3 aligned;
    - 0 clean shifts, 0 flags; the 4 breaks are all in non-core ORFs;
    - **0 clean +1 shifts at AGG**;
    - inferred intron in cob 21737–22763 (1,026 nt); two bypass-like inserts in cox1.
  - **The four homologous positions** (±60-nt windows, local alignment, checked by eye on the coding strand):
    - cox1 24711 and cob 26266: *G. yakushimaensis* lacks the extra base (120 nt aligned, 90% identity);
    - cob 26496 (AGG): lacks it ("GTATAT-AGG" against "GCGTATA-CGG"). The automatic ±3-nt test said "carries"
      because the gap is placed 4 nt away in a TA repeat; corrected by inspection;
    - cob 26477 (TGG): not comparable. *G. yakushimaensis* cob is interrupted here by the inferred intron; its
      downstream exon (from 22763) matches *G. jiainica* starting 4 nt after the TGG site;
    - so: no homologous site carries the extra base.
  - **Prediction 1 not met** (no AGG extra base at a homologous site, and < 2 clean +1 AGG shifts). By the logged
    meaning: "no shifts = the sites are probably recent in *G. jiainica*, not a refutation".
  - **Exploratory prediction 2 consistent:** its tRNA-Gly is UCC by both ARAGORN and tRNAscan-SE (score 56.1), as
    predicted for a genome that does not skip.
  - **Observation, post hoc, not interpreted further:** the *G. yakushimaensis* cob intron insertion point falls within
    a few bases of the *G. jiainica* TGG site.
- 2026-11-22 11:38: **tRNA oddity, post hoc and exploratory (user request).** results/explore/grigorovia/trna_gly_structure.md,
  fam_trna_gly.json.
  - ***Grigorovia* tRNA-Gly:**
    - anticodon GCC (tRNAscan-SE -O score 52.0; ARAGORN agrees);
    - anticodon loop 7 nt, UU-GCC-AA; **base 37 = A, the same as the relatives**;
    - the relatives' tRNA-Gly(UCC) loop is CU-UCC-AC, also 7 nt with A37;
    - no U37 and no extra loop base, so the CTG-yeast analogue (G37 proposed to pair with the skipped C) is not
      present in this form;
    - read support at anticodon position 34: GCC 1,329 reads, TCC 4, CCC 1, ACC 1.
  - **Glycine tRNAs:** each finder found exactly one tRNA-Gly in *Grigorovia* (GCC) and in each of the six relatives
    (UCC). *Grigorovia* still reads 38 GGA codons in frame; how GGA is decoded with only a G34 tRNA-Gly is unexplained
    (import or modification not tested).
  - **Across all 101 Saccharomycetaceae mitogenomes in the scan** (ARAGORN): tRNA-Gly(GCC) only in *G. jiainica*;
    the others have UCC (one has two copies).
    - Genomes with any +1 break at AGG: *G. jiainica* (3, all clean), *M. aerobia* (1 clean), and one non-clean break
      each in *T. delbrueckii*, *M. bulderi*, *S. cerevisiae* and *T. maleeae*.
    - The GCC anticodon and the clean AGG shifts co-occur only in *G. jiainica*: n = 1, possibly coincidence.
  - **Blind agent on *K. bromeliacearum*** (genome_1): no frameshifts called in the 7 core genes (VERDICT.md).
- 2026-11-22 11:50: **release draft** (user request: prepare locally, no remote, no push).
  - A curated copy is at ../release/mito-frameshift-scan; a snapshot is in release_draft/ (commit 6711345).
  - Folder name and headline, chosen after the OR703946, tRNA and respiration results: "A scan for codon-specific
    frameshifts in mitochondria". The scan's negative result leads, with its real-data sensitivity; *G. jiainica* is
    labelled exploratory.
  - Figures: the *G. jiainica* cob double site (static) and a ribosome animation labelled as a model.
  - Two fresh agents launched: a fact-check of every number and sentence against data/ and the cited sources, and a
    hostile referee on the main claim and its wording. Each fix from their reports is logged in
    release/.../data/review_log.md.
- 2026-11-22 11:54: **correction: SRA data for the *Grigorovia* subclade.**
  - What was wrong: the 10:36 search counted SRA runs by the genus taxid txid2917336 and reported "1 SRA run ... no
    data for the other Grigorovia species".
  - Cause: NCBI files *G. jiainica* under txid2917336 but *G. transvaalensis*, *G. humatica* and *G. yakushimaensis*
    under txid2862416.
  - Searching by species name gives 9 SRA entries, all DNA, none RNA:
    - *G. jiainica*: 1 (SRR16974510);
    - *G. yakushimaensis*: 2 WGS (SRR6475932, SRR16974690);
    - *G. transvaalensis*: 3 WGS (SRR6475934, SRR16974369, ERR9706995);
    - *G. humatica*: 1 WGS (SRR16974597) and 2 amplicon runs.
  - Corrected in results/explore/grigorovia/completeness_and_data.md and the release copy. No conclusion changes:
    there is still no RNA, and only *G. yakushimaensis* has a mitogenome assembly.
- 2026-11-22 11:57: **figure (user request): figures/ribosome_agg.gif**, from src/figures/ribosome_agg.py.
  - It draws from data/grigorovia/figure_site.json alone, the same file as the static figure.
  - It shows the *G. jiainica* cob site at 26496 as a labelled model: title "If the ribosome skips the A: one possible
    model", closing frame "Model, not observation ... untested".
  - Anticodons are those of the tRNAs the *G. jiainica* mitogenome encodes (tRNAscan-SE), drawn 3'->5' under each
    codon.
  - All 9 pairings were checked by the script: Watson-Crick, or G-U wobble at codon position 3 (Cys, Tyr, His).
  - Residues follow NCBI table 3. The *S. cerevisiae* Cob residues match at all 9 positions.
  - The earlier figures/ribosome.gif was replaced. The static-figure script moved to src/figures/site_figure.py.
- 2026-11-22 12:11: **corrections from the release fact-check** (two fact-check agent runs and one referee agent;
  release/.../data/review_log.md).
  1. **The 10:43 and 11:04 entries say the blinded check used "different references". That is wrong.**
     - Its references, *S. cerevisiae* NC_001224 and *Saccharomycodes ludwigii* NC_037703, are both in the scan's
       reference set.
     - The scan called the four *G. jiainica* sites against NC_037703 cox1 and cob proteins.
     - Only the aligner (miniprot) and the analysis were independent. The blinded task also asked for the reference
       residue after each shift and a "systematic or scattered" judgement, not only "find frameshifts".
  2. **The 11:36 entry says the record had not been opened "before the 11:31 entry".** The prediction entry is dated
     11:33 (commit ea426e3).
  3. **Process slips not logged at the time:**
     - an ARAGORN output-parsing error (10:12; re-parsed, nothing reported from it);
     - tRNAscan-SE waiting at an overwrite prompt (it ran from 11:24 until stopped at 11:35, then rerun);
     - a wrong DOI for Gornik et al. 2022 in the first README draft, fixed before review.
- 2026-11-22 12:11: **checks prompted by the referee (post hoc, exploratory;** results/explore/grigorovia/introns_and_congener.md).
  - **Intron positions.** In *S. cerevisiae* (NC_001224) coordinates, three of the four extra bases sit 3, 7 and 12 nt
    from intron insertion sites: cob bI2, cob bI1 and cox1 aI3. The fourth sits more than 100 nt from any intron.
    *G. yakushimaensis* has an intron at the TGG site.
  - **Congener at cox1.** *G. yakushimaensis* has AGT (Ser) where the skip-A reading of *G. jiainica* gives GGC (Gly).
    Placing the extra base at the G of AGG instead gives Ser and fits the congener. At cob 26266 the congener has the
    same two residues (Thr-Gly) as either placement.
  - **Consequences:**
    - the release text no longer presents the *G. jiainica* sites as favouring a translational skip;
    - the three readings (skip, RNA editing, disabling mutations) are given equal standing, with the respiration data,
      the glycine-tRNA problem, the intron positions and the congener's Ser stated in the opening;
    - the figure title was made neutral.
  - **Recommended by the referee and NOT done** (each needs the user's go-ahead):
    - the nuclear respiratory and translation-activator genes in GCA_030571715.1;
    - dN/dS on the *G. jiainica* branch, which needs mitogenomes assembled from the congeners' reads (the user's rule
      says no assembly from reads);
    - an intron map across the four *Grigorovia* species;
    - dropping the animation (the user asked for it, labelled as a model, so it was kept).
- 2026-11-22 12:11: note on the previous entry. The times given for the process slips ("10:12", "11:24 until stopped at 11:35") are
  approximate reconstructions; they were not recorded when the slips happened.

### Exploratory: settling *G. jiainica* computationally (user request; written 2026-11-22 12:21, before any of these tests ran)

The frozen scan result does not change. Two explanations are weighed:
- **H1 (skip):** cox1 and cob are translated, and ribosomes read past the extra bases.
- **H2 (broken):** the yeast no longer respires, cox1 and cob are dispensable, and the extra bases are disabling
  mutations, possibly scars of intron loss.

**Already seen before this entry, so not a blind prediction** (12:23 entry): the congener's codons at cox1 24711 (AGT,
Ser) and cob 26266 (CTT GGT); and the *S. cerevisiae* intron distances (3, 7, 12 nt, and > 100 nt).

**Test 1, sister codons.** At each site: the *G. yakushimaensis* codon and residue, the distant relatives' residues, and
the simplest inserted base.
- **H1 predicts** that the residue read after the skip, Gly from GGC, matches the congener at every comparable site,
  because a translated protein keeps its conserved residues.
- **H2 predicts** no such constraint. The inserted base need not be the same one each time.
- **Favours H2:** a comparable site where the congener has a residue the skip-A reading does not give, while another
  placement does (at least one such site is already known: cox1).
- **Favours H1:** Gly in the congener at all comparable sites.

**Test 2a, intron null.**
- p = the share of cox1 and cob coding positions (base boundaries) within 12 nt of an intron insertion position. It is
  computed for *S. cerevisiae* alone, then for all yeasts in the scan with annotated cox1/cob introns, mapped to *S.
  cerevisiae* protein coordinates by protein alignment.
- Reported: P(≥ 3 of 4 random sites within 12 nt) = binomial(4, p).
- **H2 (scar) predicts** P < 0.05, meaning sites cluster near intron positions.
- P ≥ 0.05 makes the intron observation uninformative. H1 makes no prediction here.

**Test 2b, per site.** The exact distance to the nearest intron junction, in the congener and in each relative with an
intron there; whether *G. jiainica* lacks an intron its relatives have at that spot; and any intron-ORF remnant within
500 nt (DIAMOND hits to intron-encoded ORFs or homing endonucleases).
- **H2 (scar) predicts** sites at or within 2 nt of the junction of an intron that *G. jiainica* lost.
- Sites 3–12 nt away, or near introns *G. jiainica* never had, do not support the scar version of H2.

**Test 3, nuclear genes.**
- Assemblies: the public nuclear assemblies of *G. jiainica* (GCA_030571715.1), *G. yakushimaensis* if one exists, and
  one relative documented to grow on glycerol or ethanol, as the control.
- **Completeness:** the share of a panel of conserved single-copy housekeeping proteins found intact.
- **Search:** *S. cerevisiae* S288C proteins with miniprot. Gene groups:
  - the cox1/cob mRNA-specific activators PET309, MSS51, CBS1, CBS2, CBP1;
  - complex III nuclear subunits COR1, QCR2, CYT1, RIP1, QCR6-10;
  - complex IV nuclear subunits COX4-9, COX12, COX13;
  - heme a synthesis, COX10 and COX15;
  - cytochrome c: CYC1, CYC7, CYC3;
  - mitochondrial translation: TUF1, MEF1 and a sample of mitochondrial ribosomal proteins.
- **Calls:**
  - *present* = the best alignment covers ≥ 70% of the *S. cerevisiae* protein with no frameshift or in-frame stop;
  - *broken* = a frameshift or in-frame stop inside the aligned region, checked against the second-best hit and,
    where needed, the reads;
  - *absent* = no alignment covering ≥ 30% of the protein;
  - otherwise *partial*.
- **H2 predicts** losses or breaks in *G. jiainica*, starting with the five activators, which exist only to make
  these two proteins. Intact genes in the control are expected.
- **H1 predicts** all intact. Intact genes in *G. jiainica* favour H1 only weakly, because a recent loss of respiration
  may not show yet.

**Test 4, congener phenotypes.** Growth on glycerol, ethanol, lactate or other non-fermentable carbon in the species
descriptions of *G. yakushimaensis*, *G. humatica* and *G. transvaalensis*, with the source and quote length.
- **Favours a *G. jiainica*-specific loss, consistent with H2:** congeners that respire while *G. jiainica* does not.
- **Weakens the link between non-growth and gene decay:** congeners that do not respire either, but have intact
  cox1/cob. Non-growth alone would then not predict decay.

**Test 5, the rest of the mitogenome** (*G. jiainica* and *G. yakushimaensis*). Every stop, indel or extra base in
atp6, atp8, atp9, cox2, cox3 and the other genes (scan events, in-frame stops, miniprot of *S. cerevisiae*
mitochondrial proteins), and the state of each tRNA.
- **H2 predicts** further disabling changes, especially in complex IV and III genes (cox2, cox3), or a decaying tRNA
  set.
- **H1 predicts** everything else intact.
- An intact remainder favours H1 weakly. ATP synthase genes may stay intact even without respiration.

**Test 6, the genus.**
- WGS run sizes for *G. transvaalensis* and *G. humatica*.
- If they total under 20 GB (user-approved), their mitogenomes are assembled de novo with SPAdes, from
  mitochondria-enriched reads, and confirmed by mapping reads back with minimap2. Then across all four *Grigorovia*
  species:
  - the cox1/cob intron map;
  - extra bases or disabling mutations in core genes (the frozen scan pipeline);
  - tRNA-Gly anticodons;
  - a branch dN/dS test (codeml branch model, per lineage), with its power stated by the number of substitutions on
    each branch.
- **H1 predicts** that AGG skips at GGC occur with a GCC tRNA-Gly wherever they occur, and purifying selection (dN/dS
  well below 1) on the *G. jiainica* branch, similar to its congeners. (This sentence completes a truncated user
  instruction.)
- **H2 predicts** extra bases only in *G. jiainica*, no tie to the tRNA, and relaxed selection (dN/dS raised toward 1)
  on the *G. jiainica* branch in cox1 and cob, and in other respiratory genes if respiration is lost.

**2026-11-22 12:55 — Exploratory test log (Grigorovia H1/H2), methods changes made while running, with reasons.**
- **Completeness (test 3).** compleasm 0.2.10 with the full saccharomycetes_odb12 references would have run about
  4 hours per genome, against the one-hour rule. It was stopped after about 5 minutes. Completeness was then measured
  with one reference protein per BUSCO gene (`data/busco_lineages_sc`). The numbers are comparative, not standard
  BUSCO scores.
- **Gene calls (test 3).** The logged miniprot rules with *S. cerevisiae* queries failed at 25–35% identity. They
  called PET309 and CBP1 "absent" in all four *Grigorovia*, though both are present.
  - Added: relatives' and congener's proteins as queries; DIAMOND blastx (ultra-sensitive, frameshift-aware) as a
    second detector; read correction of the *G. jiainica* assembly; and a local SPAdes assembly for PET122, which is
    missing from the assembly.
  - The pre-registered rule already allowed checking broken calls against the reads. Both broken calls (PET309,
    COR1) were assembly errors.
- **Introns (test 2).** The scan's insert list misplaces intron ends where an intron ORF continues the exon, and it
  missed two introns. Per-site distances were therefore taken from RefSeq annotation, a spliced nucleotide alignment
  and miniprot exon structures.
- **Tools added:** minimap2 2.31, samtools 1.25, HyPhy 2.5.104, PAML 4.10.10, sra-tools 3.4.1 (Homebrew); HMMER 3.4
  (Homebrew); SPAdes 4.3.0 (official macOS arm64 binary, sha256 2cb4a84d…; the Homebrew install failed on a leftover
  libomp directory, which was left untouched); compleasm 0.2.10 (GitHub tag, sha256 76d7bdd3…).

**2026-11-22 13:42 — Test 6 log.**
- The *G. transvaalensis* type-strain run (SRR16974369 = SRR6475934) has almost no mtDNA reads (4 of 200,000 hit
  mitochondrial proteins). It was replaced by strain CLIB 448 (ERR9706995, 0.58 GB). Total download about 1.6 GB.
- The CLIB 448 mitogenome was assembled de novo in two baiting rounds, a methods change made because the
  whole-genome assembly was fragmented. Reads were selected by similarity to the first-round contigs, then
  reassembled.
- dN/dS: codeml branch models. Added a profile-likelihood upper bound on ω for the *G. jiainica* branch, to state the
  power as the largest neutral share of the branch.
- MAFFT 7.526 installed (Homebrew).

**2026-11-22 13:48 — Local history rewritten in this session to remove a 611 MB file.**
- **What was removed.** An intermediate output of the stopped compleasm run
  (`results/explore/grigorovia/h1h2/compleasm/…/miniprot_output.gff`, 611 MB), and about 110 MB of compleasm work
  files. These were committed by mistake in 46b4e35 (was e7d72ea) and b451968 (was 05c48af).
- **How.** `git filter-branch --index-filter` over d6e06a9..HEAD, removing `results/explore/grigorovia/h1h2/compleasm`
  and `…/compleasm1`. Then the backup ref was deleted, the reflog expired and `git gc --prune=now` run; `.git`
  went from 323 MB to 17 MB.
- **What was kept.** Compact completeness tables were re-added in 68583a5, along with ignore rules. Nothing had been
  pushed. Messages, authors, and author and committer dates were preserved.
- **Cited hashes.** No hash cited in this file or in the release's `data/prereg.md` changed. Every cited commit
  predates d6e06a9. (2cb4a84d… and 76d7bdd3… in this file are file checksums, not commits.)

| old hash | new hash | date (author = committer, preserved) | commit |
|---|---|---|---|
| 841bfff | 841bfff (unchanged: no removed files) | 2026-11-22 12:24:40 | Test 1 (sister codons) |
| 4292c2f | 4292c2f (unchanged) | 12:35:57 | Test 2 (introns) |
| 14cc08a | 14cc08a (unchanged) | 12:36:45 | Test 4 (phenotypes) |
| e7d72ea | 46b4e35 | 12:40:23 | Test 2b cross-check |
| 5ab8b33 | a793b11 | 12:41:46 | Test 5 (rest of mitogenome) |
| 05c48af | b451968 | 12:55:55 | Test 3 (nuclear genes) |
| d1e4eca | 6cc342a | 13:02:13 | Test 3: spanning-read counts |
| 3f5a1fa | 0d25345 | 13:05:55 | Test 3 addendum (heterozygosity) |
| 8c64e8c | 7a7917c | 13:42:53 | Test 6 (genus) |

**2026-11-22 14:23 — Release update with the follow-up results (user request); before the second fact-check and referee
review.**
- **Numbers regenerated from final code.** Every number from the six tests was regenerated by
  `src/explore/grigorovia_h1h2.py`, `nuclear_genes.py`, `grigorovia_genus.py` and `read_haplotypes.py`, re-run end
  to end (0da95c8).
  - Values unchanged, except the read-version counts at the PET309 and COR1 assembly frameshifts. These are now
    counted in fixed windows by alignment: 21 of 23 and 23 of 25 reads in two in-frame versions, with 0 matching the
    assembly. Earlier, anchor-based counts gave 18 of 18 each.
  - Site carriage in the congeners is now called against each species' spliced CDS: 12 of 12 lack the extra base,
    including *cob* 26477, which had not been comparable across the congener's intron.
- **Corrected claim.** The tRNA-Gly(GCC) is unique among 100 other Saccharomycetaceae in the scan plus three
  congeners. The family table's 101 genomes include *G. jiainica* itself.
- **Release changes:**
  - README: the opening sentence on the candidate; the *G. jiainica* section rewritten (most likely bypassed unless
    very recently broken; caveats; tests that would settle it); the site figure kept neutral.
  - The animation (`figures/ribosome_agg.gif`, `src/figures/ribosome_agg.py`) removed from the release at the user's
    request. It stays in the working repository.
  - notes.md §8c (six-test table, false alarms, the NRRL Y-17245 caution, ploidy, downloads) and §9–11 updated.
  - `data/grigorovia/h1h2/` added; setup and run scripts updated.
  - `data/prereg.md` regenerated by `src/release/export_prereg.py`.

**2026-11-22 14:50 — Second fact-check and referee review of the release, and the changes they led to.**
- **Corrections to earlier entries.** Not edited in place:
  - **Seen-before reference.** The 12:21 entry's "(12:23 entry)" should read "(the 12:11 entries)". That is where the
    congener codons and the *S. cerevisiae* intron distances were first recorded.
  - **Identity range.** The test 3 methods entry's "failed at 25–35% identity" should read "24–32%".
  - **Hash statement.** The history-rewrite entry's "no hash cited changed" refers to hashes cited before the
    rewrite. The rewrite table itself cites new hashes.
- **Errors found in this session's follow-up work.** Logged here so that every error listed in the release notes is
  in this log:
  - an early test 1 summary put the A of AGG among the best-scoring placements at *cob* 26266 (it scores 70 against
    75);
  - pileup counts at indels suggested heterozygous frameshifts that read-level grouping showed to be alignment
    artefacts;
  - a draft called the tRNA-Gly(GCC) unique among "101 relatives", a count that included *G. jiainica* itself.
- **Checks added at the referee's request, after the tests and the first release text.** Exploratory:
  - **Codon-model and outgroup sensitivity** (`grigorovia_genus.py dnds_sensitivity`). For *cox1* + *cob*, ω = 1 is
    rejected under every setting, and the largest neutral share is 1.2–5.3%. *cox1* alone: 5–29%. The *cob*-only
    contrasts vary.
    - A first run of the no-outgroup setting used a rooted four-taxon tree, which merged the *G. jiainica* branch
      with the root branch. It was rerun unrooted.
  - **GGC codons and ancestral residues** (`ggc`). The 4 GGC codons in 1,815 are the 4 sites; Gly is ancestral at
    all four.
  - **Reads by strand and across both close *cob* sites** (`grigorovia_h1h2.py reads_strand`). Balanced; 1,516 reads
    span both sites with both bases, none without.
  - **Reciprocal best hits** (`nuclear_genes.py rbh`). All activators confirmed except PET122, which has no
    *S. cerevisiae* hit.
  - **Junctions in `compare.md`.** Refined near sites (removing a stale 10-nt value).
  - **Test 1 table.** Its "fits" column now compares residues.
  - **Test 4.** Control rows added (*M. servazzii*, *K. africana*).
  - **Strain identities** recorded in the search log (SF1S05 = CBS 10587 = ATCC MYA-4454, the type strain).
- **Release wording changed** (details in the release's `data/review_log.md`, second round):
  - **The opening sentence** now reads "survived every check that its extra bases are real". The user's wording,
    "every computational check we could run", was judged an overstatement by both reviewers; this is flagged to the
    user.
  - **The *G. jiainica* section** now says the genes were constrained until recently, and that the tests cannot tell
    a bypass from a very recent loss (including in culture). The lean toward a bypass rests on the insertion pattern,
    noticed after the fact.
  - **Caveats added:** GGA decoding, homoplasmy, orthology, and the type-strain stock.

**2026-11-22 15:34 — Final README edits (user request), after a fact-check of the edited lines.**
- **Five *G. jiainica* lines edited:** congeners at all four sites; the four not-callable genes and the controls; the
  bypass sentence; GCC pairs with GGC; the serine caveat against the logged rule.
- **Wording adjusted to the data:**
  - *cob* 26477 sits 2–3 nt from an intron, not inside one;
  - the 26266 placement qualifier;
  - "first compared", not "nearest", congener.
- **Applied after the check:** "at all four sites"; CBP2 is needed only for the bI5 intron; *M. servazzii*'s growth
  records disagree; the bypass is "proposed from DNA alone" for the CTG relatives.
- Logged in the release's `data/review_log.md`.

**2026-11-22 15:50 — Register pass on the release's README.md and notes.md with the `declauding` skill** (user request; skill
installed to the user's personal skills folder, version 0.17.0).
- **Steps:** lint; sentence and structure pass; a separate-agent structural review; two separate-agent newcomer
  reads; a claim diff against the pre-edit files (no numbers or claims lost).
- **Factual fixes found during the pass** (details in the release's `data/review_log.md`):
  - the *Verticillium* second-assembly wording (38, and 37 of 37 comparable);
  - the CTG count against P3;
  - the README's "no trace" against the notes' 1–5% bound;
  - the headline's implied weak flags;
  - the congener intron "at" rather than 2–3 nt from the TGG site;
  - the losses table's unshown categories.
- **Terms renamed and defined:** "CTG yeasts" renamed "CTG-shifting yeasts" (the "CTG clade" is a different group);
  terms defined at first use.
- **Not edited:** sentences the user dictated; the process record in notes.md.

**2026-11-22 19:41 — Release reorganised (user request): file renames only, no content changes.**
- `notes.md` → `notebook.md`; `data/` split into `protocol/` (frozen rules, corrections, design, prior art,
  accessions, scan, controls, audit) and `candidate/` (*G. jiainica* summary, funnel, follow-up tests, held-out
  congener, literature search); `src/` → `code/` (`code/plots/` for the figure script); `figures/` → `plots/`.
- **Paths updated.** Internal paths and links in the release files were updated. Exceptions, kept verbatim as
  records: `protocol/frozen_rules.md` (this log) and `protocol/design.md`. The two `notes.md` mentions in
  `protocol/prior_art.md` name a file in ry-xu/genome-research.
- **Download location.** The scripts now download the congener record to `data/raw/OR703946/`, which avoids
  clashing with the moved `candidate/heldout/`.
- **README text.** The three disclaimer bullets were replaced by a "How this was made" paragraph with the same
  content, and a credit line for the write-up format (ry-xu/genome-research) was added.
- **Rename table.** The old-to-new table is in the release's `protocol/corrections.md`.

---

## Commits touching PREREG.md or DESIGN.md (local git; hash, commit time EST, message)

| commit | time | message |
|---|---|---|
| 6ccafb5 | 2026-11-21 19:11 | Record design decisions agreed before any data or code |
| a8975e8 | 2026-11-21 20:21 | Prior-art audit answer 7 (red team), verified additions; design amendment 1; author email draft |
| bacd304 | 2026-11-21 20:28 | Remove author email draft; add standing no-outside-contact rule (CLAUDE.md); run MFannot locally |
| 6c79d41 | 2026-11-21 20:35 | DESIGN amendment 3 and PREREG section 0: user's conditions, single confirmatory test, stop rule |
| 57c8b8a | 2026-11-22 06:48 | DESIGN amendment 4 and PREREG log: user decisions after first control run |
| 942d74d | 2026-11-22 06:48 | PREREG log: correct NC_079576 identity figures |
| e878133 | 2026-11-22 06:56 | Re-run code per amendment 4: window rule with cap, split-hit shifts, 2/3/5-site trials, gates/limits, strong flags; B2 diagnosis corrected |
| 43f9003 | 2026-11-22 07:10 | PREREG log: split-hit probation rule, misplacement reporting, pool site tolerance 50 residues, >=3-gene tally (logged before re-run numbers) |
| 4e73884 | 2026-11-22 07:10 | PREREG log: correct typed timestamps to commit times |
| 9408596 | 2026-11-22 07:13 | Length filter not kept (pre-check: coverage +0.001, B9a +0.018); bars/flags for split probation and misplacement reporting |
| 8140435 | 2026-11-22 07:50 | Freeze PREREG section 1: control re-run (gates pass; B2 limit), rules, sensitivity, predictions |
| e085d3a | 2026-11-22 07:56 | PREREG log: post-freeze descriptive flag context and prediction scoring (before scan output) |
| f5fd793 | 2026-11-22 08:12 | PREREG log: post-freeze many-shift tally reporting and known-editor list (before full-scan output) |
| c404029 | 2026-11-22 08:54 | PREREG log: per-site source column for known-lineage flags (post-freeze); grebe look not started |
| dbccf9c | 2026-11-22 09:14 | PREREG log: completeness, assembly unit for strong rule, real-data sensitivity, per-flag own chance (before analysis) |
| 07327e5 | 2026-11-22 09:21 | Report additions: completeness, assembly units (strong rule), real-data lineage sensitivity, per-flag own chance, independent assemblies |
| c991c0c | 2026-11-22 09:24 | PREREG log: conditional VdLs.17 second-assembly step, write-up limits, paired-break tally, strain-only name merge (before analysis) |
| 62c56a5 | 2026-11-22 09:34 | PREREG log: per-lineage real rate with nad3-174 separate, headline excludes it; paired tally context (before analysis) |
| 076ae07 | 2026-11-22 09:41 | PREREG log: removed-break listing for sponges/Polyrhachis, out-of-sample sensitivity (before analysis); held-out run not started |
| 7f54442 | 2026-11-22 10:03 | PREREG log: Grigorovia sent to the funnel after seeing the result (exploratory); order of follow-ups |
| bd7f02b | 2026-11-22 10:31 | PREREG log: Grigorovia read pass rule (before reads counted); owed checks |
| 7164bea | 2026-11-22 10:43 | PREREG log: Grigorovia reads pass 4/4; blinded check design (decoys, refs, prompt) before running |
| 74f5edd | 2026-11-22 11:04 | Grigorovia blinded check: candidate's 4 sites and Gly residue recovered blind; decoys as expected |
| ea426e3 | 2026-11-22 11:33 | PREREG log: OR703946 (G. yakushimaensis) predictions before opening; respiration caveat |
| 313292f | 2026-11-22 11:36 | OR703946 (G. yakushimaensis) result: no extra bases at homologous sites, no AGG shifts, tRNA-Gly UCC; cob intron at the TGG site |
| 7297316 | 2026-11-22 11:38 | Post hoc tRNA-Gly checks: GCC anticodon (1,329 reads), A37, 7-nt loop; only GCC tRNA-Gly in 101 Saccharomycetaceae |
| 7e36bb9 | 2026-11-22 11:50 | PREREG log: release draft and reviews launched |
| fcc3f0c | 2026-11-22 11:54 | Correct Grigorovia SRA search (species split across two genus taxids): WGS reads exist for all three sisters, no RNA |
| 271768c | 2026-11-22 11:57 | PREREG log: ribosome_agg.gif |
| e9d9387 | 2026-11-22 12:11 | PREREG log: fact-check corrections (blind references not independent; time reference; process slips); referee-prompted checks and their consequences |
| 7be11ce | 2026-11-22 12:11 | PREREG log: mark reconstructed times as approximate |
| d6e06a9 | 2026-11-22 12:21 | PREREG exploratory: H1/H2 predictions for six Grigorovia tests (before running) |
| b451968 | 2026-11-22 12:55 | Test 3 (nuclear genes): all callable respiration and cox1/cob activator genes present and intact in G. jiainica after read correction; methods changes logged |
| 7a7917c | 2026-11-22 13:42 | Test 6 (genus): congener mitogenomes assembled (humatica; transvaalensis CLIB 448; type strain rho0-like); no extra bases in congeners; tRNA-Gly TCC in all three; strong purifying selection on G. jiainica cox1/cob (omega 0.031, neutral share <= 2.5%) |
| edadaa2 | 2026-11-22 13:48 | PREREG log: local history rewrite (611 MB file removed), old-to-new hash table, no cited hash changed |
| 1ad1f26 | 2026-11-22 14:23 | PREREG log: release update with follow-up results; prereg export script; downloads table |
| 1d516bf | 2026-11-22 14:50 | Second review round: sensitivity of the omega test, GGC pattern, reads by strand, reciprocal best hits, refined junctions; PREREG corrections and log |
| 9a7e3df | 2026-11-22 15:34 | PREREG log: final README edits after the line-level fact-check |
| 0222762 | 2026-11-22 15:50 | PREREG log: declauding register pass on the release README and notes |
| 67690f7 | 2026-11-22 19:41 | PREREG log: release reorganised into protocol/ and candidate/ (renames only); export header notes the path mapping |

## All project commits (local git)

| commit | time | message |
|---|---|---|
| 6ccafb5 | 2026-11-21 19:11 | Record design decisions agreed before any data or code |
| da34498 | 2026-11-21 19:29 | Step 1: prior art on recurring-codon mitochondrial frameshifts |
| 6f8390e | 2026-11-21 20:02 | Prior-art audit: answers 1-6, corrected PRIOR_ART.md, search log, annotation census |
| a8975e8 | 2026-11-21 20:21 | Prior-art audit answer 7 (red team), verified additions; design amendment 1; author email draft |
| bacd304 | 2026-11-21 20:28 | Remove author email draft; add standing no-outside-contact rule (CLAUDE.md); run MFannot locally |
| 6c79d41 | 2026-11-21 20:35 | DESIGN amendment 3 and PREREG section 0: user's conditions, single confirmatory test, stop rule |
| fecb8eb | 2026-11-21 20:50 | Scan core, prep, and control genome sets fixed before any control run |
| 8025cfc | 2026-11-21 20:58 | Controls script, flag rule, intron overlap rule; before the full control run |
| 972ae02 | 2026-11-21 21:32 | Full control run: bars missed (B5, B2, B9a, B9c); diagnostics and offline window-convention evaluation |
| 102b464 | 2026-11-22 06:38 | PRIOR_ART: log Xu 2026-11-21 mitochondrial codon census (CUU/CUA as Pro); adjacent, not a frameshift result |
| 57c8b8a | 2026-11-22 06:48 | DESIGN amendment 4 and PREREG log: user decisions after first control run |
| 942d74d | 2026-11-22 06:48 | PREREG log: correct NC_079576 identity figures |
| e878133 | 2026-11-22 06:56 | Re-run code per amendment 4: window rule with cap, split-hit shifts, 2/3/5-site trials, gates/limits, strong flags; B2 diagnosis corrected |
| 43f9003 | 2026-11-22 07:10 | PREREG log: split-hit probation rule, misplacement reporting, pool site tolerance 50 residues, >=3-gene tally (logged before re-run numbers) |
| 4e73884 | 2026-11-22 07:10 | PREREG log: correct typed timestamps to commit times |
| 9408596 | 2026-11-22 07:13 | Length filter not kept (pre-check: coverage +0.001, B9a +0.018); bars/flags for split probation and misplacement reporting |
| cfe26e7 | 2026-11-22 07:17 | Fix site-codon frame (from CDS) and B2 in-gene lookup; analyze: 50-residue pool sites, split regime, 3-gene tally |
| 6397e92 | 2026-11-22 07:17 | Clear first-run stage outputs before the re-run (kept in history at 972ae02) |
| 8140435 | 2026-11-22 07:50 | Freeze PREREG section 1: control re-run (gates pass; B2 limit), rules, sensitivity, predictions |
| 44a79fe | 2026-11-22 07:51 | gitignore scan hits |
| e085d3a | 2026-11-22 07:56 | PREREG log: post-freeze descriptive flag context and prediction scoring (before scan output) |
| 69cccb5 | 2026-11-22 07:58 | flag_context.py: post-freeze descriptive context per flag and prediction scoring |
| a7044a3 | 2026-11-22 08:09 | PRIOR_ART: log missed heteroloboseid RNA editing (Yang 2017), seen in the scan |
| a2e2607 | 2026-11-22 08:09 | PRIOR_ART: fix gene count in Yang 2017 entry (13, not 12) |
| f5fd793 | 2026-11-22 08:12 | PREREG log: post-freeze many-shift tally reporting and known-editor list (before full-scan output) |
| cea005b | 2026-11-22 08:14 | flag_context: many-shift tally listing with direction share, platform, related pattern, known editors (post-freeze, descriptive) |
| c404029 | 2026-11-22 08:54 | PREREG log: per-site source column for known-lineage flags (post-freeze); grebe look not started |
| d96fdd4 | 2026-11-22 08:56 | flag_context: per-site paper (PRIOR_ART only) vs record-annotation column for known-lineage flags |
| dbccf9c | 2026-11-22 09:14 | PREREG log: completeness, assembly unit for strong rule, real-data sensitivity, per-flag own chance (before analysis) |
| 07327e5 | 2026-11-22 09:21 | Report additions: completeness, assembly units (strong rule), real-data lineage sensitivity, per-flag own chance, independent assemblies |
| c991c0c | 2026-11-22 09:24 | PREREG log: conditional VdLs.17 second-assembly step, write-up limits, paired-break tally, strain-only name merge (before analysis) |
| e9a0791 | 2026-11-22 09:28 | Conditional VdLs.17 second-assembly check (tested on planted indels), paired-break tally, real-rate limits, strain-only units |
| 62c56a5 | 2026-11-22 09:34 | PREREG log: per-lineage real rate with nad3-174 separate, headline excludes it; paired tally context (before analysis) |
| 614d6b7 | 2026-11-22 09:37 | flag_context: per-lineage real rates (nad3-174 separate, headline excludes it), paired-tally totals and alignment identity |
| 076ae07 | 2026-11-22 09:41 | PREREG log: removed-break listing for sponges/Polyrhachis, out-of-sample sensitivity (before analysis); held-out run not started |
| 5f1238f | 2026-11-22 09:49 | flag_context: removed-break listing for glass sponges/Polyrhachis (TGGA signature, documented status), extra known groups, out-of-sample sites |
| 43fae2b | 2026-11-22 09:58 | Full scan (20,343 genomes) and report: 15 known-lineage flags, 1 non-known strong flag (VdLs.17) killed at second-assembly step; predictions P2, P4 fail |
| 21ba959 | 2026-11-22 09:59 | REPORT: exact identity range; add fetched Haen 2014 abstract and Oopsacas paper records |
| 7f54442 | 2026-11-22 10:03 | PREREG log: Grigorovia sent to the funnel after seeing the result (exploratory); order of follow-ups |
| bd7f02b | 2026-11-22 10:31 | PREREG log: Grigorovia read pass rule (before reads counted); owed checks |
| b775497 | 2026-11-22 10:36 | Grigorovia exploratory steps a-e, AGG usage in relatives, completeness, second tRNA finder, data search |
| 7164bea | 2026-11-22 10:43 | PREREG log: Grigorovia reads pass 4/4; blinded check design (decoys, refs, prompt) before running |
| 452342e | 2026-11-22 10:43 | Blind check key (written before the blinded agent runs) |
| 84876da | 2026-11-22 10:44 | Grigorovia: lab-ready site table (strain, coordinates, context, read support) |
| 74f5edd | 2026-11-22 11:04 | Grigorovia blinded check: candidate's 4 sites and Gly residue recovered blind; decoys as expected |
| ea426e3 | 2026-11-22 11:33 | PREREG log: OR703946 (G. yakushimaensis) predictions before opening; respiration caveat |
| 313292f | 2026-11-22 11:36 | OR703946 (G. yakushimaensis) result: no extra bases at homologous sites, no AGG shifts, tRNA-Gly UCC; cob intron at the TGG site |
| 7297316 | 2026-11-22 11:38 | Post hoc tRNA-Gly checks: GCC anticodon (1,329 reads), A37, 7-nt loop; only GCC tRNA-Gly in 101 Saccharomycetaceae |
| e375b00 | 2026-11-22 11:43 | Consolidate Grigorovia one-off follow-up commands into src/explore/grigorovia_followup.py (reproduces data files) |
| 00ccb8c | 2026-11-22 11:46 | Figures: Grigorovia cob double-site figure and ribosome animation; figure data subcommand |
| 2662381 | 2026-11-22 11:46 | Save AGG-placement residue check and relatives ranking to data files |
| 6711345 | 2026-11-22 11:50 | Release draft (before fact-check and referee review) |
| 7e36bb9 | 2026-11-22 11:50 | PREREG log: release draft and reviews launched |
| e3cbf85 | 2026-11-22 11:51 | Ignore downloaded reads; add held-out summary |
| fcc3f0c | 2026-11-22 11:54 | Correct Grigorovia SRA search (species split across two genus taxids): WGS reads exist for all three sisters, no RNA |
| 9839c5e | 2026-11-22 11:57 | Animated model of the AGG skip (src/figures/ribosome_agg.py, figures/ribosome_agg.gif) with pairing checks; static figure script moved to src/figures/ |
| 271768c | 2026-11-22 11:57 | PREREG log: ribosome_agg.gif |
| 8b8ca9f | 2026-11-22 12:07 | Post hoc: intron positions near 3 of 4 Grigorovia sites; congener has Ser at the cox1 site (referee points verified) |
| e9d9387 | 2026-11-22 12:11 | PREREG log: fact-check corrections (blind references not independent; time reference; process slips); referee-prompted checks and their consequences |
| 7be11ce | 2026-11-22 12:11 | PREREG log: mark reconstructed times as approximate |
| 40535cc | 2026-11-22 12:12 | Release draft after fact-check and referee review (see release_draft/data/review_log.md) |
| d6e06a9 | 2026-11-22 12:21 | PREREG exploratory: H1/H2 predictions for six Grigorovia tests (before running) |
| 841bfff | 2026-11-22 12:24 | Test 1 (sister codons): skip-A/GGC=Gly fits the congener at cob 26496 and 26266, not at cox1 (congener Ser); 26477 at congener intron |
| 4292c2f | 2026-11-22 12:35 | Test 2 (introns): null 0.004 (S. cerevisiae) / 0.23 (all yeasts); no site within 2 nt; congener bI2-type intron 3 nt from cob 26477; no ORF remnants |
| 14cc08a | 2026-11-22 12:36 | Test 4 (phenotypes): G. yakushimaensis and G. humatica also negative on non-fermentable carbon; transvaalensis uncertain |
| 46b4e35 | 2026-11-22 12:40 | Test 2b cross-check: miniprot exon structures; M. solicola also has a cox1 intron 9 nt from the site; congener cox1 intron far; summary corrected |
| a793b11 | 2026-11-22 12:41 | Test 5 (rest of mitogenome): all other genes intact in G. jiainica; tRNA set matches congener except Gly GCC/TCC |
| b451968 | 2026-11-22 12:55 | Test 3 (nuclear genes): all callable respiration and cox1/cob activator genes present and intact in G. jiainica after read correction; methods changes logged |
| 6cc342a | 2026-11-22 13:02 | Test 3: spanning-read counts at the PET309 and COR1 assembly frameshifts (18/18 intact each; PET309 has two in-frame repeat versions) |
| 0d25345 | 2026-11-22 13:05 | Test 3 addendum: G. jiainica nuclear genome heterozygous (~1.2%); both copies of PET309, COR1, PET111 intact by read haplotypes |
| 7a7917c | 2026-11-22 13:42 | Test 6 (genus): congener mitogenomes assembled (humatica; transvaalensis CLIB 448; type strain rho0-like); no extra bases in congeners; tRNA-Gly TCC in all three; strong purifying selection on G. jiainica cox1/cob (omega 0.031, neutral share <= 2.5%) |
| 68583a5 | 2026-11-22 13:44 | Keep compact completeness tables; ignore large intermediates (assemblies, BAMs, lineage data, compleasm work files) |
| 2ac5d8f | 2026-11-22 13:44 | Grigorovia H1/H2: summary of the six exploratory tests |
| edadaa2 | 2026-11-22 13:48 | PREREG log: local history rewrite (611 MB file removed), old-to-new hash table, no cited hash changed |
| 0da95c8 | 2026-11-22 14:19 | Consolidate every Grigorovia H1/H2 number into final scripts (grigorovia_h1h2.py, nuclear_genes.py, grigorovia_genus.py, read_haplotypes.py) and re-run end to end; intron-aware site calls (12/12 lack); write-ups cite the regenerated tables |
| 1ad1f26 | 2026-11-22 14:23 | PREREG log: release update with follow-up results; prereg export script; downloads table |
| 1d516bf | 2026-11-22 14:50 | Second review round: sensitivity of the omega test, GGC pattern, reads by strand, reciprocal best hits, refined junctions; PREREG corrections and log |
| d809257 | 2026-11-22 14:51 | Release draft after the second fact-check and referee review (snapshot of release/mito-frameshift-scan) |
| 9a7e3df | 2026-11-22 15:34 | PREREG log: final README edits after the line-level fact-check |
| a6aae2e | 2026-11-22 15:34 | Release draft: final README edits to the Grigorovia section (snapshot) |
| 0222762 | 2026-11-22 15:50 | PREREG log: declauding register pass on the release README and notes |
| 04cd84c | 2026-11-22 15:50 | Release draft after the declauding register pass (snapshot) |
| 67690f7 | 2026-11-22 19:41 | PREREG log: release reorganised into protocol/ and candidate/ (renames only); export header notes the path mapping |
