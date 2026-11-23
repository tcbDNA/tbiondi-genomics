# Prior art: frameshifts at recurring codons in mitochondrial genomes

Step 1 of the project (see `DESIGN.md`). First written 2026-11-21 19:29 EST, revised 19:57–20:10 EST
after an audit. The audit's questions and answers are in `PRIOR_ART_AUDIT.md`, the search log in
`lit/SEARCH_LOG.md`, and the corrections in the log at the end of this file. No genome sequence had been
scanned when this was written. RefSeq feature tables (annotation only) were read for the audit.

**How each source was read.**
- **F-agent:** full text read end to end by a reading agent working from the saved file; the key
  sentences quoted here were re-checked by me against the file.
- **K:** keyword passages of the full text read by me (or by the sweep agent, marked "K-sweep").
- **A:** abstract only.
- **C:** known only from a citation in another paper.

Saved sources are in `lit/`: Europe PMC metadata `*.json`, open full texts `*.xml`, and `al_arab_2017/`
for that paper's author manuscript and supplement.

## Summary

- **Four lineages are known to recode at recurring codons in mitochondria.** Each was found by hand in a
  few genomes, and in each the trigger codon is rare or unused in frame elsewhere in the genome:
  - glass sponges: UGG, +1;
  - *Perkinsus*: AGGY +1, CCCCU and TMGGY +2;
  - a clade of budding yeasts (*Kazachstania*, *Nakaseomyces*): CTGY, +1C;
  - ***Magnusiomyces***: 81 translational bypassing elements, ~50 nt each, whose "takeoff" is at unused
    codons.
- **A weaker recurring motif is described in turtles.** Shifts in *nad3*, *nad4l* and *nad4* sit at a
  CTN codon followed by a rare AGN codon (Russell & Beckenbach 2008).
- **Single or few sites are known** in birds (*nad3*-174), *Polyrhachis* ants, an oyster, venerid bivalves,
  a diatom, a dinoflagellate parasite, a strepsipteran and two chytrid fungi (the chytrids from database
  annotation only).
- **The human "−1 frameshift at AGA/AGG" is superseded:** release factors read those codons directly.
- **RefSeq already writes 685 frameshift junctions into 632 mitochondrial records** (release 238), 616 of
  them at bird and turtle *nad3*-174. The 2,695 fungal mitogenomes of Sarhan et al. 2025 carry no gene
  annotation at all.
- **Scans so far:**
  - Al Arab et al. 2017: animals only, 3,842 genomes, one shift per gene.
  - MitoZoa 2.0: animals, frameshifts flagged during curation.
  - Wolters et al. 2023: 353 budding-yeast mitogenomes annotated with MFannot and curated by hand.
    +1C frameshifts are mentioned in one sentence; bypass elements were searched for with Lang 2014's
    method. No frameshift scan is described.
  - MFannot flags "potential frameshifts" in every genome it annotates.

  None tested for shifts piling up at one codon across lineages, and none checked candidates against
  reads, a second assembly and RNA.
- **Editing-aware gene finders** (Gott 2005, Beargie 2008, Chen 2012) are computational frameshift-aware
  scans of myxomycete mitogenomes, for insertional RNA editing.
- **A fifth recurring-codon lineage, found by the red team** (20:20): *Heteranomia* (bivalves, Anomiidae)
  shifts +1 at every in-frame TAG in *nad4* and *nad6* (Lin et al. 2026, bioRxiv, August 2026), DNA only.
  Round 1 said nothing had been found from 2025–2026; that was wrong.

## 1. Lineages with recoding at recurring codons (the "prize" type, already known)

| Lineage | Trigger → effect | Sites | Evidence beyond DNA | Sources |
|---|---|---|---|---|
| **Glass sponges** (Hexactinellida) | UGG (Trp) → +1; the GGA after it is read as Gly | *Aphrocallistes vastus*: 2 (*cox3* aa 140, *nad6* aa 58), the **only two UGG codons in the genome**; all 72 other Trp are UGA. Others: 9 new insertions in 8 genomes (Haen 2014), *Lophophysema* *cob*, *Docosaccus* *nad2* aa 63 (after its only UGG). | cDNA across both *A. vastus* sites keeps the base, "unpublished data", no traces | Rosengarten 2008 (F-agent); Haen 2014 (A); Zhang 2016 (A); Kahn & Geller 2018 (K) |
| ***Perkinsus*** (Myzozoa) | AGGY → +1 (Gly); TMGGY → +2 (Gly); CCCCU → +2 (Pro); one predicted shift at TAA | 75 in 4 species; per species *cox1* 10, *cob* 8, *cox3* 1; conserved in position. "At the point of each frameshift a codon that is otherwise unused occurs in the current translation frame." | RNA-seq, transcripts and Sanger cDNA show no editing. No tRNAs are encoded. A cryo-EM structure of the frameshifted proteins exists (preprint). | Masuda 2010 (A); Zhang 2011 (K); Gornik 2022 (F-agent); Wú 2024 (A) |
| **Budding yeasts** (*Kazachstania* incl. *Huiozyma*, *Monosporozyma*; *Nakaseomyces delphensis*, *N. glabratus*) | CTGT/CTGC → +1C, read as Cys | *K. sinensis*: *cox3*, *atp6*, *cob*. *K. servazzii*: 5 sites in *atp6*, *cox2*, *cox3*, *cob*. "Several *Kazachstania* species" among 353 yeast mitogenomes (Wolters 2023). | DNA only. Proposed quadruplet reading by tRNA-Cys with A37→G. | Szabóová 2018 (A); Langkjaer 2003 (A); Wolters 2023 (F-agent) |
| ***Magnusiomyces*** (Saccharomycetes) | Translational **bypassing**: blockage at a "takeoff" codon (unused codons, not stops), skip ~50 nt, resume at a matching "landing" triplet | 81 byps in *M. capitatus* (NC_024095), each written in RefSeq as a CDS join with a 27–~60 nt gap. "Byp-like sequences" in several Saccharomycetales lack bypassing activity. | "by transcript analysis and proteomics, that byps are retained in mitochondrial mRNAs but not translated" (3 cases) | Lang 2014 (A); Nosek 2015 (A); Nosek & Tomáška 2025 (A) |
| **Anomiid bivalves** (*Heteranomia*; APPD lineage) | In-frame TAG (a stop in code 5) → +1, skipping the T. "The invariant AGN codon immediately following each frameshift suggests potential functional significance". Relatives (*Anomia*, *Enigmonia*, *Placuna*, Plicatulidae) instead reassign TAA to Tyr. | *H. aculeata* *nad4* and *nad6*, "each containing four in-frame TAG codons with four corresponding frameshifts" | DNA only for the frameshift: "only dried specimens were available, precluding transcriptomic or functional assays". The TAA→Tyr reassignment in relatives is backed by proteomics. | Lin et al. 2026 bioRxiv (K, me; found by the red team) |
| *(weak)* **Turtles** (Testudines) | +1 at CTN followed by a rare AGN (AGT/AGA) | 1–3 per genome: *nad3*-174, *nad4l*, *nad4*, *nad3*-135 (*Pelomedusa*) | *Pelusios castaneus* *nad3*: DNA and RNA reads both carry the insertion (one run each) | Russell & Beckenbach 2008 (F-agent); Andreu-Sánchez 2021 (F-agent) |

**Status in our data:**
- All four prize-type lineages are **in the scan set** (RefSeq 238): *A. vastus* NC_010769,
  *P. marinus* NC_068742, *P. olseni* NC_068743, *Huiozyma sinensis* NC_036748, *Monosporozyma servazzii*
  NC_004918, *M. unispora* NC_036378, *N. glabratus* NC_004691, *M. capitatus* NC_024095.
- At least 11 more *Kazachstania*/*Nakaseomyces*/*Monosporozyma*/*Huiozyma* genomes sit in the fungal set,
  unannotated.
- Eight more glass-sponge mitogenomes are in GenBank only (held-out set).
- Every case, with accession, position, motif and evidence level, is in
  `results/audit/controls_known_cases.tsv`.

## 2. Single or few programmed sites (controls, not prizes)

| Site | Lineage | Evidence | Sources |
|---|---|---|---|
| *nad3*-174 (+1; extra C/T/G after CTN, before AGT AGC) | Birds and turtles: 811 taxa, ~33–38 gains and 47–52 losses; absent in Passeriformes, crocodiles, squamates. 619 RefSeq records write it. | Chicken RT-PCR keeps the base (stated, no figure); others DNA only | Mindell 1998 (C); Russell & Beckenbach 2008 (F-agent); Al Arab 2017 (F-agent); Andreu-Sánchez 2021 (F-agent); Yuan 2025 (K) |
| *cob*, 1–2 sites per species (e.g. TGG AGT A) | *Polyrhachis* ants, 12 of 30 species | DNA only | Beckenbach 2005 (A) |
| *cob* | *Crassostrea virginica* (AY905542; not in RefSeq 238) | DNA only | Milbury & Gaffney 2005 (A) |
| *nad4* (3 species), *cob* (1) | Venerid bivalves (*Leukoma*, *Anomalodiscus*, *Antigona*) | DNA only; "first time ... in Heterodonta" | Wang 2020 (A) |
| *cytb* +1 inside a T run | *Arctica islandica* (NC_022709) | annotation only | RefSeq |
| *cox1*, conserved frameshift | Foraminifera, order Miliolida: all 4 species analyzed, no rotaliid | DNA only (single-cell genomics). "suggests that a mechanism for stop codon read-through or posttranscriptional mRNA modification of this codon exists in this lineage" | Macher et al. 2023 (K, me; found by the red team) |
| *cox1*, +1 at CCC/CCU | Diatom *Phaeodactylum tricornutum* (NC_016739; not annotated) | DNA only | Oudot-Le Secq & Green 2011 (A) |
| *cox1* fragment, +1 at GGC | *Amoebophrya* sp. A120 (not A25) | DNA only | Kayal & Smith 2021 (K) |
| *nad5*, −1 | Strepsipteran *Dipterophagus daci* (NC_058598; not annotated) | 6 assemblies + 5 unamplified Sanger samples | Towett-Kirui 2022 (F-agent) |
| *cob*, +1 | Chytrid fungi *Rhizophydium* sp. 136 (/ribosomal_slippage) and *Spizellomyces punctatus* (note) | annotation only; no paper found. Forget 2002, cited by the records, does not mention it in its abstract. | RefSeq; Forget 2002 (A) |
| Other single annotations | *Hackeriella* (insect) *nad5*, *Daphnia magna* *cytb* −1 (/ribosomal_slippage), *Tachybaptus* (grebe) *nad1*, *Ilyanassa* (snail) and *Argulus* (crustacean) notes | annotation only | RefSeq |
| *ATP8/ATP6* −1 at a slippery sequence | Human | ribosome profiling (second-hand) | Moran 2024 via Wakigawa 2026 (C) |

## 3. Superseded: human −1 at AGA/AGG

Temperley et al. 2010 (A) proposed a −1 shift before the AGA/AGG ends of human *MT-CO1* and *MT-ND6*.
Later work shows that the release factors mtRF1 and mtRF1L recognize AGA/AGG directly:
- ribosome profiling, *in vitro* assays and cryo-EM, reviewed by Krüger et al. 2024 (F-agent);
- knockouts in Wakigawa et al. 2026 (F-agent).

Krüger et al. note that 3,497 of 7,188 vertebrate records lack the U that a −1 shift would need. **AGA/AGG
gene ends are not a control, and shifts at gene ends are excluded from the codon test.**

## 4. Earlier scans, tools, and how databases write frameshifts

**Scans:**
- **Al Arab et al. 2017** (F-agent; author manuscript and supplement in `lit/al_arab_2017/`). This is the
  method behind MITOS2's training alignments.
  - **Data:** "The data of the 3842 complete metazoan mtDNA sequences that are contained in RefSeq release
    63 serve as training data set."
  - **Detection:** each gene is translated with every single and double nucleotide deleted. Variants are
    scored against a Metazoa HMM, and outliers kept by a Grubbs test with Δs ≥ 25.
  - **Results:** all 200 annotated RefSeq frameshifts were found ("98.5%" in the summary). 85 literature
    cases were recovered. 36 new shifts turned up in 12 genes, "not phylogenetically conserved ...
    Therefore we cannot determine whether real frameshifts have occurred in different sites in the tree,
    or the reading frame change is caused by sequencing or annotation errors."
  - **Limits:** "our method only searches for cases with a single frameshift"; 3′-end shifts are hard to
    detect.
  - **Not done:** no codon test beyond *nad3*-174, no pooling, no reads or RNA, no non-animals.
  - **Published list:** accession, gene and family only (`supplements/12/fsnewtax.dat`).
- **MitoZoa 2.0** (D'Onorio de Meo 2012, K): "198 CDS with frameshift sites". Most were judged sequencing
  errors near homopolymers longer than 8 nt; others "cannot be easily explained".
- **Wolters et al. 2023** (F-agent; methods re-read in full by me, 20:08). What it did, exactly:
  - **Genomes:** it built and annotated 353 Saccharomycotina mitogenomes: 110 curated GenBank references
    plus mitochondrial contigs from 332 nuclear assemblies and reassemblies of 288 species.
  - **Annotation:** MFannot was used for new assemblies and for references without annotations. All
    annotations were then curated by hand.
  - **Frameshifts:** the paper **does not describe a frameshift scan**. The +1C frameshifts are reported in
    one sentence of the annotation methods, with no method, count or species list.
  - **Bypass elements:** these were searched for in the other species' coding sequences "using methods
    previously described (Lang et al., 2014)", with no further detail.

  The passage (Materials and methods, "Mitochondrial genome mining, assembly, and annotation", paragraphs
  1 and 3; `lit/wolters_2023.xml`):
  > "When annotations were not available, new annotations were generated using MFANNOT (Lang et al.,
  > 2007)." ... "All new assemblies, as well as existing mtDNAs that were not annotated, were annotated from
  > scratch; all genome annotations, including published ones, were curated for consistency and to improve
  > accuracy as described below." ... "Finally, all annotations (for new and existing mtDNAs) were further
  > manually curated to eliminate truncated genes, annotations split across contigs, and annotations
  > containing large extensions due to misannotated introns or readthroughs." ... "We identified several
  > Kazachstania species with frameshifts consistent with the +1C frameshift mechanism previously described
  > (Szabóová et al., 2018). To match the formatting in GenBank for those references, we encoded these as
  > single-bp introns, but these were excluded from all intron analyses. We did not observe any byp
  > elements, as described in Magnusiomyces capitatus, in the coding sequences of other species using
  > methods previously described (Lang et al., 2014)."
- **Lang et al. 2014** (A): byp-like sequences detected in "mtDNAs of several Saccharomycetales".
- **Andreu-Sánchez et al. 2021** (F-agent): *nad3* only, 87,707 sequences (10,309 taxa). It compared *nad3*
  annotations of 101 bird species between GenBank and independent B10K assemblies, and checked one turtle
  site against DNA and RNA reads.
- **Adjacent, not mitochondrial:** Poehls et al. 2025 (bioRxiv, A), frameshift sites across bacteria
  (SLIPPERRS).
- **Contested:** Seligmann (2011–2019; A). Computational predictions of "tetracodons" (four-base codons)
  in mitochondrial genes, in more than one taxon:
  - Seligmann 2012: "Occasional single tetracodons within regular mitochondrial genes coevolve
    positively/negatively with antisense tRNAs with predicted reduced/expanded anticodons (depending on
    taxon)".
  - Seligmann & Labra 2013: across Lepidosauria.
  - Seligmann & Warthi 2019: peptides in human mass-spectrometry data.

  No independent replication was found. The first version of this file called it "single species", which
  understated it.
- **Editing-aware scans outside animals** (red team):
  - Gott, Parimi & Bundschuh 2005 (*NAR* 33:5063; K, me): "our new algorithm Predictor of Insertional
    Editing (PIE) can be used to locate genes whose mRNAs are subjected to multiple frameshifting events"
    in *Physarum* mitochondria.
  - Beargie et al. 2008 (A): genome annotation "in the presence of insertional RNA editing".
  - Chen et al. 2012 (K, me): editing-site prediction in *Didymium*.

  These are frameshift-aware computational scans of non-animal mitogenomes (one or two genomes, for
  editing rather than translational shifts).

**Tools** (code and papers checked by me, 2026-11-21):
- **MFannot** (code `BFL-lab/Mfannot`, commit 0990c85, 2023-10-02) flags frameshifts. When two BLAST
  matches on one strand lie ≤ 142 nt apart, it prints "Potential frameshift in <gene>" and writes
  "potential frameshift" into its annotation. It also refuses to translate gene models whose length is not
  a multiple of 3. Lang et al. 2023 (K): MFannot will "identify and flag potential frameshifts indicative of
  either a pseudo-gene or sequencing errors".
- **MITOS2** (code `gitlab.com/Bernt/MITOS`, commit 0f7f916, 2025-12-20) has no frameshift flagging in the
  annotation path. A frame difference is computed only in `scripts/analyse.py`, which compares two
  annotations. Donath et al. 2019 (K): "the computational prediction of frameshifts close to the 3′ end of
  PCGs is not a trivial task".
- **MitoHiFi** (K) prefers final contigs whose "annotation includes no genes that contain a frameshift".
- **MitoFinder** (K): in its ant case study the authors excluded "Sequences with internal stop codons ...
  potentially corresponding to nuclear mtDNA segments (NUMTs)". That was an analysis step, not tool
  behavior; the first version of this file misattributed it.

**How RefSeq and GenBank write a known shift.** Counted from RefSeq release 238, 17,648 non-land-plant
records (`code/audit/annotation_census.py`, `code/audit/known_cases.py`):
- **CDS join with 1–2 nt gap or overlap:** 632 records, 685 junctions; 616 at *nad3*. This is also Al Arab's
  definition: "Frameshift are indicated in RefSeq as 'joined' parts of annotated CDS separated by 1 or 2 nt."
- **/ribosomal_slippage:** 6 records (*A. vastus*, *P. marinus*, *P. olseni*, *Rhizophydium*,
  *Daphnia magna*, *Sternotherus*).
- **Notes:** /note mentioning a frameshift in 567 records; notes only, with no join, in 9.
- **Short introns:** intron features ≤ 10 nt in 2 records (*Glarea*, tRNA; *Pyramimonas parkeae*, a 6-nt
  "intron" in *cox1*). The CTGY yeasts in RefSeq are written as joins or notes, not 1-bp introns.
- **Bypass elements:** written as CDS joins with 27–~60 nt gaps (*M. capitatus*, 81).
- **Placement varies.** The *nad3*-174 base is placed at CDS 174 in 341 records and at 175 in 221: the same
  site, so the codon "holding" it is CTC in some records and AGT in others.
- **The fungal set** (BK072095–BK074789) has 0 CDS and 0 intron features in all 2,695 records.

**Database audits** (sources of false frameshifts):
- Sangster & Luksenburg 2021 (K): "Of 1,876 complete or partial mitogenomes of birds ... the authenticity
  of 1,559 could be assessed with sequences of conspecifics. Of these, 78 (5.0%) were found to be
  problematic", 18 of them from sequencing errors or nuclear copies.
- Hassanin 2010 (A): a 1,881-nt nuclear copy inside goat RefSeq NC_005044.
- Kneubehl 2022 (K): in nanopore tick mitogenomes, frameshifts were "almost exclusively single-base
  deletions in long A/T homopolymer regions (> 8 bp)".
- Parson 2013 (K): with Ion Torrent, "About two thirds of those [differences] were observed in or around
  homopolymeric stretches".

## 5. Things that make DNA look frameshifted without a translational shift

1. **Sequencing and assembly errors**, mostly homopolymer indels.
2. **Nuclear copies** inside mitogenome records.
3. **Incomplete stop codons** completed by polyadenylation.
4. **Introns** unannotated or mis-annotated (and frameshifts stored as 1-bp introns).
5. **Insertion/deletion RNA editing:**
   - *Physarum* (now *Badhamia polycephala*, NC_002508, in the scan set; edits not annotated): 54 C
     insertions in *atpA*, 772 further insertions (Mahendran 1991, A; Bundschuh 2011, K);
   - *Didymium* (C);
   - *Diplonema*: ~220 U insertions (Moreira 2016, K);
   - kinetoplastid *coxII*: four non-encoded nucleotides (Benne 1986, A; maxicircles in GenBank, not
     RefSeq);
   - **calcaronean sponges, inside animals** (*Sycon*, *Leucosolenia*; not in RefSeq 238): "single or double
     uridylate (U) insertions in pre-existing short poly(U) tracts" (Lavrov et al. 2016, A; found by the
     red team). These sit in homopolymers, so a homopolymer kill rule would discard real editing as
     sequencing error.
6. **Gene ends**: AGA/AGG stops, out-of-frame read-through. **Gene starts**: *Pagurus* and *Daphnia pulex*
   *cox1* "shifts" at CDS position 4 are start-codon cases.
7. **Translational bypassing** (*Magnusiomyces*): a 27–60 nt insert kept in the RNA and skipped by the
   ribosome. Not an intron and not a 1-nt shift.

## 6. Mechanisms proposed (background for interpreting a hit)

- **Hungry codon plus P-site slippage:** a rare or unread A-site codon stalls the ribosome
  (*Polyrhachis*, turtles; Ty1/Ty3 analogy).
- **Out-of-frame A-site pairing** by an abundant tRNA (glass sponges).
- **Missing tRNAs for unused codons** (*Perkinsus*).
- **Quadruplet reading** by a modified tRNA (yeast tRNA-Cys G37).
- **Bypassing** after a stall at an unused takeoff codon (*Magnusiomyces*).

The thread common to all prize-type cases is that **the trigger codon is almost or entirely absent in
frame elsewhere in the genome.**

## 7. What a new result would have to add

**A positive result** must show all of the following.
1. **A lineage outside every case above:** glass sponges, *Perkinsus*, CTGY yeasts, *Magnusiomyces* byps,
   anomiid bivalves (*Heteranomia*, TAG), the turtle CTN-AGN sites, venerid bivalves, Miliolida
   foraminifera, the chytrid *cob* annotations, and the single sites in section 2.
2. **Pile-up:** shifts at one codon in ≥ 2 genes of one genome, or at several distinct orthologous sites
   across a lineage.
3. **That codon rare in frame.**
4. **Not already written into RefSeq** (632 records) **or the literature.**
5. **Evidence beyond one assembly:** raw reads, a second assembly or run, and RNA keeping the base.

**A clean negative** must beat the earlier scans on scope and sensitivity:
- 17,648 RefSeq records outside land plants plus 2,695 fungi, against Al Arab's 3,842 animal genomes;
- several shifts per gene, and 3′ ends;
- codon and lineage tests;
- per-phylum alignment coverage, showing where a code like *Perkinsus*'s could not have been seen.

**Smaller additions, labeled exploratory:**
- the extent of CTGY shifts across the 2,695 fungi;
- UGG shifts across all glass sponges;
- bypass-like inserts outside *Magnusiomyces*;
- a codon census of the 685 annotated junctions.

## 8. Consequences for the design

See `PRIOR_ART_AUDIT.md`, answer 6.

## Corrections log

- **2026-11-21 19:57–20:10, after the audit:**
  - Added *Magnusiomyces* bypassing (Lang 2014, Nosek 2015, Nosek & Tomáška 2025), missed in round 1
    because those papers say "bypassing" and "jumping", not "frameshift". Found through the Atkins 2016
    review.
  - Added venerid bivalves (Wang 2020), chytrid *cob* annotations, *Arctica*, other RefSeq-annotated single
    sites, Seligmann's claims, and Poehls 2025.
  - Corrected: *Perkinsus* *is* in RefSeq (NC_068742, NC_068743). Round 1 said GenBank only.
  - Corrected: the MitoFinder "internal stop" exclusion was a case-study step, not tool behavior.
  - Corrected: sources read by the sweep agent as keyword windows were labeled "F"; they are now "K-sweep"
    or re-read by me ("K").
  - Corrected: *Physarum* (*Badhamia polycephala*) is in the scan set.
  - Corrected: the CTGY yeasts in RefSeq are written as joins, not 1-bp introns (Wolters' statement
    concerns their own and the GenBank reference formatting).
  - Added MFannot and MITOS2 behavior from their code, and RefSeq annotation counts.
- **2026-11-21 20:25, after the red team** (`PRIOR_ART_AUDIT.md`, answer 7). Every quote below was
  verified by me against the fetched source.
  - Added Lin et al. 2026 (*Heteranomia* TAG +1), which corrects the "nothing from 2025–2026" line.
  - Added Macher et al. 2023 (Miliolida *cox1*). It had appeared in a round-1 result list (Europe PMC query
    19) and was not read.
  - Added Lavrov et al. 2016 (U-insertion editing in calcaronean sponges, inside poly(U) runs).
  - Added Gott 2005, Beargie 2008 and Chen 2012 (editing-aware scans of myxomycete mitogenomes).
  - Widened Seligmann (more than one taxon).
- **2026-11-21 20:11:** the Wolters 2023 entry was rewritten to say exactly what the paper did, with the
  methods passage quoted. It annotated with MFannot and curated by hand, mentions +1C frameshifts in one
  sentence, and searched other species for byps with Lang 2014's method. It does not describe a frameshift
  scan. The summary's "Scans so far" list and the verdict on claim (a) rest on that reading.
- **2026-11-22, new adjacent work (not a frameshift result).** Xu R (AI agents), "Yeasts that read CUU and
  CUA as proline", github.com/ry-xu/genome-research/tree/main/proline-codon-yeasts, first committed
  2026-11-21 20:44 UTC. Read in full (README and notes.md, K). What it is:
  - A codon-meaning census of about 22,700 mitochondrial genomes: 18,015 of 18,031 RefSeq, plus one per
    additional species from GenBank. Translated with declared codes, aligned to profile models, codons
    placed at profile columns. notes.md: "The census was built to look for new stop codons; this came
    out of it."
  - Result: a sense-to-sense reassignment, CUU/CUA as Pro, in *Wickerhamiella siamensis*, *W. infanticola*,
    *W. verensis* and *Sugiyamaella pinicola*. All four are from Sarhan 2025 and are in our fungal set
    (BK074538, BK074532, BK074542, BK074262), read with table 4.
  - Overlap with this project: same data, and it borders our exploratory "possible code change" tally
    (stops only). It does not test for frameshifts or for shifts recurring at a codon, so claims (b) and
    (c) are unchanged.
  - Effect on our scan: CUU/CUA at Pro sites become Leu-for-Pro mismatches, which cost much less than a
    frameshift (F = 15). They should not create shift calls. Any flag in these four will be checked
    against this.
  - The repo's other posted results (cgg-skip-tanb77, ricin-loop-zymomonas, tag-bacteria,
    tta-myxobacteria) are bacterial or not about frameshifts.
- **2026-11-22 08:09, missed prior art found during the scan** (a peek at the code-1 batch with the frozen rules).
  Yang J, Harding T, Kamikawa R, Simpson AGB, Roger AJ (2017) "Mitochondrial genome evolution and a novel
  RNA editing system in deep-branching heteroloboseids", *Genome Biol Evol* 9:1161–1174 (PMID 28453770);
  record NC_034794 (*Eukaryota* sp. BB2), as cited in the record (R).
  - The scan sees it as 57 clean −1 shifts at scattered codons in 13 core genes, the same signature as
    *Badhamia*/*Physarum* insertion editing.
  - It belongs in section 5 (things that make DNA look frameshifted). Heterolobosea editing was not in the
    round-1 or red-team searches.
  - It does not flag: no codon pile-up.

## References

Verification as in the key above; "not saved" means no copy in `lit/`.

- Al Arab M et al. (2017) *Mol Phylogenet Evol* 106:209–216. doi:10.1016/j.ympev.2016.09.024 (F-agent, author manuscript; published version 403)
- Andreu-Sánchez S, Chen W, Stiller J, Zhang G (2021) *GigaScience* 10:giaa161. doi:10.1093/gigascience/giaa161 (F-agent)
- Atkins JF et al. (2016) *NAR* 44:7007. doi:10.1093/nar/gkw530 (K)
- Beckenbach AT, Robson SKA, Crozier RH (2005) *J Mol Evol* 60:141–152. doi:10.1007/s00239-004-0178-5 (A)
- Benne R et al. (1986) *Cell* 46:819. doi:10.1016/0092-8674(86)90063-2 (A)
- Bundschuh R et al. (2011) *NAR*. doi:10.1093/nar/gkr180 (K)
- D'Onorio de Meo P et al. (2012) *NAR* 40:D1168. doi:10.1093/nar/gkr1144 (K)
- Donath A et al. (2019) *NAR*. doi:10.1093/nar/gkz833 (K; MITOS2)
- Beargie C et al. (2008) *Bioinformatics* 24:2571–2578. doi:10.1093/bioinformatics/btn487 (A)
- Chen C et al. (2012) *PLoS Comput Biol*. doi:10.1371/journal.pcbi.1002400 (K)
- Gott JM, Parimi N, Bundschuh R (2005) *NAR* 33:5063. doi:10.1093/nar/gki820 (K)
- Lavrov DV et al. (2016) *Curr Biol* 26:86–92. doi:10.1016/j.cub.2015.11.043 (A)
- Lin Y-T et al. (2026) bioRxiv. doi:10.64898/2026.08.14.744788 (K; full-text XML in `lit/lin_2026.xml`)
- Macher JN et al. (2023) *mBio* 14:e0030223. doi:10.1128/mbio.00302-23 (K)
- Seligmann H (2012) *Biosystems* 110:84–106. doi:10.1016/j.biosystems.2012.09.002 (A)
- Seligmann H, Labra A (2013) *Biosystems* 114:155–163. doi:10.1016/j.biosystems.2013.09.002 (A)
- Forget L et al. (2002) *Mol Biol Evol* 19:310–319. doi:10.1093/oxfordjournals.molbev.a004084 (A)
- Gornik SG et al. (2022) *Mol Biol Evol* 39:msac191. doi:10.1093/molbev/msac191 (F-agent)
- Haen KM, Pett W, Lavrov DV (2014) *Gene* 535:336–344. doi:10.1016/j.gene.2013.10.041 (A)
- Hassanin A et al. (2010) *Mitochondrial DNA*. doi:10.3109/19401736.2010.490583 (A)
- Kahn AS, Geller JB (2018) *Genome Announc*. doi:10.1128/genomeA.00234-18 (K)
- Kayal E, Smith DR (2021) *Mol Biol Evol*. doi:10.1093/molbev/msab041 (K)
- Kneubehl AR et al. (2022) *Sci Rep*. doi:10.1038/s41598-022-23393-5 (K)
- Krüger A, Kovalchuk D, Shiriaev D, Rorbach J (2024) *Hum Mol Genet* 33:R42–R46. doi:10.1093/hmg/ddae032 (F-agent)
- Lang BF et al. (2014) *PNAS* 111:5926. doi:10.1073/pnas.1322190111 (A; Europe PMC full text failed)
- Lang BF et al. (2023) *Front Plant Sci*. doi:10.3389/fpls.2023.1222186 (K; MFannot)
- Langkjaer RB et al. (2003) *NAR*. doi:10.1093/nar/gkg423 (A; full text failed)
- Mahendran R, Spottswood MR, Miller DL (1991) *Nature* 349:434. doi:10.1038/349434a0 (A)
- Masuda I, Matsuzaki M, Kita K (2010) *NAR* 38:6186–6194. doi:10.1093/nar/gkq449 (A)
- Milbury CA, Gaffney PM (2005) *Mar Biotechnol* 7:697–712. doi:10.1007/s10126-005-0004-0 (A)
- Mindell DP et al. (1998) *Mol Biol Evol* 15:1568 (C)
- Moreira S et al. (2016) *NAR*. doi:10.1093/nar/gkw188 (K)
- Nosek J, Tomaska L, Burger G, Lang BF (2015) *Trends Genet* 31:187–194. doi:10.1016/j.tig.2015.02.010 (A)
- Nosek J, Tomáška Ľ (2025) *Cold Spring Harb Perspect Biol*. doi:10.1101/cshperspect.a041849 (A)
- Oudot-Le Secq MP, Green BR (2011) *Gene*. doi:10.1016/j.gene.2011.02.001 (A)
- Parson W et al. (2013) *Forensic Sci Int Genet*. doi:10.1016/j.fsigen.2013.06.003 (K)
- Poehls J et al. (2025) bioRxiv. doi:10.1101/2025.10.02.679949 (A)
- Rosengarten RD et al. (2008) *BMC Genomics* 9:33. doi:10.1186/1471-2164-9-33 (F-agent)
- Russell RD, Beckenbach AT (2008) *J Mol Evol* 67:682–695. doi:10.1007/s00239-008-9179-0 (F-agent)
- Sangster G, Luksenburg JA (2021) *Genome Biol Evol*. doi:10.1093/gbe/evab210 (K)
- Sarhan MS et al. (2025) *Sci Data*. doi:10.1038/s41597-025-06447-x (K; the data source)
- Seligmann H, Warthi G (2019) *Comput Struct Biotechnol J*. doi:10.1016/j.csbj.2019.08.006 (A)
- Szabóová D, Hapala I, Sulo P (2018) *FEMS Yeast Res* 18:foy028. doi:10.1093/femsyr/foy028 (A)
- Temperley R et al. (2010) *Science* 327:301. doi:10.1126/science.1180674 (A)
- Towett-Kirui S, Morrow JL, Riegler M (2022) *Sci Rep* 12:477. doi:10.1038/s41598-021-04398-y (F-agent)
- Uliano-Silva M et al. (2023) *BMC Bioinformatics*. doi:10.1186/s12859-023-05385-y (K; MitoHiFi)
- Allio R et al. (2020) *Mol Ecol Resour*. doi:10.1111/1755-0998.13160 (K; MitoFinder)
- Wakigawa T et al. (2026) *Nat Commun* 17:6128. doi:10.1038/s41467-026-75248-6 (F-agent)
- Wang Y, Yang Y, Liu H, Kong L (2020) *Zool Scr*. doi:10.1111/zsc.12454 (A, via OpenAlex)
- Wolters JF et al. (2023) *Front Microbiol* 14:1268944. doi:10.3389/fmicb.2023.1268944 (F-agent)
- Wú F et al. (2024) bioRxiv. doi:10.1101/2024.05.25.595893 (A)
- Yang J, Harding T, Kamikawa R, Simpson AGB, Roger AJ (2017) *Genome Biol Evol* 9:1161–1174. PMID 28453770 (R; record citation)
- Yuan B et al. (2025) *Genes* 16:1030. doi:10.3390/genes16091030 (K)
- Zhang H et al. (2011) *PLoS ONE*. doi:10.1371/journal.pone.0019933 (K)
- Zhang et al. (2016) *Mitochondrial DNA A*. doi:10.3109/19401736.2014.945552 (A)
