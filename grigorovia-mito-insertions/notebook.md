# Notes: a scan for codon-specific frameshifts in mitochondria

The full write-up.
- **Numbers.** Every number here comes from a file in `protocol/` or `candidate/`, named in the text.
- **Rules and changes.** The rules, predictions and every change are in `protocol/frozen_rules.md`, with times and commit hashes
  from a local git repository. These were not externally timestamped.
- **Novelty.** "As far as we can find" applies to every statement that something is new.
- **Review.** Changes made after the fact-check and the referee review are listed in `protocol/corrections.md`.

**Approach.** It follows [ry-xu/genome-research](https://github.com/ry-xu/genome-research) (Ryan Xu): a
frameshift-aware scan, rules fixed before the data, controls with numerical bars, then reads, relatives and a blinded
re-check. The code is our own.

## 1. Question

Do any mitochondria shift reading frame at one recurring codon, as glass sponges do at UGG? Known cases
(`protocol/prior_art.md`, section 1):

- glass sponges, +1 at UGG;
- *Perkinsus*, +1 at AGG and +2 at other motifs;
- *Kazachstania*-group yeasts, +1 at CTG;
- *Heteranomia* (anomiid bivalves), +1 at in-frame TAG (Lin et al. 2026, bioRxiv);
- a recurring CTN-AGN motif in turtles;
- translational bypassing of about 50-nt inserts in *Magnusiomyces*, a different mechanism.

Single programmed sites, such as *nad3* in birds and turtles, are not the target. They were used as controls.

## 2. Data

**RefSeq mitochondrion release 238: 18,220 records.**
- The 17,648 that are not land plants form the main scan.
- A separate stratum of land plants was planned but not run.

**The 2,695 fungal mitogenomes of Sarhan et al. 2025** (GenBank TPA BK072095–BK074789).
- The deposited records carry no gene features.
- The authors annotated them with MFannot, but that annotation is not public (`protocol/design.md`, amendment 1).

**Genetic code.** Each genome is translated with its own mitochondrial code: the record's CDS table, or NCBI
Taxonomy's mitochondrial code where there are no CDS.

**Composition of the 20,343 scanned genomes** (`protocol/scan/scan_composition.tsv`):

| group | genomes |
|---|---|
| animals | 16,075 |
| fungi | 3,582 |
| protists and algae | 686 |

**Accession lists:** `protocol/accessions.tsv`.

**Completeness:** all 20,343 were scored (`protocol/scan/completeness.md`).

## 3. Method (the frozen rules; `protocol/frozen_rules.md` section 1)

**Alignment.**
- DIAMOND 2.2.9 blastx in frameshift mode (`--sensitive -F 15 --range-culling -k 10 -e 1e-5`), each genome under its
  own code.
- References: 53,641 RefSeq mitochondrial proteins, one genome per family. NC_079576 was removed as a mislabeled
  record.
- Hits to the query's own family are dropped, then the hits are culled to the best per locus.

**Shift events.** Frameshift operations inside one alignment, folded to +1 or −1.

**Clean shifts.** A shift is clean unless any of these applies:
- paired with an opposite shift within 150 nt;
- within 30 nt of an inferred intron or insert edge;
- at a gene start or end;
- short of 15 aligned columns at ≥ 20% identity on either side;
- near a run of one base (≥ 5 for +1, ≥ 4 for −1).

**Codon.**
- Each shift has a placement window: all positions within ±9 nt that give the same protein.
- It votes for the old-frame codon at each end of the window, if the window is ≤ 4 nt wide (a cap set from the planted
  data).

**Expectation.**
- Random insertions and deletions were planted into the core genes of 620 genomes, three copies each, and run through
  the same pipeline.
- Their clean shifts give each codon's vote rate, by genetic code and GC half.

**Single-genome test** (confirmatory).
- In core genes, ≥ 2 clean shifts of one direction vote for the same codon, in ≥ 2 genes.
- The binomial tail against the planted rate must be ≤ 1e-3 for sense codons, ≤ 3.2e-3 for stop codons.

**Pool test.** Family and order pools. Sites count once if within 50 residues. Flag at ≤ 2.5e-6.

**Strong flag.**
- Frozen rule: expected by chance less than once across the whole scan (E < 1), or the same codon flagged in a related
  genome.
- Added after the freeze, a conservative change: records from one assembly or strain count as one genome.
- Only strong flags outside the listed known lineages go to the funnel (reads, relatives and a blinded re-check), and
  the funnel is not started without the user. If there are none, the result is a short negative.
- The known-lineage list works by family and order names, not by codon.

## 4. Controls and missed bars

**First run** (commit 972ae02; numbers in the `protocol/frozen_rules.md` log). Every bar below was missed, and no α met the
threshold rule:

| bar | result | target |
|---|---|---|
| two shifts planted at one codon flagged | 0.28 of trials | 0.90 |
| RefSeq-annotated frameshifts detected, *nad3* | 0.756 | 0.80 |
| RefSeq-annotated frameshifts detected, other | 0.788 | 0.80 |
| annotated fungal introns inferred | 0.62 | 0.80 |
| shifts planted near intron edges recovered | 0.28 | 0.80 |

**Decisions before the single allowed re-run** (`protocol/design.md` amendment 4; `protocol/frozen_rules.md` log).
- one codon rule, with a window cap;
- the two-site bar scored when both sites are counted, with the unconditional rate as the headline;
- α fixed at the values the pre-written rule gave;
- the fungal intron blind spot stated as a limit;
- a strong-flag rule.

**Two changes failed their own tests and were dropped.**
- A reference-length filter raised near-edge coverage from 0.438 to 0.439. The rule logged beforehand needed +0.05.
- Frame changes between two hits: no planted shift was recovered only that way, so they became a tally only.

**Re-run** (`protocol/controls/BARS.md`).

- **Gates for starting the scan: all passed.**
  - 0 flags among 360 null sequences: 120 genomes, plus copies with homopolymer and random sequencing-error indels.
  - False clean shifts at intron edges: 0.0018 (n = 17,742; bar ≤ 0.02).
  - All five named positives flag: *Aphrocallistes vastus* (TGG), *Perkinsus marinus* and *P. olseni* (AGG),
    *Huiozyma sinensis* and *Monosporozyma servazzii* (CTG).
  - *Magnusiomyces* shows 65 bypass-like inserts.
- **Bars passed.**
  - Planted single shifts: clean recall 0.901.
  - Two sites at one codon, when both are counted: flagged 1.000 (n = 76; stop codons 0.976, n = 42).
- **Bar missed, kept as a limit: annotated RefSeq frameshifts.**
  - Detected at the site: 0.756 for *nad3* (n = 619) and 0.788 for others (n = 66).
  - Detected anywhere in the gene: 0.989 and 0.894. Counted at the right codon: 0.590 and 0.409.
  - The 155 sites found in the gene but not at the site are placed 10–150 nt off in poorly conserved stretches
    (median 75 nt; 153 of 155 within 150 nt).
  - The rate difference implies about 151 such sites. Four sites detected at the site have no matching CDS and so count
    as not "in gene".
  - My first diagnosis, that genes were split between two hits, was wrong. It was corrected in the log before the
    re-run.
- **Stated limits: fungal introns.**
  - Inferred introns 0.619.
  - Recall of shifts planted near intron edges 0.279.
  - Per-genome core coverage, median 0.795.

## 5. Results

**Flags** (`protocol/scan/flags_summary.md`, `protocol/scan/flag_context.md`).
- 143 genome-directions have ≥ 2 clean core shifts.
- The expected number of chance flags under the rule is 2.5, an upper bound.
- There were 16 flags:

| flag | genomes | status |
|---|---|---|
| +1 TGG | *Aphrocallistes vastus* | listed (glass sponges) |
| +1 AGG | *Perkinsus marinus*, *P. olseni* | listed |
| +1 CTG | *Huiozyma sinensis*, *Monosporozyma servazzii*; unannotated records of *H. naganishii*, *Jamesozyma lodderae*, *Kazachstania viticola*, *M. aerobia*, *M. solicola* | listed (CTG-shifting yeasts; family Saccharomycetaceae) |
| +1 AGT | *Podiceps cristatus*, *Tachybaptus novaehollandiae* (grebes); *Graptemys ouachitensis*, *Malaclemys terrapin* (turtles) | listed |
| +1 AGG | *Grigorovia jiainica* (nom. inval.) | listed only through its family; see section 8 |
| −1 CAC | *Verticillium dahliae* VdLs.17 | **strong, outside listed lineages; killed at the second-assembly step** |

**The bird and turtle flags.** Each is *nad3*-174 plus a second documented AGT site.
- Turtle *nad4l* follows the CTN-AGN motif (Russell & Beckenbach 2008).
- Grebe *nad1* is annotated in one RefSeq record (*T. novaehollandiae*), with no paper in our prior-art file. It is not
  annotated in *P. cristatus*, though the scan finds it there too.

***K. viticola*.** Its count of 6 includes one duplicated event, so it has 5 distinct sites. This is the only duplicate
in the scan (`protocol/scan/scan_report.md`).

***Verticillium dahliae* VdLs.17** (a 2011 Broad genome-project assembly in three scaffolds).
- The flagged scaffold has 29 clean core shifts, all −1. Its chance of flagging somewhere, given 29 random shifts, is
  0.021.
- The second-assembly step was pre-approved and ran only because the flag was strong. It tested all 38 clean shifts on
  the three scaffolds: 31, 6 and 1, including 2 in non-core genes.
- *V. dahliae* NC_008248 (Pantou et al. 2006) carries the base VdLs.17 lacks at all 38 positions. *V. nonalfalfae*
  NC_029238 carries it at 37; the 38th is not comparable (`protocol/scan/second_assembly_VdLs17.tsv`).
- So the flag is an assembly error.

**Pool flags.** Saccharomycetaceae and Saccharomycetales, +1 CTG, 6 of 9 sites (listed).
- The shuffled-pool control flagged 2, 3 and 4 family pools, against the design's "≤ 1".
- Each of those came from one CTG-yeast genome with ≥ 3 CTG sites landing in a random pool.
- In this scan the pool test only repeats single-genome signals.

**Exploratory tallies** (descriptive, never flags).
- **Many-shift tally** (`protocol/scan/many_shift_tally.md`): 23 genomes have clean shifts in ≥ 3 genes.
  - The top two are known RNA editors: *Eukaryota* sp. BB2 (Yang et al. 2017) and *Badhamia* (*Physarum*).
  - Three tilapia (*Oreochromis*) records share a cluster of +1 shifts in *nad2*.
  - The mark "possible editing or noisy assembly", as worded, also catches the known CTG-shifting yeasts.
- **Paired breaks recurring at one codon** (`protocol/scan/paired_break_tally.md`): 37 genome-codon rows.
  - Both *Perkinsus* appear (AGG), so the tally catches the known AGG case.
  - The rest are mostly AT-rich codons in arthropods, nematodes, cnidarians and others, usually at low alignment
    identity. Some rows are other codons at identity up to 0.93.

## 6. Sensitivity of a negative result

**Planted** (`protocol/controls/BARS.md`).
- A genome with 2, 3 or 5 planted same-codon shifts is flagged 0.317, 0.483 or 0.717 of the time (sense codons).
- A single random planted shift is counted at its codon 0.553 of the time: animals 0.564, fungi 0.528, protists 0.553.
- Losses, as a share of planted shifts:

| reason | share |
|---|---|
| homopolymer filter | 0.149 |
| placed at another codon | 0.100 |
| not aligned | 0.092 |
| flank filter | 0.045 |
| aligned but no shift written | 0.041 |

Smaller categories make up the rest (`protocol/controls/BARS.md`).

**Real documented sites** (`protocol/scan/lineage_sensitivity.md`; the definitions were added after the freeze).
- **Headline 0.29:** the unweighted mean of 8 per-lineage rates. Bootstrap over lineages 0.09–0.54; pooled 0.32 over 84
  sites.
- ***nad3*-174 is excluded and given separately: 0.61 over 598 sites.** It is one homologous position that was gained
  and lost many times (Andreu-Sánchez et al. 2021).
- **Implied flag chance** for a genome with 2, 3 or 5 documented same-codon sites: 0.08, 0.20 and 0.44. Ranges:
  0.01–0.29, 0.02–0.56, 0.07–0.86.
- **Out of sample:** 7 sites set no rule (not annotated and not in a named positive). 3 of them were detected as clean
  shifts; their codon could not be checked.
- **By lineage, sites counted at their own codon:**
  - *Perkinsus* 12 of 38; the pair filter removed 14, because +1 and +2 shifts sit close together;
  - CTG-shifting yeasts 4 of 10. The summary table lists 11 sites; the 11th is a gene-level site whose codon cannot be checked;
  - birds 2 of 14 and turtles 7 of 16, beyond *nad3*-174;
  - the remaining lineages are in `protocol/scan/lineage_sensitivity.md`.

**Scope.** Most scanned genomes are animals; 686 are protists or algae. Only core genes were tested, recall near fungal
introns is 0.28, and listed lineages were not tested for new codons.

**Limits.**
- A negative rules out lineages only at these real-data rates, and only within this sample.
- The pair filter removes real sites in lineages that mix +1 and +2 shifts.
- The chance figures assume errors land at random codons. Assembly errors cluster, so a noisy genome can be strong
  without being real, as VdLs.17 was.

## 7. Predictions (written before the scan; `protocol/scan/predictions.md`)

| | | result | count |
|---|---|---|---|
| P1 | named positives flag | pass | 5/5 |
| P2 | another glass sponge flags at +1 TGG | **fail** | 0/3 |
| P3 | another CTG-shifting yeast flags at +1 CTG | pass | 4 of 16 in the four named genera (and *Jamesozyma lodderae*) |
| P4 | no bird or turtle flag involves *nad3* | **fail** | 4 flags |
| P5 | flags outside listed lineages ≤ chance + 5 | pass | 1 ≤ 2.5 + 5 |
| P6 | *Magnusiomyces* in the bypass tally | pass | 65 inserts |

**Why P2 failed** (`protocol/scan/removed_breaks_known_groups.md`). None of the three other glass sponges has a break at
TGG, so no filter removed one; the sites are absent or unknown.
- *Oopsacas minuta*: its own genome paper reports no +1 frameshifting (Jourda et al. 2015).
- *Caulophacus*: no documented sites.
- *Vazella pourtalesii*: unknown. Its source paper is Haen et al. 2014, which reports insertions in its new genomes, but
  we could read only the abstract, and that names no species or positions.

**Why P4 failed.** The prediction was wrong. I predicted that birds and turtles could not flag because they have one
site each, but some have a documented second AGT site.

## 8. *Grigorovia jiainica* (exploratory; followed up after the result)

**How it was chosen.**
- The frozen rule counted this flag as a listed lineage, because its family, Saccharomycetaceae, is on the list for the
  CTG-shifting yeasts.
- Its codon, AGG, is not the CTG mechanism. The user sent it to the funnel after seeing the result, and logged that it
  was a post-result decision.
- Everything in this section is exploratory.
- Steps a, b and f had stop rules written before they ran. The others were descriptive (`protocol/frozen_rules.md`; files in
  `candidate/funnel/`).

**The name.**
- Described as *Kazachstania jiainicus* (Lee et al. 2008).
- Placed in *Grigorovia* in 2020 (NCBI: "nom. inval.").
- Combined as *Pachytichospora jiainica* in a 2023 preprint (Malimas et al. 2023).
- The strain sequenced is NRRL Y-48843.

### 8a. Sites, reads and the blinded re-check

**Sites** (`candidate/funnel/step_a_sites.md`). Four clean +1 shifts:
- *cox1* 24711;
- *cob* 26496, 26477 and 26266.

In each, the extra base can be placed just before a GGC codon:
- at *cob* 26496 (AGG) and 26477 (TGG), the placement window is a single position;
- at *cox1* (GGA·AGG) and *cob* 26266 (ACT·AGG), a neighboring placement gives the same protein.

Read straight on, each gene meets a stop 3 to 9 codons later. Read with the extra bases removed, *cox1* runs 450
codons to the end of its alignment without a stop. *cob* runs on for 6, 70 and 172 codons between and after its three
sites.

**Reads** (`candidate/funnel/step_f_reads.md`). Rule fixed before the count: ≥ 10 spanning reads, ≥ 90% with the extra base.
- **Passed at all four sites: 1,440, 1,718, 1,733 and 1,785 reads with the base, none without.**
- The run is SRR16974510, the Illumina run the assembly was built from. That excludes assembly and read errors in
  this sample. It is one strain, one culture and one library, and read support from the same run is necessary but not
  sufficient (`protocol/design.md`, decision 13).
- Other deposited strains (CBS 10587, ATCC MYA-4454, TL1903) have marker sequences only.

**Blinded re-check** (`candidate/funnel/blind/`).
- **Setup.**
  - Three coded genomes: *G. jiainica*, its closest relative (*Kazachstania bromeliacearum*, no shifts in the scan) and
    a CTG-shifting yeast (*M. solicola*).
  - The agent used miniprot and the core proteins of *S. cerevisiae* and *Saccharomycodes ludwigii*.
  - **Both protein sets are also in the scan's references.** The scan called the four *G. jiainica* sites against the
    *S. ludwigii* proteins. Only the aligner was independent.
  - The agent was asked to find and describe frameshifts, give the reference residue after each, and say whether they
    looked systematic.
- **Result.**
  - The same four sites in *G. jiainica*, with glycine at each under its preferred placement. At 26266 miniprot's own
    placement faces Thr.
  - None in *K. bromeliacearum*.
  - Five +1 CTG sites, read as Cys, in *M. solicola*, which has no documented sites; the scan counted 3.
- **What it shows.** Another aligner finds the same bases. The check cannot tell the possible causes apart.

### 8b. Evidence on the cause before the follow-up tests

**Fits a translational skip.**
- The pattern: one extra base, each time just before GGC.
- Glycine in distant reference proteins at the GGC position (`candidate/funnel/step_b_residue.md`): 8 of 10, 10 of 10, 10 of 10. At
  *cob* 26266 it is 6 of 9, under the AGG placement; the scan's own placement faces Thr.
- The GGC is derived: in the two relatives the blinded check compared, these positions are GGA or GGT.

**Does not fit, or fits other causes as well.**
- **Respiration** (`candidate/search_log/grigorovia_literature.md`).
  - The species description lists only glucose and galactose as assimilated: glycerol, ethanol, DL-lactate, succinate
    and citrate are negative (Lee et al. 2008).
  - A trait table agrees; it is curated from *The Yeasts* (2011) (Opulente et al. 2018).
  - The congener *G. yakushimaensis* is also negative there.
  - One stock of the congener *G. transvaalensis*'s type strain was reported not to grow on non-fermentable carbon
    (Fekete et al. 2007), while Y1000+ records slow glycerol growth for CBS 2186 (Opulente et al. 2024). The public
    reads of NRRL Y-17245 contain almost no mtDNA (§8c).
  - These records were the basis of the "disabling mutations" reading. §8c tests it.
- **Glycine decoding** (`candidate/funnel/trna_gly_structure.md`, `candidate/funnel/agg_usage_relatives.md`).
  - Its only mitochondrial tRNA-Gly has anticodon GCC (1,329 reads; 4 for TCC). By standard pairing it cannot read the
    38 GGA codons in its aligned core genes.
  - The relatives' tRNA-Gly is UCC.
  - The anticodon loop is an ordinary 7 nt with A37, so nothing like the G37 proposed for the CTG-shifting yeasts.
- **Introns** (`candidate/funnel/introns_and_congener.md`).
  - In *S. cerevisiae* coordinates, the extra bases sit 7 nt from intron bI1 and 3 nt from bI2 in *cob*, and 12 nt
    from aI3 in *cox1*. The fourth is more than 100 nt from any *S. cerevisiae* intron.
  - The congener *G. yakushimaensis* has an intron 2–3 nt from the TGG site (§8c).
  - Imprecise intron loss is one way to leave single extra bases. §8c maps introns across the genus.
- **The congener** (`candidate/heldout/summary.md`, `candidate/funnel/introns_and_congener.md`).
  - *G. yakushimaensis* (OR703946) lacks the extra base at the three comparable positions. §8c compares the fourth
    across its spliced gene, where it is absent too. The written prediction,
    that it would share sites if the trait ran through the lineage, failed.
  - At the *cox1* site it has Ser (AGT), not Gly. Removing the G of AGG instead of the A gives Ser in *G. jiainica*
    too.
- **The six closest relatives** in the scan (75.1–76.2% protein identity) lack all four extra bases: 22 of 24
  comparisons lack them, 2 are not comparable. So the sites are probably recent.
- **AGG use is not informative.** *G. jiainica* has 0 AGG in 1,748 aligned core codons, but five of its six closest
  relatives also have 0. An extra A before GGC creates an AGG anyway.

**Literature** (`candidate/search_log/grigorovia_literature.md`).
- No report of these sites under any name of the species, and no published +1 at AGG or TGG in a yeast mitochondrion.
- In a fungal mitochondrion, the nearest case is a database annotation in the chytrid *Rhizophydium* *cob*. It is a
  different event: the C after AGG is skipped and AGG stays Arg.
- AGG followed by C is a +1 site in the yeast Ty1 retrotransposon, as CUU-AGG-C (Belcourt & Farabaugh 1990;
  Farabaugh et al. 2006). None of the four sites has the CUU.

**Other public data** (`candidate/funnel/completeness_and_data.md`).
- SRA has DNA for all four *Grigorovia* species, and no RNA.
- Of the three congeners, only *G. yakushimaensis* had a mitogenome assembly. We assembled the other two from reads (§8c).
- Proteomics repositories were not searched.

**Reading before the follow-up tests.**
- The extra bases are real in the sequenced culture, and they are patterned.
- Three readings remained: a +1 translational skip, removal by RNA editing, or disabling insertions in genes no
  longer needed.
- The respiration records, the glycine-tRNA problem, the intron positions and the congener's Ser at *cox1* made the
  last look at least as likely as the first.
- §8c tests this.

### 8c. Six follow-up tests of bypass against loss of function (exploratory; `candidate/followup/`)

Two readings were weighed:
- **H1 (bypassed):** *cox1* and *cob* are translated, and the extra bases are read past.
- **H2 (broken):** the yeast no longer respires; the genes are dispensable and the extra bases disable them, possibly
  as scars of intron loss.

What each predicted, and which result would favor which, was logged in `protocol/frozen_rules.md` before each test ran. Some
inputs had already been seen: the congener's codons at *cox1* and *cob* 26266, and the *S. cerevisiae* intron
distances. The H1 prediction for test 6 completes a truncated instruction and was written by the agent. Methods changes
made while running are logged there too, as are the checks added after the referee review. The summary is
`candidate/grigorovia_jiainica.md`.

| test | H1 predicted | H2 predicted | result | favors |
|---|---|---|---|---|
| 1. Congener codons (`candidate/followup/test1_*`) | Gly in the congener at every comparable site | no constraint | *cob* 26496: Gly, matching the A-skip. *cob* 26266: Thr-Gly, as under every top placement. *cox1*: *G. yakushimaensis* has Ser, which skipping a G reproduces exactly. *cob* 26477: not comparable (congener intron) | H2 at *cox1* by the logged criterion. But that Ser is derived: the other two congeners and the three outgroups have Gly (`candidate/followup/genus/ggc_sites.md`, added after review) |
| 2. Introns (`candidate/followup/test2_*`) | no prediction | sites cluster near intron positions (P < 0.05), within 2 nt of junctions of introns *G. jiainica* lost | 3 of 4 within 12 nt of *S. cerevisiae* intron positions (P = 0.0042), but 44% of coding positions lie within 12 nt of some budding-yeast intron (P of ≥ 3 of 4 = 0.23). The nearest is the congener's bI2-type intron, 2 or 3 nt from 26477 (tied; 3 under the group I splice rule); the others are 7–20 nt away or near none. *G. transvaalensis* lacks that intron and has no extra base there. No intron-ORF remnants | no precise scar; a weak link to intron positions (the *S. cerevisiae* null was met) |
| 3. Nuclear genes (`candidate/followup/test3_*`) | all intact | losses, starting with the *cox1*/*cob* activators | no broken copy of any callable gene, including the activators PET309, MSS51, CBS2 and CBP1. Each activator's best hit in the *S. cerevisiae* proteome is the expected gene, except PET122, which finds none (`candidate/followup/nuclear/rbh.md`, added after review) | H1, weakly (a recent loss would not show yet) |
| 4. Congener phenotypes (`candidate/followup/test4_*`; Mikata et al. 2001; Opulente et al. 2018, 2024) | — | congeners respire while *G. jiainica* does not | *G. yakushimaensis* and *G. humatica* are also negative on non-fermentable carbon, yet keep intact genes. *G. transvaalensis* is variable, and its type stock has lost its mtDNA | undercuts H2's premise |
| 5. The rest of the mitogenome (`candidate/followup/test5_*`) | intact | further damage, especially in *cox2*/*cox3* | protein genes full length without frameshifts or stops; *VAR1* open; rRNAs and RNase P RNA present (presence only); tRNA set identical to the congener's except tRNA-Gly | H1, weakly |
| 6. The genus (`candidate/followup/test6_genus.md`, `candidate/followup/genus/`) | AGG skips with a GCC tRNA-Gly; purifying selection on the *G. jiainica* branch | extra bases only in *G. jiainica*; no tie to the tRNA; relaxed selection | extra bases only in *G. jiainica* (12 of 12 congener comparisons lack them); tRNA-Gly(GCC) only in *G. jiainica*. That outcome was near-certain given 100 relatives with UCC, and *G. jiainica*'s GCC was known in advance. *cox1* + *cob* ω = 0.031 (background 0.033; congeners 0.013–0.044); ω = 1 rejected (p = 1.4 × 10⁻⁷). The 95% upper bound limits a neutral period to 1.2–5.3% of the branch across codon models and outgroup choices (2.5% under F3x4) (`candidate/followup/genus/dnds_sensitivity.md`, added after review) | H1, or H2 only if very recent |

**Checks added after the referee review** (`candidate/followup/genus/dnds_sensitivity.md`, `candidate/followup/genus/ggc_sites.md`, `candidate/followup/reads_strand.md`,
`candidate/followup/nuclear/rbh.md`).
- **Codon models and outgroups.** ω = 1 is rejected for *cox1* + *cob* under F1x4, F3x4 and F61, and without the
  outgroups. The largest neutral share is 1.2% (no outgroups), 2.5% (F3x4), 4.9% (F61) or 5.3% (F1x4); for *cox1*
  alone, 5–29%.
  - The *cob*-specific results vary with the model: its dS on the branch ranges 0.71–2.05, and the p-value for "more
    constrained than background" ranges 0.022–0.46.
  - The share is measured in synonymous change. N·dS on the branch is about 2,000, so each 1% of the branch adds only
    about 20 expected nonsynonymous changes. A loss in the last fraction of a percent of the branch, or in culture,
    would add too few to detect.
- **GGC.** *G. jiainica*'s seven core genes have 4 GGC codons in 1,815 (5.4% G or C at third codon positions), all at
  the extra bases. The congeners and outgroups have GGT, GGA or GGG there, and Gly is ancestral at all four positions.
- **Reads.**
  - The extra bases are on both strands in similar numbers (e.g. *cox1* 639 forward, 801 reverse).
  - 1,516 reads span both close *cob* sites with both bases, and none lack either.
- **Orthology.** Searched against all 6,733 reviewed *S. cerevisiae* proteins, the *G. jiainica* and congener
  activators find their expected genes as best hits. The exception is PET122, which finds none. Synteny was not
  checked.

**Reading after the tests.**
- An old loss of function is ruled out: *cox1* and *cob* were constrained for nearly all of the branch, and none of
  the other nuclear or mitochondrial genes we checked is broken.
- These tests cannot tell a bypass (a translational skip, or RNA editing) from a loss within roughly the last 1–5% of
  the branch, the range of the bound: both predict the same results. A loss in the last fraction of a percent, or in
  culture, would leave too few changes to detect at all.
- *Superseded on 2026-11-23: the cause is now stated as unresolved (see `protocol/corrections.md`, "Changes after
  an outside review"). The original reading follows.* The lean toward a bypass rests on the pattern of the
  insertions:
  - each is an extra A or T beside a new GGC;
  - the four GGCs are the only GGC codons in the core genes;
  - the only tRNA-Gly is GCC.
- That pattern was noticed after the case was chosen.
- DNA cannot settle it; RNA, peptides, a respiration test or a PCR of an independent stock of the type strain could.

**Methods in brief** (`code/explore/grigorovia_h1h2.py`, `nuclear_genes.py`, `grigorovia_genus.py`,
`read_haplotypes.py`).
- **Introns.** RefSeq annotation, a single-intron spliced nucleotide alignment, and miniprot 0.18 exon structures.
- **Nuclear genes.**
  - Searched with miniprot and DIAMOND 2.2.9 (ultra-sensitive, frameshift-aware), using *S. cerevisiae*, relatives'
    and congener proteins as queries.
  - Completeness by compleasm 0.2.10 (BUSCO saccharomycetes_odb12, one reference protein per gene; `candidate/followup/completeness/`).
  - Read checks with minimap2 2.31 and samtools 1.25.
  - A local assembly with SPAdes 4.3.0.
- **Mitogenomes.**
  - SPAdes 4.3.0 on WGS reads: *G. humatica*, all reads; *G. transvaalensis* CLIB 448, then two rounds of baiting
    mitochondrial reads and reassembling.
  - Confirmed by mapping all reads back (minimap2 2.31; 100% of the main contigs' positions at ≥ 10×).
- **Codon alignments.** MAFFT 7.526.
- **dN/dS.** codeml (PAML 4.10.10), yeast mitochondrial code, branch models with each *Grigorovia* lineage as
  foreground, and a profile-likelihood bound.
- **tRNAs.** tRNAscan-SE 2.0.12 and ARAGORN 1.2.41.
- **Intron type.** Rfam group I/II models (Infernal 1.1.5): no significant hit.
- **After review.** codeml under F1x4/F61 and without outgroups; reciprocal best hits (DIAMOND blastp, UniProt
  *S. cerevisiae* reference proteome); read counts by strand.

**False alarms in test 3** (`candidate/followup/nuclear/`).
- **First pass (the logged rule: miniprot with *S. cerevisiae* queries).** It called PET309, CBP1 and CBS2 absent in
  all four *Grigorovia*. Their proteins are only 24–32% identical to the relatives'.
  Relatives' and congener proteins with DIAMOND found all three in every species.
- **Congener queries.** These called COR1 and PET309 "broken" (one frameshift each) and COX9, PET111 and PET122
  "absent" in *G. jiainica*. In the read check:
  - the two frameshifts disappear in the read consensus;
  - COX9 and PET111 are in the assembly but beyond miniprot at 31–65% identity;
  - PET122 sits in an assembly gap, and the 203 read pairs that hit it assemble into an intact copy.
- **Heterozygosity scan.** A per-position pileup view suggested heterozygous frameshifts or a stop-gain in COR1, PET111
  and PET309.
  - Grouping the reads that span each spot by their exact sequence showed both copies open and in frame at every spot.
  - At the two assembly frameshifts, no spanning read matched the assembly: 0 of 25 at COR1, 0 of 23 at the PET309
    repeat.

**A caution for users of *G. transvaalensis* NRRL Y-17245 (= CBS 2186, the type strain).**
- Its public WGS reads contain essentially no mtDNA: 4 of 200,000 reads match mitochondrial proteins, against 991–2,933
  in the other strains (`candidate/followup/genus/mtreads.md`). No mitochondrial contig assembles (`candidate/followup/genus/pick.md`).
- SRR16974369 and SRR6475934 are the same dataset, listed twice.
- This matches Fekete et al. (2007): their collection stock already did not grow on non-fermentable carbon.
- Y1000+ records slow glycerol growth for CBS 2186, while the Fekete stock did not grow, so the stocks differ.
- Phenotypes recorded from the mtDNA-less stock describe that culture, not the species. Strain CLIB 448 has mtDNA
  (37 kb assembled, incomplete).

**Ploidy** (`candidate/followup/nuclear/heterozygosity.md`).
- The *G. jiainica* nuclear genome is heterozygous: 11.95 sites per kb have a second base in 25–75% of reads, with
  allele fractions peaking at 0.4–0.5. That fits a diploid or higher ploidy.
- Heterozygous sites in every panel gene were therefore screened, and each flagged spot was resolved from the reads.
- The mitogenome is unaffected: each extra base is in 100% of mitochondrial reads.

**Downloads** (`candidate/followup/genus/downloads.md`).
- 1.61 GB of public reads: *G. humatica* SRR16974597, *G. transvaalensis* SRR16974369 and ERR9706995.
- Six nuclear assemblies (0.02 GB) and the BUSCO lineage data (0.25 GB).
- Under the 20 GB limit set for the project.

## 9. Changes after the freeze (logged with times in `protocol/frozen_rules.md`)

1. **The *Grigorovia* override.** Sent to the funnel after the result; stated as such.
2. **Descriptive context for every flag:** sites, shift totals and direction mix, platform and assembly, related
   genomes at the same codon, and second assemblies. The prediction-scoring criteria were fixed before the results.
3. **The many-shift tally:** a full listing with known-editor citations and the "possible editing or noisy assembly"
   mark.
4. **A per-site source column** for flags in listed lineages: paper, or database record.
5. **Completeness reporting**, with "not scored" never reported as zero.
6. **Assembly units in the strong rule** (conservative). A first version merged by bare organism name and was corrected
   before the analysis ran.
7. **Real-data sensitivity per listed lineage**, with *nad3*-174 separate and an out-of-sample rate.
8. **Each flag's own chance of flagging**, by Monte Carlo from the planted null.
9. **The pre-approved, conditional second-assembly step for VdLs.17.** It was approved after a partial look in which
   VdLs.17 had E = 0.21, and ran only because it came out strong.
10. **The removed-break listing** for glass sponges and *Polyrhachis*.
11. **The paired-break tally** (a safety net for the pair filter).
12. **Refinements logged at 09:21.**
    - The second-assembly length threshold.
    - Species-level searches.
    - Documented sites in non-core genes left out of the sensitivity counts.
13. **The *Grigorovia* steps, all post hoc:** steps a–g, the blinded re-check, the held-out congener, the tRNA checks,
    and the intron and congener checks prompted by the referee.
14. **Text changes after the fact-check and the referee review** (`protocol/corrections.md`).
15. **The six follow-up tests (§8c).**
    - Predictions were logged before each test.
    - Methods changes made while running are logged with reasons:
      - a one-reference-per-gene completeness run (the full run would have taken about 4 hours per genome);
      - relatives' and congener queries, DIAMOND, read checks and a local assembly for the nuclear genes;
      - intron positions from annotation and spliced alignments, because the scan's insert ends are unreliable where
        an intron ORF continues the exon;
      - CLIB 448 instead of the type strain, which has no mtDNA, plus read baiting for its mitogenome;
      - a profile-likelihood bound on ω;
      - site calls in the congeners against their spliced CDS;
      - read versions counted in fixed windows from alignments.
    - All numbers were then regenerated from the final scripts.
    - After the referee review:
      - codon-model and outgroup sensitivity of the ω test;
      - the GGC count with ancestral residues;
      - read counts by strand;
      - reciprocal best hits.
    - **Forking paths.**
      - Anomalies that favored H2 (apparent losses and frameshifts in nuclear genes) were pursued until resolved.
      - Results that favored H1 were stress-tested (codon models, orthology) only after the referee asked.
      - Test 1's outcome at *cox1* favors H2 by its logged criterion, and is reported that way in §8c.

**Some of these definitions were written after partial looks at the data.** These included peeks at codes 1–4, the
vertebrate batch with its grebe and turtle flags, and a partial real-site rate. None changes the frozen single-genome
flags.

## 10. Errors made and corrected during the work

All are in the `protocol/frozen_rules.md` log.
- **The B2 diagnosis.** I first said the missed annotated shifts were genes split between two hits. They were shifts
  misplaced in poorly conserved stretches.
- **Timestamps.** Several log times were typed ahead of the clock. They were corrected to commit times in a later entry,
  without editing the original entries.
- **An empty file read as zero.** A mid-scan peek read a still-empty results file as "0 shifts" for two glass sponges.
  That line was withdrawn, and the completeness rule was added.
- **The SRA search.** It first counted *Grigorovia* runs by one genus taxid and missed the three congeners, which NCBI
  files under another taxid. It was corrected.
- **The blinded check's references.** It was first logged as using "different references". In fact both protein sets
  are in the scan's references.
- **Process slips.** A tRNA-finder parsing error, a tool waiting at an input prompt, and a wrong DOI in a draft. All
  were caught before release.
- **In the follow-up tests:**
  - The first intron mapping used the scan's insert ends, which can be hundreds of nt off, and missed two introns.
    It was replaced (§8c).
  - An early test-1 summary put the A of AGG among the best-scoring placements at *cob* 26266. It scores lower, though
    it gives the same protein.
  - An earlier draft said the GCC tRNA-Gly was unique among "101 relatives". The 101 included *G. jiainica* itself;
    the comparison is 100 Saccharomycetaceae in the scan plus three congeners.
  - Pileup counts at indels suggested heterozygous frameshifts that read-level grouping showed to be alignment
    artifacts (§8c).
- **A 611 MB intermediate file was committed** to the local repository by mistake. The local history after commit
  d6e06a9 was rewritten to remove it, with messages, authors and dates preserved. No hash cited before the rewrite
  changed; the old-to-new table is in the `protocol/frozen_rules.md` log.
- **Review findings in the follow-up section:**
  - a stale intron distance in one table;
  - a test 1 row claiming "the same protein" where two tied placements differ;
  - "intact" used where only presence was tested;
  - claims stronger than the evidence (listed in `protocol/corrections.md`).

## 11. Limitations

- **The sensitivity figures.**
  - The real-data sensitivity rests on 8 lineages and 7 out-of-sample sites. The planted rates overstate it.
  - Placement errors in poorly conserved stretches probably cost more than the filters do.
- **The design.**
  - The listed-lineage exemption is by family name. It let *G. jiainica* through as "known", and would do the same for
    any other off-codon flag in a listed family. Birds (1,107 genomes), turtles, Saccharomycetaceae (101) and glass
    sponges were never tested for new codons.
  - The pool test adds nothing here.
- **For *G. jiainica*.**
  - No RNA, protein or respiration measurement. A loss of function more recent than a few percent of its branch
    (about 1–5% depending on the codon model and the outgroups, measured in synonymous change), including one that
    arose in culture, cannot be excluded.
  - One strain and one sequencing run. Its nuclear genome is heterozygous.
  - Collection stocks can lose function in culture. One congener's type stock has lost its mtDNA.
  - The tRNA-Gly(GCC) is one lineage, noticed after the case was chosen.
  - The ω bound assumes a uniform synonymous clock along the branch and depends on the codon model. *cob*'s
    synonymous divergence on that branch (dS 1.82 under F3x4, 0.71 under F1x4) is model-dependent and near saturation
    under F3x4.
  - Nuclear gene calls depend on detecting divergent proteins (24–32% identity for some activators). Orthology rests
    on reciprocal best hits, without synteny. PET122 finds no *S. cerevisiae* match, and CBS1 was found nowhere,
    controls included.
  - 100% of reads shows the culture's mtDNA is uniform; it does not date the bases.
  - Several literature sources were read only as abstracts or short excerpts (`protocol/search_log/` and `candidate/search_log/`).

## 12. References

Short form. Full entries, with how deeply each was read, are in `protocol/prior_art.md`, `protocol/search_log/` and `candidate/search_log/` and
`candidate/followup/test4_phenotypes.md`.

- Al Arab M et al. (2017) *Mol Phylogenet Evol* 106:209–216. doi:10.1016/j.ympev.2016.09.024
- Andreu-Sánchez S et al. (2021) *GigaScience* 10:giaa161. doi:10.1093/gigascience/giaa161
- Beargie C et al. (2008) *Bioinformatics* 24:2571–2578. doi:10.1093/bioinformatics/btn487
- Belcourt MF, Farabaugh PJ (1990) *Cell* 62:339–352. doi:10.1016/0092-8674(90)90371-K
- Chen et al. (2012): editing-site prediction in *Didymium* (see `protocol/prior_art.md`)
- Farabaugh PJ, Kramer E, Vallabhaneni H, Raman A (2006) *J Mol Evol* 63:545–561. doi:10.1007/s00239-005-0311-0
- Fekete V et al. (2007) *FEMS Yeast Res* 7:1237–1247
- Gornik SG et al. (2022) *Mol Biol Evol* 39:msac191. doi:10.1093/molbev/msac191
- Gott JM, Parimi N, Bundschuh R (2005) *Nucleic Acids Res* 33:5063. doi:10.1093/nar/gki820
- Haen KM, Pett W, Lavrov DV (2014) *Gene* 535:336–344. doi:10.1016/j.gene.2013.10.041
- Jourda C et al. (2015) *Genome Announc* 3:e00823-15. doi:10.1128/genomeA.00823-15
- Langkjaer RB et al. (2003) *Nucleic Acids Res*. doi:10.1093/nar/gkg423
- Lee CF et al. (2008) *FEMS Yeast Res* 8:114–118 (description of *Kazachstania jiainicus*)
- Lin et al. (2026) bioRxiv (*Heteranomia* TAG frameshifts; see `protocol/prior_art.md`)
- Malimas et al. (2023) Jxiv preprint 334 (*Pachytichospora*)
- Masuda I, Matsuzaki M, Kita K (2010) *Nucleic Acids Res* 38:6186–6194. doi:10.1093/nar/gkq449
- Mikata K, Ueda-Nishimura K, Hisatomi T (2001) *Int J Syst Evol Microbiol* 51:2189. doi:10.1099/00207713-51-6-2189
- Mindell DP et al. (1998) *Mol Biol Evol* 15:1568–1571
- Opulente DA et al. (2018) *BMC Biol* (yeast trait data)
- Opulente DA et al. (2024) *Science* 384:eadj4503. doi:10.1126/science.adj4503 (Y1000+ phenotype data)
- Pantou MP, Kouvelis VN, Typas MA (2006) *Curr Genet* 50:125–136. doi:10.1007/s00294-006-0079-9
- Rosengarten RD et al. (2008) *BMC Genomics* 9:33. doi:10.1186/1471-2164-9-33
- Russell RD, Beckenbach AT (2008) *J Mol Evol* 67:682–695. doi:10.1007/s00239-008-9179-0
- Sarhan MS et al. (2025) *Sci Data*. doi:10.1038/s41597-025-06447-x
- Szabóová D, Hapala I, Sulo P (2018) *FEMS Yeast Res* 18:foy028. doi:10.1093/femsyr/foy028
- Wolters JF et al. (2023) *Front Microbiol* 14:1268944. doi:10.3389/fmicb.2023.1268944
- Yang J et al. (2017) *Genome Biol Evol* 9:1161–1174
