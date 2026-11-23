# Review log: fact-check and referee review of the release draft (2026-11-22)

Three fresh agents reviewed the draft. None of them edited any file.
- **Two fact-check runs.** Every number and sentence in README.md and notebook.md was checked against `protocol/` and `candidate/` and the
  cited sources (DOIs via Crossref; abstracts and some full texts).
- **One hostile referee review** of the main claims and their wording.

Each change is listed below. A snapshot of the draft before the changes is kept in the working repository's history.
The same items are logged in `protocol/frozen_rules.md`.

## Factual corrections (fact-check)

| # | Where | Was | Now | Source |
|---|---|---|---|---|
| 1 | README, notebook, frozen-rules log | The blinded re-check used "other" / "different" reference proteins | Both of its reference sets (*S. cerevisiae*, *S. ludwigii*) are in the scan's references; only the aligner was independent | `candidate/funnel/blind_key.json`; scan events against NC_037703 |
| 2 | README | "The only earlier frameshift scan ... was of animals (Al Arab 2017)" | Earlier searches listed: Al Arab 2017; editing-aware gene finders (Gott 2005, Beargie 2008, Chen 2012); MFannot flags; the *nad3* survey (Andreu-Sánchez 2021). None tested for a codon pile-up | `protocol/prior_art.md` §4 |
| 3 | README, notebook | Known lineages omitted *Heteranomia* | Added (Lin et al. 2026, bioRxiv); *Magnusiomyces* bypassing mentioned in notebook | `protocol/prior_art.md` §1 |
| 4 | notebook | "only it has clean AGG shifts" | "in core genes": *M. aerobia* has one clean AGG shift in a non-core ORF | `candidate/funnel/fam_trna_gly.json` |
| 5 | README | "*nad3*-174, one inherited site" | "one homologous position, gained and lost many times" | `protocol/prior_art.md` §2 |
| 6 | notebook | "For *Vazella*, the sites are in Haen 2014" | Status unknown: the source paper's abstract names no species or positions | `protocol/scan/removed_breaks_known_groups.md` |
| 7 | README, notebook | "Sister species" used for two different species | "Congener" | no phylogeny in `protocol/` and `candidate/` |
| 8 | README, notebook | Name status not given; *Pachytichospora* called a former name | Described as *K. jiainicus* (2008); *Grigorovia* 2020 (nom. inval.); *Pachytichospora* in a 2023 preprint | NCBI Taxonomy; search log |
| 9 | README, notebook | "Assembled by Sarhan et al. without gene annotation" | Annotated with MFannot by the authors; the annotation is not public, and we did not run MFannot | `protocol/design.md` amendment 1; Sarhan methods |
| 10 | README | "Neither has [the congener]" | The congener lacks the extra base at the three comparable positions; the fourth is at its intron; the written prediction failed | `candidate/heldout/summary.md` |
| 11 | notebook | "They are recent and unique" | "probably recent" | `candidate/heldout/summary.md` |
| 12 | README, notebook | Blinded task "told only to find frameshifts"; "the known CTG shifts"; "glycine at each" | It was also asked for the residue and a systematic/scattered judgment. The CTG yeast has no documented sites (the agent found 5, the scan 3). Glycine holds under its preferred placement only; at 26266 miniprot faces Thr | `candidate/funnel/blind/TASK.md`, `candidate/funnel/blind/VERDICT.md` |
| 13 | README | "In the DNA, not the assembly" | Read support from the same run excludes assembly and read errors in this sample; it is necessary, not sufficient | `protocol/design.md` decision 13 |
| 14 | README | "flagged every time (76 of 76)" when every site is "detected" | When every site is *counted* at its codon: 76 of 76 with two sense sites; stop codons 0.976 (n = 42) | `protocol/controls/BARS.md` |
| 15 | README | "would very likely catch a *Perkinsus*-like lineage (7 clean shifts)" | Removed. The text gives the computed chances with ranges, the scope, and the fact that listed lineages were exempt | `protocol/scan/lineage_sensitivity.md` |
| 16 | notebook, scan report | Paired-break rows "AT-rich codons in divergent arthropod genomes"; identity "0.31–0.77" | Mostly AT-rich, in arthropods, nematodes, cnidarians and others; identities up to 0.93 (correction note added to the report) | `protocol/scan/paired_break_tally.md` |
| 17 | README | *Rhizophydium* the "nearest AGG case" (unhedged) | Hedged; noted as a different event (the C after AGG skipped, AGG stays Arg) | NC_003053.1 |
| 18 | README | Haen 2014 cited for UGG | Rosengarten 2008 for UGG; Haen 2014 cited for "gained and lost repeatedly" under relaxed selection (abstract) | abstracts |
| 19 | README | "No RNA or protein ... in public databases" | SRA holds only DNA; proteomics repositories were not searched | `candidate/funnel/completeness_and_data.md` |
| 20 | README, notebook | Process claims not documented in `protocol/` and `candidate/` | This file and the frozen-rules log now document them | — |
| 21 | notebook | "each with a stop rule written before it ran" | Steps a, b and f only; the others were descriptive | `protocol/frozen_rules.md` 10:03 |
| 22 | notebook | Frozen strong rule included "from a different assembly" | Marked as a post-freeze, conservative addition | `protocol/frozen_rules.md` 09:14, 09:24 |
| 23 | notebook | VdLs.17 "29" and "38" side by side; "both carry the base at 38" | 29 core shifts on the flagged scaffold; 38 clean shifts on three scaffolds; NC_008248 carries 38 of 38, NC_029238 37 of 37 comparable | `protocol/scan/second_assembly_VdLs17.tsv` |
| 24 | notebook, scan report | CTG yeasts "4 of 10" against "4 of 11" | Explained: 11 sites, 10 checkable at their codon | `protocol/scan/lineage_sensitivity.md` |
| 25 | notebook | *K. viticola* duplicate cited to the wrong file | `protocol/scan/scan_report.md` | — |
| 26 | notebook | Incomplete references (Farabaugh, Lee, Fekete; Pantou, Malimas, Belcourt missing) | Completed from Crossref metadata | Crossref |
| 27 | README, notebook | "cob 172" | *cob* reads on 6, 70 and 172 codons between and after its sites | `candidate/funnel/step_a_sites.md` |
| 28 | notebook | Section 9 omitted the 09:21 refinements | Added | `protocol/frozen_rules.md` 09:21 |
| 29 | README | "about 29%" without qualification | "unweighted mean of 8 lineages; pooled 32%" | `protocol/scan/lineage_sensitivity.md` |
| 30 | README, notebook | "686 protists" with no table in `protocol/` and `candidate/` | `protocol/scan/scan_composition.tsv` added | — |
| 31 | data | `step_c` (2,012 codons, AGA 38) against `agg_usage` (1,748, AGA 33) | Note added: all kept alignments against core genes only | — |
| 32 | data | SF1S05 type status "not checked" | Type strain (SF1S05T = CBS 10587T; Lee et al. 2008) | Lee 2008 |
| 33 | data | Search log said the OR703946 annotation reads cox1/cob straight through | Clarified: the authors' figshare annotation files, not the GenBank record, which has no CDS | OR703946 |
| 34 | notebook | B2: 155 misplaced against about 151 implied by the rates | Explained: 4 sites detected at the site have no matching CDS | recount |

## Changes from the referee review

- **README opening, rewritten.**
  - It separates the confirmatory result (the only strong flag outside listed lineages was an assembly error) from the
    exploratory *G. jiainica* case.
  - It gives the negative its scope (mostly vertebrates and arthropods, core genes, listed lineages not tested for new
    codons) and its intervals.
- ***G. jiainica*, reframed as "cause unknown".** The three readings (translational skip, RNA editing, disabling
  mutations) have equal standing, and the opening now states:
  - the growth data;
  - the congeners' phenotypes;
  - that the GCC tRNA-Gly cannot read GGA by standard pairing.
- **The evidence list, split** into "the bases are real" and "what bears on the cause".
  - The blinded re-check moved under "real", with its limits.
  - The tRNA-Gly(GCC) coincidence moved to open points, described as n = 1 and noticed after the case was chosen.
- **New post hoc checks** (`candidate/funnel/introns_and_congener.md`).
  - Three of four sites lie 3–12 nt from *S. cerevisiae* intron insertion sites.
  - The congener has Ser at the *cox1* site, which an alternative placement of the extra base reproduces.
- **Placement stated precisely.** It is fixed at the two *cob* sites and ambiguous at the other two.
- **Ty1.** Its site needs CUU before AGG, which none of the four sites has.
- **Static figure.** Title made neutral ("Two extra bases in cob (DNA only)"); the footer states that the cause is
  unknown; "Without the skips" became "Read straight on".
- **Process.** The README says some post-freeze definitions came after partial looks at the data. Notes §9 says the
  VdLs.17 step was approved after a partial look, at E = 0.21.

## Recommended by the referee and not done in the first round (for the user to decide)

The nuclear-gene check, the dN/dS test and the intron map were done afterwards as follow-up tests (second round, below).

- Check the nuclear genes for respiration and for *cox1*/*cob* translation (PET309, MSS51, CBS1, CBS2, CBP1, and
  others) in the nuclear assembly GCA_030571715.1.
- dN/dS on the *G. jiainica* branch. This needs mitogenomes assembled from the congeners' public reads, and the
  project's rules exclude assembly from reads.
- An intron map across the four *Grigorovia* species.
- A neutral figure with congener rows; dropping the animation. The animation was kept because the user asked for it,
  labeled as a model. (Second round: the animation was removed from the release at the user's request.)
- Possibly cutting the *G. jiainica* section to a short "candidate, cause unknown" paragraph.

# Second round (2026-11-22, after the six follow-up tests on *G. jiainica*)

Two fresh agents reviewed the revised text. Neither edited any file.
- **A fact-check** of every changed passage (README opening, figure caption, the *G. jiainica* section, open points;
  notebook §8b–8c, §9–11, references; `protocol/README.md`; the new `protocol/frozen_rules.md` log entries) against `protocol/` and `candidate/`.
- **A hostile referee review** of the new *G. jiainica* section, whose claim had become stronger.

## Fact-check findings and changes

| # | Where | Finding | Change |
|---|---|---|---|
| 1 | README; this file | The README said both reviews were run again, but this log covered only the first round | This section added |
| 2 | notebook §8c, test 1; `candidate/followup/test1_sister_codons.tsv` | "Same protein as every top placement" at *cob* 26266 was wrong (two tied placements give Met). The TSV's "fits" column required Gly | Text: "Thr-Gly, as under every top placement". The column now compares residues (`scan_placement_residue_matches_congener`) |
| 3 | `candidate/followup/genus/compare.md`; `candidate/followup/test2_summary.md`; notebook | The intron distance at *cob* 26477 appeared as 10, 2 and 3 nt | `candidate/followup/genus/compare.md` now refines junctions near sites (2 nt, tied 2–3). Text everywhere: "2 or 3 nt (tied; 3 under the group I rule)" |
| 4 | notebook §10 | "All are in the frozen-rules log": two errors were not | Logged in `protocol/frozen_rules.md` |
| 5 | notebook §8b | "Only *G. yakushimaensis* had a mitogenome assembly" | "Of the three congeners, only…" |
| 6 | `protocol/frozen_rules.md` | A 12:21 entry refers to a "12:23 entry" (should be the 12:11 entries); a later entry says "25–35%" (should be 24–32%) | Corrected by new log entries (the log is not edited in place) |
| 7 | README opening | "Survived every computational check we could run" overstates: test 1 favors H2 at *cox1* by its logged rule, and the held-out congener prediction failed | Now: "survived every check that its extra bases are real". The original was the user's wording and is flagged to the user |
| 8 | notebook §8c test 1 verdict | "Mixed" against the logged criterion | "H2 at *cox1* by the logged criterion", with the derived Ser noted |
| 9 | notebook §8c test 2 verdict | "Against the scar version" was stronger than `candidate/followup/test2_summary.md` | "No precise scar; a weak link to intron positions" |
| 10 | notebook §8b | "More than 100 nt from any intron" (true only for *S. cerevisiae*) | "From any *S. cerevisiae* intron" |
| 11 | README nuclear bullet | "Probably diploid"; translation was only sampled; CBP2, PET54 and AEP2 not mentioned | "Heterozygous (diploid or higher)"; "a sample of mitochondrial translation genes"; CBP2, PET54 and AEP2 listed as not callable |
| 12 | README | "The two respiring controls" lacked support in `candidate/followup/test4_phenotypes.md` | Control rows added to `candidate/followup/test4_phenotypes.md`; the README now says "the two controls" |
| 13 | README clue | The GGC placement is not the best-scoring at *cox1* either | Rewritten (GGC pattern below; 26266 noted) |
| 14 | notebook §8b | Type-strain stock "does not grow", while Y1000+ records glycerol growth for CBS 2186 | Both records given; "stocks evidently differ" |
| 15 | README, notebook | "Intact" for RNA genes and *VAR1*, where only presence and an open frame were tested | "Shows no damage", with each kind of evidence named |
| 16 | README | "No frameshift or stop": the table counts clean shifts and stop columns | "No clean frameshift or stop" |
| 17 | README, notebook | "Predictions logged before each test" without the seen-before caveat | Caveat added |
| 18 | notebook §10, `protocol/frozen_rules.md` | "No hash cited changed": rewritten hashes are now cited | "No hash cited before the rewrite changed" |
| 19 | notebook references | Mikata 2001 and Opulente 2024 were not cited in the text; Opulente 2024 was incomplete | Cited in §8c; full reference |
| 20 | notebook §8b | "Probably compiled from that description" | "Curated from *The Yeasts* (2011)" |
| 21 | README | "Two literature searches" | "Three" |
| 22 | notebook §9 | Item 15 omitted two methods changes | Added |
| 23 | `code/README.md`, `site_figure.py` | Stale references to the removed animation | Removed, including unused GIF code in `site_figure.py` |
| 24 | notebook §8c | "Ordinary mitogenome" for CLIB 448 | "Has mtDNA (37 kb assembled, incomplete)" |
| 25 | `candidate/README.md`; `candidate/followup/test3_nuclear_genes.md` | Unlisted files; 372 + 2,056 = 2,428 against 2,319 | Files listed; "2,428 groups in the reference file; compleasm scores 2,319" |
| 26 | README; `candidate/followup/test3_nuclear_genes.md` (found while applying the fixes) | "CBP2, PET54 and AEP2 were found in no *Grigorovia*": DIAMOND finds partial AEP2 matches in two species and a weak PET54 match in one | "Not found by miniprot; only partial DIAMOND matches"; listed gene by gene in `candidate/followup/test3_nuclear_genes.md` |

## Referee findings and changes

| # | Finding | Change |
|---|---|---|
| R1 | Tests 3, 5 and 6 predict the same under a bypass and under a recent loss; "the genes appear to be in use" overstates | Heading now: "The genes were under constraint until recently, and nothing else looks broken". The text says the tests cannot tell a bypass from a very recent loss (including in culture), and that the lean toward a bypass rests on the insertion pattern |
| R2 | The ω bound depends on the codon model (F1x4, F61, no outgroups), and the *cob* result is model-dependent | Rerun in `candidate/followup/genus/dnds_sensitivity.md`. Text: "at most about 1–5% of the branch, depending on the codon model and the outgroups", measured in synonymous change; the *cob* result is called model-dependent |
| R3 | The tRNA-Gly(GCC) clue did not draw its own implication; the test 6 prediction was near-certain | GGA-decoding caveat added (either another decoder, or impaired translation of every gene, which fits a loss); the near-certainty is noted in §8c |
| R4 | The *cox1* Ser is derived: other congeners and outgroups have Gly | Recomputed (`candidate/followup/genus/ggc_sites.md`) and stated |
| R5 | The strongest pattern was unreported: the 4 GGC codons are the 4 sites | Recomputed (`candidate/followup/genus/ggc_sites.md`) and added as the stated basis for the lean |
| R6 | Orthology of low-identity activators not established | Reciprocal best hits added (`candidate/followup/nuclear/rbh.md`): all expected, except PET122 (no hit), stated. Synteny not done |
| R7 | CBS 10587 and ATCC MYA-4454 are the type strain, not "another strain" | Corrected, and recorded in `candidate/search_log/grigorovia_literature.md`. A PCR of an independent type-strain stock is suggested |
| R8 | 100% of reads does not date the bases | Caveat added |
| R9 | Read counts at the two close *cob* sites; strand balance | Recomputed (`candidate/followup/reads_strand.md`): balanced strands; 1,516 reads span both sites with both bases, none without |
| R10 | Forking paths: H2-favoring anomalies were pursued, H1-favoring results were not stress-tested; the test 6 H1 prediction was partly written by the agent | Stated in notebook §9 and §8c |
| R11 | An indel-to-substitution analysis would put a number on "very recently" | Not done; listed as open |
| R12 | RNA editing that deletes a base may lack precedent in fungal mitochondria | Not adopted: we did not check the literature for it, so skip and editing are still listed without ranking |

## Final edits to the *G. jiainica* section (user request, after the second round)

- **The user's five edits.**
  - Congeners lack the extra bases at all four sites; *cob* 26477 sits 2–3 nt from an intron in two congeners.
  - CBS1, CBP2, PET54 and AEP2, gene by gene, with the controls.
  - The bypass sentence.
  - "A GCC anticodon pairs exactly with GGC, the codon each skip produces (at *cob* 26266, under the placement in
    AGG)".
  - The serine caveat, now stated against the logged rule.
- **Three phrases adjusted to match the data before checking:**
  - the 26477 site sits 2–3 nt from an intron, not inside one;
  - the 26266 placement qualifier;
  - "the congener first compared", not "the nearest congener", because *G. jiainica* is the earliest-diverging
    *Grigorovia*.
- **A fact-check of these lines only** found edits 1, 4 and 5 supported, edit 2 supported with two gaps, and edit 3
  overstated.
- **Applied, at the user's request:**
  - "at all four sites", since all 12 site-by-congener comparisons were made (`candidate/followup/genus/compare.md`);
  - CBP2 "needed only to splice the *cob* bI5 intron, which no *Grigorovia* *cob* has", because two congeners do have
    a *cob* intron (bI2-type) and still lack CBP2;
  - the controls sentence notes that *M. servazzii*'s growth records disagree (`candidate/followup/test4_phenotypes.md`);
  - the bypass sentence now reads "a +1 shift during translation, like the one proposed from DNA alone for its CTG
    relatives (sites in *cob*, *cox2*, *cox3* and others), or RNA editing". The CTG-yeast evidence is DNA only
    (`protocol/prior_art.md`).

## Register pass with the declauding skill (user request)

**Scope.** README.md and notebook.md were edited with the `declauding` skill (oaustegard/claude-skills, version 0.17.0).
`protocol/frozen_rules.md` (a verbatim log) and this file were not edited.

**Steps run:** lint; sentence and structure pass; a structural review of headers, openings and closers by a separate
agent; two newcomer reads by separate agents; and a claim diff against the pre-edit files.

**Register changes:**
- Coy or verdict headers became labels ("Scan results", "Sensitivity of a negative result", "Prior work",
  "Limitations", "Sites, reads and the blinded re-check").
- Staged closers and "X, not Y" reversals were rewritten as statements.
- A fragment run was merged.
- Absolute negatives were scoped to what was checked ("none of the other genes we checked").

**Clarity changes from the newcomer reads:**
- Terms defined at first use: clean shift, strong flag, core genes, the known-lineage list, bypass.
- The aligner named (DIAMOND; miniprot for the blinded check).
- "CTG yeasts" renamed "CTG-shifting yeasts", because "CTG clade" already names *Candida* yeasts.
- The *G. transvaalensis* strain named (CLIB 448).
- The AGG/GGC reading explained.
- P3's outcome and named positives stated.
- Outgroups named.

**Factual fixes found during the pass:**

| # | Where | Was | Now | Source |
|---|---|---|---|---|
| 1 | README, scan results | "Independent assemblies carry the missing base at all 38" | One *V. dahliae* assembly at 38; *V. nonalfalfae* at 37, the 38th not comparable | `protocol/scan/second_assembly_VdLs17.tsv` |
| 2 | README | "6 of the 18" against P3's "4 of 16", unexplained | The 2 named positives stated; P3 passed | `protocol/scan/predictions.md` |
| 3 | README | "a loss … would leave no trace", while the notebook says "the last few percent" | The tests cannot exclude a loss in the last 1–5%; one in the last fraction of a percent leaves too few changes to detect | `candidate/followup/genus/dnds_sensitivity.md` |
| 4 | README headline | "the only strong flag outside" implied weak flags existed | "The one flag outside … was an assembly error, so the test found no new lineage" | `protocol/scan/flags_summary.md` |
| 5 | notebook §8b | Congener intron "at the TGG site" | "2–3 nt from the TGG site" | `candidate/followup/genus/refine_cob_junction.md` |
| 6 | notebook §6 | The loss table did not add up; only some lineages were named | Smaller categories and the other lineages pointed to | `protocol/controls/BARS.md`, `protocol/scan/lineage_sensitivity.md` |

**Left as written, with reasons:**
- **Sentences the user dictated** in the *G. jiainica* section.
- **The process record in notebook.md:** "I", "the user", "the agent", the first control run, dropped approaches, and
  "added after review". The notebook is a record of the work by design, and post hoc analyses must be marked.
- **Length of the *G. jiainica* section and its gene-by-gene list:** user-specified.
- **The author-name placeholder** and the missing issue link: for the user at publication. (The name was filled in
  when the repository was prepared: Thomas Biondi.)

2026-11-23 07:08 EST, user request: the README's GGC bullet now states that the four GGC codons are the only GGC codons in *cox1* and
*cob* (one in *cox1*, three in *cob*, each at an extra-base site; none in the other five core genes). Checked against
`candidate/followup/genus/ggc_sites.md` and a recount of the *G. jiainica* rows of the codon alignments
(`codon_alignments/`): *cox1* 529 codons, 1 GGC; *cob* 384 codons, 3 GGC; the other five genes 902 codons, 0 GGC.
