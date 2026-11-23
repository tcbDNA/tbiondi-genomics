# A scan for codon-specific frameshifts in mitochondria

Ribosomes read a gene three bases at a time. A few mitochondria break this rule at one particular codon. This scan
found no new case under its frozen rules; a follow-up outside those rules turned up one exploratory candidate, the
yeast *Grigorovia jiainica*. The known cases:
- glass sponges skip a base at UGG;
- the oyster parasite *Perkinsus* skips at AGG (and at CCC);
- some budding yeasts of the *Kazachstania* group (here called the CTG-shifting yeasts) skip at CTG;
- in a 2026 preprint, an anomiid bivalve skips at TAG.

Their genes look broken in the DNA, but ribosomes read past the shifts (for the yeasts, this is proposed from DNA
alone). We asked whether other mitochondria do the same.

We aligned 20,343 mitochondrial genomes (17,648 from RefSeq, all but land plants, plus 2,695 fungal genomes) to
reference proteins with a frameshift-aware aligner (DIAMOND). The test applied to each genome (defined under Scan
results) was fixed in advance.

**Under that test, the one flag outside the lineages listed in advance as known was an assembly error, so the test
found no new lineage.** The test is not sensitive. In the genomes sampled, which are mostly vertebrates and
arthropods, and in the core genes (the 14 standard mitochondrial protein genes):
- documented real frameshift sites are detected as clean shifts at their own codon about 29% of the time (mean of 8
  lineages; bootstrap 9–54%; pooled 32%);
- so a genome with 2, 3 or 5 such sites would be flagged about 8%, 20% or 44% of the time (ranges 1–29%, 2–56%,
  7–86%). Planted shifts are flagged more often (see Sensitivity).

The list of known lineages names families and orders with documented sites (among them glass sponges, *Perkinsus*, the
Saccharomycetaceae for the CTG-shifting yeasts, turtles and anomiid bivalves), plus all birds (`protocol/frozen_rules.md`).
Listed lineages were not tested for new codons.

One yeast candidate, found outside the frozen rules, survived every check that its extra bases are real; it is
described below as exploratory.

![Part of the cob gene of Grigorovia jiainica: two extra bases (an A and a T), each just before a GGC codon. Removing them restores the reading frame and the residues of the S. cerevisiae protein; read straight on, the gene stops nine codons later. DNA only; the cause of the extra bases is unknown.](plots/grigorovia_cob_site.png)

*Two of the four extra bases, in* cob *of* Grigorovia jiainica *(GenBank BK073075.1). DNA only: the cause of the
extra bases is not known.*

---

**How this was made.** The scan, the analysis and the write-up were done by a combination of Claude Fable 5.2 and
Claude Opus 5.5 agents (Anthropic), steered by Thomas Biondi. The approach follows
[ry-xu/genome-research](https://github.com/ry-xu/genome-research) (Ryan Xu): a frameshift-aware scan, rules fixed
before the data, controls with pass thresholds ("bars"), then read checks and a blinded re-check. The checking was
also done by AI; no human expert has reviewed the work. It included a blinded re-check of the *G. jiainica* sites with
another aligner (miniprot; its reference proteins overlapped with the scan's), three literature searches (prior art,
*G. jiainica* under all its names, and the growth records of its congeners, the other *Grigorovia* species), and a
fact-check and a hostile referee review, each run on the whole draft and again on the *G. jiainica* section
([`protocol/corrections.md`](protocol/corrections.md)). The work is computational only, with no lab work: SRA holds
only DNA for *G. jiainica* and its three congeners, no RNA, and proteomics repositories were not searched. The
*G. jiainica* follow-up used about 1.6 GB of public reads and six nuclear assemblies.

- **Rules fixed in advance.** The single-genome test and the scan's rules were fixed in a local git repository
  (`protocol/frozen_rules.md`). They were not externally timestamped. Many descriptive additions came after the freeze,
  some after partial looks at the data; all are logged with times.
- Corrections and criticism are welcome as GitHub issues.

### Scan results

| | genomes | result |
|---|---|---|
| Scanned (RefSeq release 238 without land plants; Sarhan et al. 2025 fungi) | 20,343 | all scored |
| Flagged in listed lineages | 15 | glass sponge *Aphrocallistes vastus* (TGG); both *Perkinsus* (AGG); 7 CTG-shifting yeasts (CTG; 5 of them records without gene annotation); 2 grebes and 2 turtles (AGT); and *G. jiainica* (AGG), counted as known only because its family is listed |
| Flagged outside listed lineages, strong | 1 | *Verticillium dahliae* VdLs.17: an assembly error. An independent *V. dahliae* assembly carries the missing base at all 38 of its shift positions; a *V. nonalfalfae* assembly carries it at 37, and the 38th cannot be compared |
| Expected by chance under the single-genome test | 2.5 (upper bound) | |

The single-genome test needs at least two clean shifts at the same codon, in different genes, more often than randomly
placed shifts would give. A shift is clean when it is not near a gene end, an intron edge, a run of one base or an
opposite shift, and is flanked by enough aligned sequence. A flag is strong when it is expected by chance less than
once across the scan, or the same codon flags in a related genome. Of the listed lineages that could be checked:
- 1 of the 4 glass sponges in the scan flagged, so prediction P2 (another glass sponge flags) failed;
- 6 of the 18 scanned genomes in the four genera of CTG-shifting yeasts flagged at CTG (4 of 16 besides the two
  named positives, *Huiozyma sinensis* and *Monosporozyma servazzii*, so prediction P3 passed), and so did
  *Jamesozyma lodderae*.

Details: [`protocol/scan/flags_summary.md`](protocol/scan/flags_summary.md),
[`protocol/scan/scan_report.md`](protocol/scan/scan_report.md), [`protocol/scan/predictions.md`](protocol/scan/predictions.md).

### Sensitivity of a negative result

- **Planted shifts.**
  - A genome with 2, 3 or 5 planted shifts at one codon is flagged 32%, 48% and 72% of the time.
  - When every planted site is counted at its codon, a genome with two sense-codon sites is flagged 76 of 76 times
    (stop codons: 0.976, n = 42). So the shortfall comes from detecting single sites.
- **Real sites.** Documented sites in known lineages are counted less often: 0.29 against 0.55 for planted shifts.
  - That figure is an unweighted mean of 8 lineages; the pooled rate is 0.32.
  - Only 7 sites played no part in setting the rules. 3 of those were detected, and for them the codon could not be
    checked.
  - The bird and turtle *nad3*-174 position (598 annotated sites; gained and lost many times) is counted 61% of the
    time.
- **Causes of missed sites.**
  - filters for runs of one base (homopolymers), where sequencing errors pile up;
  - shifts placed at the wrong codon in poorly conserved stretches;
  - genes that do not align;
  - a pair filter that removes real sites where +1 and +2 shifts sit close together (14 of 38 *Perkinsus* sites).
- **Scope.** Most scanned genomes are animals; protists and algae contribute 686. Only core genes were tested, and
  recall of shifts near fungal introns is 0.28.

Source: [`protocol/scan/lineage_sensitivity.md`](protocol/scan/lineage_sensitivity.md),
[`protocol/controls/BARS.md`](protocol/controls/BARS.md).

### *Grigorovia jiainica* (exploratory; found outside the frozen rules)

The frozen rule counted this flag as known, because its family (Saccharomycetaceae) is listed for the CTG-shifting
yeasts; its codon, AGG, is not theirs. We followed it up after seeing the result, so everything here is exploratory.
The species is recorded as not growing on non-fermentable carbon such as glycerol (Lee et al. 2008), so a loss of
function in *cox1* and *cob* was a live alternative to a bypass (the extra bases being skipped in translation or
removed from the RNA).

The predictions for each follow-up test were logged before it ran
([`candidate/grigorovia_jiainica.md`](candidate/grigorovia_jiainica.md)). Some inputs had already been seen: the codons of the
congener *G. yakushimaensis* at two sites, and the distances from the sites to *S. cerevisiae* intron positions.

**Four extra bases.**
- They sit in *cox1* (position 24711) and *cob* (26496, 26477 and 26266) of BK073075.1, strain NRRL Y-48843. Read without
  skipping, three form an AGG codon (the extra A plus the first two bases of the GGC after it) and one, *cob* 26477, a
  TGG.
- Each is in 100% of the reads that span it, from the run the assembly was built from:
  - 1,440 to 1,785 reads per site (1,151 to 1,381 with PCR duplicates removed), on both strands in similar numbers,
    none without;
  - 1,516 reads (1,220 with duplicates removed) span both *cob* 26496 and 26477 (19 nt apart) and carry both bases
    ([`candidate/funnel/step_f_reads_dedup.md`](candidate/funnel/step_f_reads_dedup.md)).
- Removing each one restores the reading frame.
- All three congeners lack them, at all four sites:
  - *G. yakushimaensis*, deposited as OR703946;
  - *G. humatica* and *G. transvaalensis* (strain CLIB 448), whose mitogenomes we assembled from public reads.
- One site, *cob* 26477, sits 2–3 nt from an intron in *G. yakushimaensis* and *G. humatica*. Compared across their
  spliced genes, they lack it too.
- None of the three has a clean frameshift or a stop in any core gene.

**Nuclear genes, selection and the rest of the mitogenome.** The genes were under constraint for all but possibly the
last few percent of the branch, and none of the other genes we checked looks broken.
- **Nuclear genes.** We found no broken copy of the nuclear genes for making and using Cox1 and Cob. The nuclear
  genome is heterozygous (diploid or higher); we screened its heterozygous sites and resolved each flagged spot from
  the reads. The genes checked:
  - the translational activators that exist only to make these two proteins (PET309 and MSS51 for Cox1; CBS2 and
    CBP1 for Cob). For each, the best hit among all *S. cerevisiae* proteins is the expected gene;
  - the nuclear subunits of complexes III and IV;
  - heme a synthesis and cytochrome c;
  - a sample of mitochondrial translation genes.
- **Genes that could not be called** (none is missing in *G. jiainica* alone). The controls are *M. servazzii* and *K.
  africana*, each recorded growing on a non-fermentable carbon source in at least one record (*M. servazzii*'s records
  disagree).
  - **CBS1** (a COB activator): missing in all four *Grigorovia* species and in both controls.
  - **CBP2:** missing in all four *Grigorovia* and in *M. servazzii*; present in *K. africana*. In *S. cerevisiae*
    CBP2 is needed only to splice the *cob* bI5 intron, which no *Grigorovia* *cob* has (*G. jiainica*'s has no
    introns at all).
  - **PET54:** missing in all four *Grigorovia* (apart from one weak partial match in *G. jiainica*) and in
    *M. servazzii*; present in *K. africana*.
  - **AEP2:** not found in any of the six genomes by miniprot, the controls included. There are partial DIAMOND
    matches in *G. jiainica* and *G. humatica*; DIAMOND was not run on the controls.
- **Selection.** *cox1* and *cob* are under strong purifying selection on its branch: dN/dS (ω) 0.031, like its
  congeners (0.013–0.044). The tree includes three outgroups: *M. servazzii*, *Arxiozyma bovina* and *A. telluris*.
  - Neutral evolution (ω = 1) is rejected under every codon model tried.
  - Any period without constraint is at most about 1–5% of the branch, depending on the codon model and the outgroups
    (*cox1* alone: 5–29%).
  - That share is measured in synonymous change, not in years.
- **The rest of the mitogenome.** We found no damage:
  - *atp6*, *atp8*, *atp9*, *cox2* and *cox3* are full length, without frameshifts or stops;
  - *VAR1* is an open frame;
  - both rRNAs and the RNase P RNA are present;
  - the tRNA set matches *G. yakushimaensis*'s, except tRNA-Gly.

**So the cause is unresolved: the extra bases may be bypassed, or the genes may have broken very recently, perhaps
even in culture.**
- A bypass could be a +1 shift during translation, like the one proposed from DNA alone for its CTG relatives (sites
  in *cob*, *cox2*, *cox3* and others), or RNA editing.
- The tests above cannot tell a bypass from a loss within roughly the last 1–5% of the branch: both predict intact
  machinery and a low ω. A loss in the last fraction of a percent, or in culture, would leave too few changes to
  detect at all.
- The pattern of the insertions (below) is the main argument for a bypass, but part of it follows from how the sites
  were found. The glycine tRNA raises a decoding problem that cuts the other way (below).

**The pattern of the insertions (noticed after the fact).**
- **GGC.** *G. jiainica*'s seven core mitochondrial genes (*cox1*, *cob*, *cox2*, *cox3*, *atp6*, *atp8*, *atp9*; the
  scan found no *nad* genes in it) contain 4 GGC codons in 1,815 (5.4% G or C at third codon positions), and all four
  sit at the extra bases. These are the only GGC codons in *cox1* and *cob*: one in *cox1* and three in *cob*, each at
  an extra-base site; the other five genes have none. An extra base can be placed right before the GGC at every site,
  but against *G. yakushimaensis* that placement scores best only at *cob* 26496 (tied). At *cob* 26266 the best
  placements are three bases earlier (score 75; the placement before GGC scores 70; both give Thr-Gly). At *cox1*,
  removing a G instead scores 85 against 80 and gives the congener's serine (see Caveats). *cob* 26477 cannot be
  compared ([`candidate/followup/test1_sister_codons.md`](candidate/followup/test1_sister_codons.md)).
- **How the sites were found.** The scan flagged *G. jiainica* because shifts recurred at AGG, so at three sites the
  extra A and the GG after it are part of the signal that selected them. The part the scan did not select for is the C
  that completes GGC. C is the third base of 3.7% of codons in the seven genes (67 of 1,815), and of none of the other
  105 glycine codons ([`candidate/followup/genus/codon_positions.md`](candidate/followup/genus/codon_positions.md)).
  Two sites are 19 nt apart, so the four are not independent, and no null model for the pattern has been run.
- **Congeners and outgroups.** They have GGT, GGA or GGG at these four codons, so each event pairs an extra A or T
  with a change of the next codon to GGC, a codon the genes use nowhere else. The same pairing at all four sites
  points to a process tied to these places.
- **tRNA-Gly.** Its only mitochondrial glycine tRNA has anticodon GCC (1,329 reads in an exact-match count; the 4
  reads with TCC are sequencing errors, see Glycine decoding). That is unique among the four *Grigorovia* species and
  the 100 other Saccharomycetaceae in the scan, which all have UCC. A GCC anticodon pairs exactly with GGC, the codon
  each skip produces (at *cob* 26266, under the placement in AGG).

**Glycine decoding (noticed after the fact).**
- By standard pairing a GCC anticodon reads GGC and GGT, not GGA, and the seven genes have 38 GGA codons
  ([`codon_positions.md`](candidate/followup/genus/codon_positions.md)).
- There is no second, TCC tRNA-Gly in the mitogenome. ARAGORN, even at a relaxed threshold, and tRNAscan-SE each find
  one tRNA-Gly, with GCC. Of 1,249 reads at its anticodon, 4 read TCC (0.32%), a share matched or exceeded by the
  C-to-A error share at 74 of 564 other C positions nearby. Three of the four are contradicted by another read of the
  same DNA molecule (the mate, or PCR copies), which reads GCC; the fourth has 8 mismatches. No copy of the
  mitochondrial tRNA-Gly was found in the nuclear assembly
  ([`candidate/funnel/trna_gly_tcc.md`](candidate/funnel/trna_gly_tcc.md)).
- So either something else decodes GGA in the mitochondrion (for example a tRNA imported from the cytoplasm), which
  would weaken the tie between the GCC tRNA and the sites, or nothing does, which would impair translation of every
  mitochondrial gene and fit a loss of function. The growth records cannot decide this: two congeners with intact
  genes and a UCC tRNA-Gly also fail to grow on non-fermentable carbon (Caveats).

**Caveats.**
- **Serine at the *cox1* site.** By the rule logged in advance, the *cox1* site counts against a bypass, because the
  congener first compared (*G. yakushimaensis*) has serine there. Skipping a G instead of the A would reproduce it.
  But that serine is unique to *G. yakushimaensis*: the other congeners and all three outgroups have glycine.
- **Age of the bases.** 100% of reads does not date them, because a culture's mtDNA becomes uniform quickly.
- **Nuclear evidence.** Intact nuclear genes show that the machinery exists, not that Cox1 and Cob are made. A search
  with PET122 (needed for Cox3) finds no *S. cerevisiae* match, so its orthology rests on relatives' proteins.
- DNA only: no RNA, protein or respiration measurement.
- It was found by overriding the family rule after seeing the flag.
- **Growth records may describe damaged stocks.**
  - The species, and two congeners with intact genes (*G. yakushimaensis*, *G. humatica*), are recorded as not
    growing on glycerol, ethanol or lactate.
  - One type-strain stock (*G. transvaalensis* NRRL Y-17245) has lost its mtDNA: 4 of 200,000 reads match
    mitochondrial proteins, against 991–2,933 in the other three strains.

**What would settle it:**
- RNA across the four sites;
- COX1 or COB peptides spanning them;
- an oxygen-uptake or petite test on NRRL Y-48843;
- a PCR across the sites in an independent stock of the type strain SF1S05 (CBS 10587, ATCC MYA-4454). Whether NRRL
  Y-48843 is the same isolate is not recorded in what we found, and this would also test for a change that arose in
  culture.

Files: [`candidate/grigorovia_jiainica.md`](candidate/grigorovia_jiainica.md) (the six follow-up tests and the checks
added after review), [`candidate/funnel/`](candidate/funnel/) (sites, reads, blinded re-check, tRNAs) and
[`candidate/heldout/`](candidate/heldout/) (the deposited *G. yakushimaensis* mitogenome, run as a held-out relative). Coordinates
and context for lab tests:
[`candidate/funnel/sites_for_lab.md`](candidate/funnel/sites_for_lab.md).

### Prior work

**Codon-specific mitochondrial frameshifts.**
- glass sponges at UGG ([Rosengarten et al. 2008](https://doi.org/10.1186/1471-2164-9-33)). Across glass sponges such
  insertions were "gained and lost repeatedly", which the authors read as tolerated under relaxed selection
  ([Haen et al. 2014](https://doi.org/10.1016/j.gene.2013.10.041), abstract).
- *Perkinsus* at AGG and CCC ([Masuda et al. 2010](https://doi.org/10.1093/nar/gkq449);
  [Gornik et al. 2022](https://doi.org/10.1093/molbev/msac191)).
- *Kazachstania*-group yeasts, the CTG-shifting yeasts, at CTG ([Langkjaer et al.
  2003](https://doi.org/10.1093/nar/gkg423); [Szabóová et al. 2018](https://doi.org/10.1093/femsyr/foy028); [Wolters
  et al. 2023](https://doi.org/10.3389/fmicb.2023.1268944)).
- *Heteranomia*, an anomiid bivalve, at in-frame TAG (Lin et al. 2026, bioRxiv).
- single programmed sites such as *nad3* in birds and turtles (Mindell et al. 1998; Russell & Beckenbach 2008).

**Earlier published searches** (`protocol/prior_art.md`).
- frameshift-aware annotation of animal mitogenomes (Al Arab et al. 2017; at most one shift per gene);
- editing-aware gene finders for myxomycete mitogenomes (Gott et al. 2005; Beargie et al. 2008; Chen et al. 2012);
- annotation tools that flag frameshifts (MFannot);
- a survey of *nad3* across birds and turtles (Andreu-Sánchez et al. 2021).

None we found tested for shifts piling up at one codon across many genomes.

**This species.**
- Described as *Kazachstania jiainicus* (Lee et al. 2008). Placed in *Grigorovia* in 2020, a name NCBI marks "nom.
  inval.", and combined as *Pachytichospora jiainica* in a 2023 preprint.
- We found no report of frameshifts in it under any name.
- Its mitogenome was assembled and annotated with MFannot by Sarhan et al. 2025. Only the sequence was deposited; the
  annotation is not public, and we did not run MFannot on it.

**Nearest cases at AGG elsewhere.**
- In a fungal mitochondrion we found a database annotation, with no paper describing it, in the chytrid *Rhizophydium*
  *cob*. It is a different event: the C after AGG is skipped and AGG keeps its arginine.
- In the yeast nuclear genome, the Ty1 retrotransposon shifts +1 at CUU-AGG-C (Farabaugh et al. 2006). None of the four
  *G. jiainica* sites has the CUU.

### Open points

- **The cause of the *G. jiainica* extra bases.** No experiment so far separates a translational skip, RNA editing and
  a very recent loss of function.
- **Glycine decoding.** How *G. jiainica* reads GGA with only a GCC tRNA-Gly (see Glycine decoding above; the GCC
  tRNA-Gly is seen in one lineage and was noticed after the case was chosen). Whether a nuclear tRNA-Gly is imported
  cannot be checked from DNA data; it would need RNA from purified mitochondria.
- **How recent a loss could be.** An indel-to-substitution comparison in neutral yeast mtDNA could put a number on
  "very recently". It was not done.
- **Limits of the scan:**
  - real-site sensitivity rests on 8 lineages and 7 out-of-sample sites;
  - the known-lineage list works by family, not by codon, so it counted *G. jiainica* as "known";
  - the pool test (genomes pooled by family and order) only repeated single-genome signals.

### Files

- `notebook.md`: the running notebook. Method, controls and missed bars, the frozen rules, results, failed predictions,
  every change after the freeze, errors made and corrected, what could be wrong, references.
- `protocol/`: the frozen rules and their timestamped log (`frozen_rules.md`), the public corrections log
  (`corrections.md`), the design, prior art and accessions, and the scan's results, controls and audit. See
  [`protocol/README.md`](protocol/README.md).
- `candidate/`: the *G. jiainica* follow-up. It holds the summary (`grigorovia_jiainica.md`), the funnel steps and
  blinded re-check (`funnel/`), the six follow-up tests (`followup/`), the held-out congener (`heldout/`) and its
  literature search. See [`candidate/README.md`](candidate/README.md).
- `plots/`: the site figure.
- `code/`, `scripts/`: runnable with `uv run`. `scripts/setup_tools.sh` installs the external tools,
  `scripts/setup_data.sh` fetches the public data by accession, and `scripts/run_all.sh` reruns everything.

The format of this write-up (frozen rules, running notebook, public corrections log) follows ry-xu/genome-research
(https://github.com/ry-xu/genome-research).
