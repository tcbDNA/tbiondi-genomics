# Pan-metazoan check of tRNA-Tyr(GUA) position 33 and mtRF1a against mitochondrial UAA=Tyr

**DRAFT, 2026-11-24. The independent cold review is pending and will be run separately.** The original brief and PREREG §7 ask for that review before this report is written; the user asked for the draft first, so the order is a recorded deviation. Figures and full outlier dossiers are not made yet.

Evidence labels used throughout: **[C]** confirmatory (PREREG plus committed Gate 3, non-panel genomes only), **[E]** exploratory, **[S]** sensitivity, **[P]** post-hoc. "Panel" genomes are the 360 calibration genomes; they never enter confirmatory counts.

## Abstract
We tested two claims about why some animal mitochondria read UAA as tyrosine: that tRNA-Tyr(GUA) carries a U33C change (Lin et al. 2026; Zhao et al. 2026), and that nuclear mtRF1a carries T191S (Zhao et al. 2026; human T216). We screened 32,982 animal mitogenomes with rules fixed in `PREREG.md` before the run. **[C]** One non-panel genome met the Tyr rule, *Rotylenchulus reniformis*, and it is blocked from counting as a new clade by the locus audit of its tRNA-Tyr. **[C]** 17 of 28,583 readable UAA=Stop species (0.059%, 95% CI 0.037–0.095%) have a non-U33 tRNA-Tyr. After post-hoc audits, 11 of the 17 sit in lineages with support, 4 are unconfirmed and 2 were dropped **[P]**. C33 is therefore not sufficient for UAA=Tyr. **[C]** None of 77 readable starfish is C33. **[S]** *Rotylenchulus* mtRF1a has the T216S change Zhao described in other lineages, confirmed in independent reads, but its Tyr call is sensitivity-labelled. S216 recurs in three unrelated Tyr lineages, and neither U33C nor T216S is necessary.

## 1. Introduction
Zhao et al. 2026 (BMC Biology, `10.1186/s12915-026-02767-5`) proposed that tRNA-Tyr(GUA) U33C and mtRF1a T191S drive UAA reassignment in several lineages. Lin et al. 2026 (bioRxiv `10.64898/2026.08.14.744788`) reported the tRNA change as exclusive to reassigned lineages among the species they examined (they call it "U36C"; canonical position 33). Li et al. 2019 reported a C in starfish tRNA-Tyr. This study checks the tRNA claim across all public mitogenomes and, for flagged species, the release factor.

Prior art that bears on *Rotylenchulus*: Jacob et al. 2009 (read) established UAA=Tyr in *Radopholus similis* and reported it in *R. arabocoffeae*. Vovlas et al. 2026 (Plants 15:1423, `10.3390/plants15101423`; full text read) state that translating COI of *R. borealis* showed the use of the flatworm mitochondrial code (NCBI table 14, in which UAA is tyrosine). They infer that two internal TAGs in a short haplotype "probably" code for tyrosine, "as does TAA stop codon in Radopholus similis", and they say the same occurs in other *Rotylenchulus* sequences. They also consider heteroplasmy or a pseudogene. They cite Palomares-Rius et al. 2021 and Torres et al. 2022, which I could not read. **COI-barcode work therefore recognised the unusual code in *Rotylenchulus* before this study.** What this study adds is mitogenome-wide, read-confirmed internal UAA and the release-factor genotype. A GenBank `transl_table` 14 on a record is **not evidence** of reassignment (§8).

## 2. Methods, briefly
Full rules are in `PREREG.md` (commit `5aa5774`); code and outputs are in the commits listed in §10.
- **Mitogenomes:** 32,982 included (RefSeq, DToL, INSDC-only and other assemblies); 98 excluded with reasons. Confirmatory set 32,623. 8,189 UNVERIFIED-flagged and 5,635 other incomplete INSDC-only species were not eligible.
- **UAA call:** frame-aware counts of internal TAA at profile-aligned codons and conserved-Tyr columns, a log-likelihood, and Codetta (`-m`). Tyr needs ≥5 internal TAA in ≥3 genes, ≥3 at conserved-Tyr columns, LLR ≥3 and Codetta Y.
- **tRNA-Tyr position 33:** Infernal with a tRNA-Tyr model, corroborated by ARAGORN or tRNAscan-SE, then a locus audit (≥4 of 5 stem pairs, ≤10 nt overlap, agreement with the annotation). Calls that fail only the stem rule are "stem_audit_review".
- **mtRF1a:** human `NP_061914.3` positions 203 and 216; annotated proteins, then TSA, then assemblies (miniprot), then reads; paralogs and bacterial RF1 excluded by blastp and tree.

### Detection floor and what the locked rules can miss
- In-sample recall is 20 of 20 positive genomes, but that is about 5 origins and 15 distinct sequences. Leave-one-lineage-out recall falls to 3 of 5 and 2 of 5 for *Pseudotanais* when positives set the thresholds.
- The planted-null recall is 0, 24, 90, 95 and 97% when 5, 10, 25, 50 and 100% of conserved Tyr codons are replaced. A partly reassigned lineage is likely to be missed.
- 0 of 340 panel negatives were called Tyr. The 95% upper bound on that rate is about 1.1%, so up to ~370 false Tyr calls are possible among ~33,000 genomes. The sweep made one Tyr call.
- The tRNA-Tyr call is a borderline pass: 9.5% of panel negatives with an annotated trnY failed extraction, and 12.1% with wrong-locus calls added. Cnidaria, Ctenophora and Chaetognatha have 0 readable tRNA-Tyr.

## 3. Pre-registration adherence, deviations and disclosures
`PREREG.md` has not been edited since the lock. Everything else is listed here.
1. **Interim looks.** Four interim looks at the sweep at the user's request (looks 1–4; look 3 reversed an earlier instruction that day), and a continuous 5-minute watcher for the last part of the sweep. The committed classifier was used; no rule or threshold was changed because of them. The watcher was also run on the UNVERIFIED set. Looks are logged in `results/gate3/interim/LOG.md`.
2. **Post-lock bug fix, adopted.** Origin-spanning annotation intervals were stored as covering the whole genome, so good trnY calls looked wrong-locus. Found after look 2, adopted as primary on the user's instruction, with the without-fix run as sensitivity. The candidate list and the Tyr and weak calls are the same with and without it.
3. **Follow-up code.** The weak audit, A1/B1, T6, post-hoc stem and candidate code were written after the lock; some after look 2.
4. **Order of the report.** The draft precedes the cold review (see the top).
5. **Not run:** A2, B2, the §5b exploratory items (incomplete T/TA stops, backup TAG, extra tRNAs), the strand-skew analysis and the read search for mtRF1a in *Microdevario* and Serrivomeridae. B1 was not testable. B3 and B4 are descriptive.
6. **Corrections logged during the work**
   - The first candidate read check tested the **wrong alternative base** for minus-strand genes; it was corrected before the numbers reported here.
   - The conspecific search first used ungated calls (16S fragments and junk scores), later redone with the gate.
   - *Albula*: the Gate 3 report says no candidate has its own nuclear assembly. That is wrong for *A. glossodonta* (three assemblies). `AP002973` is not a replicate of `CM165726` (87% identity) and its species label is doubtful; it is counted once, as a named extra.
   - Two *Serrivomer*-family records are mislabelled (`OR546184` is a frogfish; `OR482581` is *S. sector*).
   - `GCA_056783515.1`, an accession supplied by the user, is *Draco maculatus*; the bonefish is `GCA_056783665.1`.
   - The Gate 4 size estimate in Gate 3 followed the wider original brief, not PREREG §2.4.
   - `transl_table` 14 was first read as evidence of reassignment; it is not (§8).
   - The five *Rotylenchulus* transcript sets are one project (PRJNA286314), so they count once.
   - The UNVERIFIED weak-audit script first omitted the Gate 3 Tyr genomes; the fix changed only *Pratylenchus goodeyi*.
   - The `PV962` fish batch is a suspected bad submission (§7).
   - Smaller items: an early 10.3% trnY failure figure omitted wrong-locus calls (12.1% with them); `gate3_posthoc_stem.py` totals included the named extra; a hard-coded wrong taxid fixed before the lock; one read-result file was overwritten and restored from the log; the COX1 fragment downloads were added to the ledger late; Gate 2 used data from outside the frozen panel in several logged cases (GATE_2 A4).

## 4. Confirmatory results [C]
### 4.1 UAA
| Call | Confirmatory (n = 32,623) |
|---|---|
| Stop | 32,510 |
| ambiguous / insufficient | 48 / 61 |
| weak | 3 (all rejected by the pre-registered audit) |
| **Tyr** | **1** |

The one Tyr genome is *Rotylenchulus reniformis* `CM003310.1`. **No new reassigned clade is claimed under the locked rules**: the trnY of `CM003310` is stem_audit_review (stem 3/5, score 12.8), so the locus audit blocks it.

**Caveats on that null.** (i) The eligibility filter excluded 8,189 UNVERIFIED species, a flag that in-frame stops can cause (the known reassigned *Anomia* and *Placuna* records carry it), so the confirmatory set is biased against reassigned lineages; the exploratory run (§7) only partly answers this. (ii) The detection floor and leave-one-lineage-out recall in §2. (iii) Most genomes are arthropods or chordates; *Nematoda* 394, Bryozoa 44, Acanthocephala 41 (incomplete tRNA coverage).

### 4.2 Contingency table (UAA by tRNA-Tyr position 33, confirmatory)
| UAA | U33 | non-U33 | stem_audit_review | unresolved |
|---|---|---|---|---|
| Stop | 28,566 | **17** (14 C33, 3 G33) | 2,018 | 1,909 |
| Tyr | 0 | 0 | 1 | 0 |
| weak | 3 | 0 | 0 | 0 |

Non-U33 among resolved Stop: 17 of 28,583, 0.059% (Wilson 95% CI 0.037–0.095%). By origins (NCBI taxonomy, minimum changes): 11 independent origins of non-U33 among Stop, 1 origin of Tyr, so **no 2×2 test was run** (pre-registered minimum of 3 new Tyr origins). T2 precision 0 of 17, recall undefined. By phylum: Chordata 13 of 10,795; Nematoda 2 of 352; Porifera 2 of 143; all others 0. With the 5-nt and 0-nt overlap cut: 17 of 26,784 and 15 of 21,376.

### 4.3 Coverage of tRNA-Tyr (resolved / genomes; bird-heavy stem reviews excluded from "resolved")
Arthropoda 14,498 / 15,047; Chordata 10,799 / 12,956 (1,959 stem_audit_review, mostly birds); Mollusca 1,492 / 1,625; Platyhelminthes 427 / 486; Annelida 463 / 477; Nematoda 353 / 394; Echinodermata 323 / 326; Porifera 144 / 186; Bryozoa 25 / 44; Acanthocephala 20 / 41; Rotifera 8 / 21; Cnidaria 0 / 892; Ctenophora 0 / 17; Chaetognatha 0 / 6. Others are in `results/gate3/tests_fixed.txt`.

### 4.4 Named tests
T4: 77 of 77 resolved asteroids are U33 ("Zhao-consistent"), including *Asterias rubens* `LR745847.1`; 0 non-U33 in 323 resolved echinoderms. T6: lineage-level (Albuliformes). A1: 1 of 1 origin has U33-majority sisters. B1 not testable.

## 5. Claims and their limits
### a. C33 is not sufficient for UAA=Tyr
**[C]** 14 non-panel Stop species with C33 (the 17 minus the 3 with G33) meet the locked definition of a sufficiency counterexample; `results/gate3/tests_fixed.txt` lists all 17. **[P]** Audits after the sweep, by lineage:

| Lineage | Support | Limits |
|---|---|---|
| *Microdevario* (2 species) | all 9 confident calls in the genus C33; reads 572 C33 to 1 U33 and 10 to 0; independent records of *M. kubotai* | no nuclear assembly, so NUMT not checked |
| Albuliformes (4 genomes; 1 panel, 1 named extra) | all C33; independent *A. vulpes* mitogenome; reads 469 to 0; 0 nuclear copies over the trnY site in 2 *A. glossodonta* assemblies (hap2 not searched; *A. vulpes* has none); a planted 10%-diverged copy was found | small; named-extra species label doubtful; *Pterothrissus* is one 2004 record, not audited further |
| Serrivomeridae (4 species) | independent specimens of *S. sector*, *S. beanii*, *Stemonidium*; reads 363 to 0, 272 to 0, 85 to 0; *S. lanceolatoides* is its own specimen only (419 to 0) | no nuclear assembly; NUMT not checked |
| *Lycenchelys verrillii* | 4 records, reads 622 C33 to 1 U33, 0 nuclear copies | an 82 kb mitochondrial-like contig was excluded |

The eels hold **3 or 4 separate origins** of C33 (4 by taxonomy, 3 on a weakly supported gene tree); Albuliformes is a further lineage whose monophyly the tree does not support. **Unconfirmed**, one specimen with nothing independent: *Coloconger*, *Ophichthus asakusae*, *Xyrias*, *Toxocara* (and the post-hoc crane). **Dropped** after audit: *Globodera ellingtonae* (chromosome II has another GUA tRNA-like locus with U33, score 11.1 against 10.0) and *Bufotes pewzowi* (no mitochondrial reads in its runs, no independent record).

### b. Craniella is a non-U33 Stop lineage (G33)
*C. zetlandica* and *C. infrequens* are UAA=Stop with G33 (4/5 stem, two callers) **[C]**. **[P]** Support: identical trnY in both; DNA reads 55 to 0 and 40 to 0; RNA 35 to 0, 66 to 0 (a separate study) and 2 to 0; U33 in *Cinachyrella* (3), *C. arb* (216 reads to 0) and *C. villosa*; the two *Craniella* mitogenomes are 97.1% identical, so this is most simply one origin. mtRF1a is canonical (Q203, T216) **[E]**. The two records come from the same Sanger pipeline (different specimens).

**Nuclear-copy caveat:** the assemblies hold 217 and 322 nuclear mitochondrial segments; all 4 and 38 that cover the site also read G. No U33 read was seen, but the source copy of the record cannot be proven from these data.

### c. Neither U33C nor T216S is necessary for UAA=Tyr (calibration lineages)
Panel genomes, so descriptive.
- **U33C not necessary:** *Pseudotanais* sp. (OZ281197, OZ281199) are resolved U33 and Tyr; two more are unresolved. *Radopholus* is UAA=Tyr; its U33 tRNA-Tyr is not confidently called by our rules (score 16.5), and Zhao states that *R. similis* lacks the change.
- **T216S not necessary:** *Radopholus* (K203, T216, a 2-residue insertion) and Plicatulidae (*Plicatula muricata*, T216; *P. plicata* T216 from 3 clean reads of `SRR3350473`).
- *Pseudotanais* mtRF1a was not tested by us (Zhao lists it among S191 lineages). Its U33 result does not speak to T216S.
- Our *Plicatula muricata* trnY is C33 in raw data but not confidently called (stem 3/5).

### d. S216 recurs in three unrelated Tyr lineages
Anomiidae + Placunidae, *Cephalodiscus* (both Zhao's, re-derived here in Gate 2) and ***Rotylenchulus*** **[S]**.
- All 48 *R. reniformis* gene models (2 genome loci; 40 from TSA sets of one project) give Q203 and S216, loop `EGGVHRVQRLPFTEKTGRVHSSTVSVA`. Reads from two other projects: 25 (RNA-seq) and 6 (genome) with codon TCC, 0 with any Thr codon. The gene groups with nematode RF1.
- Stop relatives read T216 and Q203: *Pratylenchus vulnus* (confirmatory Stop), *Meloidogyne hapla* (exploratory Stop) and *Rotylenchus robustus* (Hoplolaimidae; **no mitogenome call**, Stop inferred from 6 COX1 fragments without internal UAA and a congener). *Hoplolaimus columbus* was a bacterial contaminant and is untested.
- **Limits.** *Rotylenchulus* is sensitivity-labelled: its trnY fails the stem audit (the trnY is U33, with a weak score) and it is not a counted new clade. Whether S216 is sufficient is untested: **no Stop species with S216 was found**, and the Tier B C33/G33 Stop species with nuclear data all have T216, so the combination C33 or G33 plus S216 is untested. *Microdevario* and the Serrivomeridae are untested (no nuclear assembly).

### e. At least two origins of UAA=Tyr in tylenchid nematodes
*Radopholus* (Pratylenchidae; K203 with a 2-residue insertion, T216) and *Rotylenchulus* (Hoplolaimidae; Q203, S216, no insertion) have Stop relatives in each family (by NCBI taxonomy) and different mtRF1a changes, so at least two origins are implied **[S]**. **[E]** In 10,294 tylenchid COX1 fragments, 91 of 92 *Rotylenchulus* fragments in 10 or more species have internal UAA. No other genus shows a fixed pattern: *Hoplolaimus* (about 220) and *Rotylenchus* (about 60) have none; *Helicotylenchus* and *Scutellonema* are mixed and unconfirmed; *Radopholus* shows UAA in 1 of 12 fragments, so the fragment test has low power and a negative family is weak evidence. *Pratylenchus goodeyi* is a supported-weak call only (§7). Rotylenchulus internal TAGs: COX3 column 239 is at a Tyr site (8 of 8 readable relatives), COX1 404 and CYTB 62 are not; only COX1 404 appears in RNA-seq. UAG reading is unresolved.

### f. The starfish C33 claim is not supported
**[C]** Li et al. 2019 report C33 in starfish tRNA-Tyr. None of 77 resolved asteroids, and none of 323 resolved echinoderms, has a non-U33 call; *Asterias rubens* is U33. Limits: one asteroid is unresolved, and the species in Li et al. were not identified or fetched, so the claim is not contradicted for their specimens in particular.

### g. No new reassigned clade under the locked rules
See §4.1. The sufficiency counterexamples mean that the pre-registered **null result** (no new clade *and* no counterexample) does not hold; the clade half does. Read the null with the UNVERIFIED-filter caveat, the detection floor and the leave-one-lineage-out results in §2.

## 6. mtRF1a by species [C for PREREG §2.4 scope, S/E otherwise]
| Species | UAA | 203 | 216 | Status |
|---|---|---|---|---|
| *Anomia* ×2, *Enigmonia*, *Placuna vitream* | Tyr | Q | S | Gate 2, panel |
| *Cephalodiscus* ×2 | Tyr | Q | S | Gate 2, panel and reads |
| *Plicatula muricata*, *P. plicata* | Tyr | Q | T | Gate 2, Gate 4 |
| *Radopholus* | Tyr | K | T | Gate 2 |
| *Rotylenchulus* | Tyr (blocked) | Q | S | **[S]** |
| *Rhabdopleura annulata*, *Heteranomia*, *Pododesmus*, 36 background species | Stop | K / I / Q / Q | T | Gate 2 |
| *P. vulnus*, *M. hapla* (*R. robustus*: no call) | Stop | Q | T | **[S]** |
| *Albula*, *Lycenchelys*, *Craniella* ×2 and 5 sisters (*Megalops*, *Anguilla*, *Lycodopsis*, *Zoarces*, *Cinachyrella*) | Stop | Q | T | **[E]** Tier B |

B3 (S216 by origin): 3 of 5 Tyr origins with loop data carry S216 with *Rotylenchulus*, 2 of 4 without. 0 Stop species with data carry it. Descriptive only.

## 7. Exploratory, sensitivity and post-hoc
- **[E] UNVERIFIED run:** 7,707 of 8,189 species (480 for ambiguous bases, 2 excluded by PREREG). 0 Tyr; 12 weak (11 rejected by the locked audit, 1 supported weak); 1 non-U33 Stop, *Desmacella* (A33, judged a repeat artefact); 1 non-U33 of 6,506 resolved Stop. Phylum coverage is in GATE_3B. Reads contradict the UAA of 3 of the rejected weak records.
- **[E] Suspected bad submission: PV962 fish batch:** 216 records, 20 odd calls (9.3%), against 0 of 304 other bony fish.
- **[S] *Pratylenchus goodeyi*:** supported weak; 35 internal UAA, 9 of 28 conserved-Tyr columns, LLR −2.68. No other specimen or read exists in NCBI. The ARAGORN UUA-anticodon locus in `ON525593` is not evidence (29 of 444 other nematode genomes have one).
- **[P] Bird and turtle stem re-inclusion:** the Aves (1,575 calls) and Testudines (161) mismatch in the trnY stem is shared; with them T1 is 18 of 30,307 (the extra is a crane, *Antigone rubicunda*, which has one record and a duplicate, a 3/5 stem and no reads). The criteria were fixed after look 2, so this analysis is not blind.
- **[P]** strand-skew analysis: not run.

## 8. Limitations
1. Detection floor and LOLO (§2). UAA thresholds were set with positives in view.
2. The sweep cannot see partly reassigned lineages or mitogenomes excluded by completeness (§4.1).
3. tRNA-Tyr is unreadable in Cnidaria, Ctenophora and Chaetognatha, and under-sampled in Nematoda, Acanthocephala and Bryozoa.
4. Non-U33 candidates are single species; taxonomy stands in for phylogeny; the eel tree is weak.
5. The non-U33 and Craniella audits are post-hoc and were run after seeing the calls.
6. Nuclear copies (NUMTs) could not be checked in *Microdevario*, the Serrivomeridae, the unconfirmed singletons or *Albula vulpes*.
7. `transl_table` 14 on a GenBank record is filing practice: *Ditylenchus* and *Radopholus* fragments without UAA carry it too. The "41 of 49 *Rotylenchulus* fragments" earlier treated as support is therefore not support by itself; the internal UAA codons are.
8. The *Rotylenchulus* fragments may include pseudogenes or heteroplasmy (Vovlas et al. raise both). The genome and read evidence are independent.
9. Panel genomes support claims c and part of d, so those claims are not new findings for those lineages.
10. The release factor in *Pseudotanais*, *Placuna placenta*, *Microdevario* and Serrivomeridae is untested.

## 9. Discussion
The two-part model, in which U33C and T216S both matter, is not contradicted or supported by this sweep. C33 is not sufficient: four lineages with support (Microdevario, Albuliformes, Serrivomeridae, Lycenchelys) have C33 with UAA=Stop, and Craniella has G33 with UAA=Stop. Neither change is necessary in the calibration set. S216 appears in three unrelated Tyr lineages, and this fits a convergent change linked to the code, but Plicatulidae and *Radopholus* reach UAA=Tyr without it, and no Stop species with S216 was found. The *Rotylenchulus* result is the closest this study has to a new data point, and it rests on a genome that the pre-registered rules do not count.

## 10. Files and commits
PREREG lock `5aa5774` (SHA256 `e23bc13f…`); Gate 3 `7e0541f`; exploratory `418b0ba`; Gate 4 `243c24f`. Gate reports: `GATE_0_VERIFY.md` to `GATE_4_MTRF1A.md`. Data: `results/mitogenome_master.tsv`, `results/gate3/`, `results/gate3_unverified/`, `results/gate4/`. Downloads are in `data/manifests/download_ledger.tsv` (88.37 GB in 613 rows). The report is a draft and uncommitted.

## 11. Summary of cold review
Pending. No critique has been received.
