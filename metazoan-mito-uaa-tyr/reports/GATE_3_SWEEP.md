# Gate 3: Pan-metazoan mitogenome sweep

Date: 2026-11-24. **Stopped before Gate 4.** Sweep run 2026-11-23 16:51 to 2026-11-24 08:39; all 32,982 genomes have one valid result. Raw output: `results/gate3/run/out.jsonl`; master table `results/mitogenome_master.tsv`; tests `results/gate3/tests_fixed.txt` (primary) and `tests_nofix.txt`; follow-ups in `results/gate3/post.txt`, `posthoc_stem.txt`, `candidate_*.txt`, `roty_*.txt`, `notes_albula_records.md`. PREREG.md is unchanged from `5aa5774`.

## Read this first: what weakens or qualifies these results
1. **Post-lock bug fix.** Origin-spanning annotation intervals made good trnY calls look wrong-locus (GATE_2 A15). Adopted as primary on the user's instruction after interim look 2; the without-fix run is the sensitivity check. The candidate list and the Tyr/weak calls are **identical** with and without it.
2. **Four interim looks and a continuous watcher were run during the sweep at the user's request**, after the code was committed but before it finished; logged in `results/gate3/interim/LOG.md`. The follow-up code (weak-tier audit, A1/B1, T6, post-hoc stem analysis) was written after the lock; the stem criteria and candidate-check code after look 2. No PREREG rule changed. The manifest order was changed at 07:53 (non-arthropods first); order cannot affect results.
3. **The INSDC completeness filter is a bias toward missing reassigned animals.** It excluded 8,189 species flagged `UNVERIFIED` and 5,635 others (mostly partial). GenBank uses UNVERIFIED for records it could not verify, which in-frame stops can cause: *Anomia* and *Placuna* are themselves UNVERIFIED (Zhao says UAA codons were removed in them), and among relatives of the new candidate, *Rotylenchus*, *Hoplolaimus*, *Nacobbus* and *Pratylenchus goodeyi* are too. **A null for reassignment is therefore weaker than the 32,623 tested genomes suggest.** See "Scope of the null".
4. **Tested genomes are mostly vertebrates and arthropods** (28,000 of 32,623). Minor invertebrate phyla have 1–1,600 genomes.

## Results (confirmatory = non-panel, n = 32,623 representatives; 32,982 included in total, 98 excluded with reasons)
| UAA class | n |
|---|---|
| Stop | 32,510 |
| ambiguous | 48 |
| insufficient coverage | 61 |
| weak | 3 |
| Tyr | **1** |

trnY position 33: **28,647 resolved** (with the fix), 2,020 stem_audit_review, 1,956 unresolved (1,034 gate fail, 363+342 wrong locus, 217 no call). Locus audit passes at 0 / 5 / 10 nt: 21,421 / 26,846 / 28,647.

### T1. Non-U33 among confident UAA=Stop species
**17 of 28,583 (0.059%, 95% CI 0.037–0.095%)**: 14 C33, 3 G33. By phylum: Chordata 13/10,795; Nematoda 2/352; Porifera 2/143; all other phyla 0. By provenance: RefSeq 8/13,387; INSDC-only 7/11,289; other assemblies 2/1,816; DToL 0/2,091. Restricting to 5-nt or 0-nt overlap: 17/26,784 and 15/21,376. Fisher test not run (pre-registered minimum 3 independent new Tyr origins; there is 1). Independent origins of non-U33 among Stop: 11 (parsimony, NCBI taxonomy).

Candidates (**sufficiency counterexamples to "C33 makes UAA read as Tyr"**, strongest first):

| Species | Clade | Evidence |
|---|---|---|
| *Microdevario kubotai* (+ *M. nanus*) | fish, Danionidae | C33, loop `CCGUAGC`, 5/5 stem, 2 callers. **All 9 confident calls in the genus are C33**, across independent records of *M. kubotai* (isolates kub1–kub4, one study; AP011403, AB937094) and *M. nanus*. Read support: see §Candidate checks |
| *Albula vulpes* (panel), *A. glossodonta*, *Pterothrissus*, plus `AP002973` | Albuliformes | C33 in every confident Albuliformes call (4 of 4 incl. the extra). An **independent *A. vulpes* mitogenome** (`OZ556087.1`) is C33; 469 reads of the record's own run carry C33, 0 carry the U33 k-mer |
| *Serrivomer* ×3, *Stemonidium*, *Coloconger*, *Ophichthus asakusae*, *Xyrias revulsus* | eels | C33; each a single record; independent support not yet checked. 213 other resolved eels are U33 |
| *Lycenchelys verrillii* | fish (Zoarcidae) | C33, 4/5 stem, 2 callers, single record |
| *Toxocara apodemi*; *Globodera ellingtonae* | nematodes | C33; *Globodera* is at the gate minimum (score 10.0, 1 caller) and has a multi-chromosome mitogenome — treat as weak |
| *Craniella zetlandica*, *C. infrequens* | sponges | **G33**, both species, same loop `CGGUAAA` (4/5 stem); the family member *Cinachyrella apion* is U33. Other-assembly provenance |
| *Bufotes pewzowi* | toad | **G33**, 4/5 stem, 1 caller; no independent conspecific record (the only other record is the RefSeq source `MK855100`), and 6 other *Bufotes* confident calls are U33. **Unconfirmed** |

**Reading:** C33 is not sufficient for UAA=Tyr. In Microdevario and Albuliformes it is a fixed, replicated feature of UAA=Stop lineages. Lin et al. reported the substitution "exclusively" in TAA=Tyr lineages among the species they examined. That still holds in molluscs here (0 non-U33 among 1491 resolved UAA=Stop molluscs) but does not extend to Metazoa. C33 remains rare (0.05%), so it could still be a contributing condition.

### T0, T2, T3, T5: new UAA=Tyr
**One new primary-Tyr genome: *Rotylenchulus reniformis* `CM003310.1`** (Nematoda, Hoplolaimidae): 32 internal TAA, 6 genes, 15 at conserved-Tyr columns, LLR 67.4, Codetta Y, plus 3 internal TAG.
- *Independent evidence searched (rule a, b, c):* **(c) satisfied**: in PRJNA1026095 (separate BioProject, 2023 WGS) all 35 unique-k-mer TAA/TAG sites carry the same codon at read depth 70–474 (Tyr alternatives 0–5 per site). **(a) supported**: RNA-seq of PRJEB8325 (low mitochondrial depth) shows UAA at 16 of 32 TAA sites in 55 reads, none reading UAU/UAC. **(b)** no proteomics found. **Independent specimens:** 49 partial *COX1* records (Japan 15, Colombia 10, Mexico 5, Turkey 5, USA 5, Spain 3, Vietnam 2, others) carry TAA at positions 16031, 16163 and 16322 in 47–49 of 49 and at 16235 in 41 of 49 (TAG in the other 8). At 16010 (genome TAA, conserved-Tyr column) the 10 records that cover it read **TAT**, so the column accepts both Tyr codons and TAA across populations. 41 of 49 of those records are annotated **with translation table 14, in which UAA = Tyr**, so other authors appear to have treated this code as unusual; I found no paper that states it.
- *Limits:* **its trnY is not resolved** (stem 3/5, score 12.8; U33, `UUGUAAA`; the panel already reads Radopholus similarly), so under the locked audit it does **not** count as a new clade. The 3 TAG sites are real in the reads (227–434) and one (COX3 col 239) sits in a conserved-Tyr column, so UAG may also be read as sense there, or an alignment artefact — unresolved. Relatives in the family that I could test (*Rotylenchus goodeyi*, *Hoplolaimus galeatus*, *H. columbus*; UNVERIFIED records, post-hoc run of the locked pipeline) are UAA=Stop, so the reassignment is limited to *Rotylenchulus* on current data.
- **Status: single new UAA=Tyr genome with independent evidence, trnY unresolved — reported as a candidate, not as a new clade.** Its U33 (if real) would be a third necessity counterexample after *Radopholus* and *Pseudotanais*.
- *Pratylenchus goodeyi* (UNVERIFIED `ON525593.1`, post-hoc): 35 internal TAA in 7 genes, 9 at Tyr columns, Codetta Y, but LLR −2.7, so only "weak"; unverified record, not audited further.
- T2 precision 0/17, recall undefined (descriptive). T3: 1 new independent origin. B1 not testable (no Tyr species with ≥3 called genes); B2 was not run (exploratory).
- Weak tier: 3 (*Teleopsis* fly, *Pteronarcys* stonefly, *Ornithodoros* tick); all **rejected** by the pre-registered weak audit (0 supported sites each).

### T4. Starfish
All 77 resolved Asteroidea are U33 (1 further unresolved); *Asterias rubens* `LR745847.1`: U33. All 5 echinoderm classes: 0 non-U33 among 323 resolved. **Reading: Zhao-consistent** (≥80% of ≥5 asteroid calls U33); Li's C33 in starfish is not supported by any confident call. Li's own species were not identified or fetched; they may differ.

### T6. Albuliformes and Elopomorpha
Reading **lineage-level** (≥1 other Albuliformes species confident C33): 4 confident Albuliformes, all C33. Anguilliformes 213 U, 7 C, 6 unresolved, 1 stem review; Elopiformes 5 and Notacanthiformes 13, all U. No second trnY locus and 13 intact PCGs in all four. `AP002973` is a **named extra**, excluded from counts, and is *not* the same species as `CM165726` by sequence (87.0% identity; see `notes_albula_records.md`).

### Auxiliary checks (confirmatory)
A1 (U33 in the nearest Stop relatives of Tyr origins): 1 of 1 origin (order-level sisters: 334 U, 2 C). B1: not testable. B3 not run (Gate 4). Post-hoc stem re-inclusion (not confirmatory): **Aves** (1,575 calls, shared U·U mismatches at stem positions 1–2) and **Testudines** (161) meet all fixed criteria; with them, T1 is 18/30,307 (the extra is the crane *Antigone rubicunda*, C33, from stem_audit_review). The criteria were fixed after look 2 and are not blind.

## Scope of the null
- 2,033 genomes in stem_audit_review (mainly birds) and 1,956 unresolved cannot show non-U33 either way; Cnidaria, Ctenophora and Chaetognatha have no resolved trnY.
- Nematoda 89.6%, Porifera 77%, Bryozoa 57%, Acanthocephala 49% resolved.
- **8,189 UNVERIFIED-excluded species and the 5,635 other ineligible INSDC-only records were not tested.** 366 molluscs and 46 nematodes are among them. A post-hoc run of the same pipeline on the 8,189 (about 4 h on 8 workers, no new rules) is the single most useful check and is **recommended, pending your approval**. It would be exploratory.
- 2 module errors (Myxobolus, Enteromyxum: codon O/E, zero codons), descriptive only.

## Candidate checks (post-sweep; corrected read-support code after a strand bug in the first version, which tested the wrong alternative allele)
- *Microdevario kubotai*: record k-mer (C33) 572 reads, U33 k-mer 1 read in WGS run SRR39401262; amplicon run SRR4556446: 10 record reads, 0 U33. The WGS run is a different run and specimen from the RefSeq source (AP011403); specimen identity was not traced further. Together with the 7 independent *M. kubotai* records, a sequencing error is effectively excluded. The single U33 read is consistent with noise or a low-level NUMT.
- *Albula vulpes*: own sequencing run SRR22396674 (same specimen): 469 reads carry C33, 0 the U33 k-mer. This excludes a sequencing error in that specimen, not a NUMT.
- *Bufotes pewzowi*: RAD-seq runs carry neither k-mer (not informative); conspecific mitogenome: none independent. **Unconfirmed.**
- No species tested has a nuclear assembly of its own in NCBI among the three, so NUMT was not searched. PCG intactness (13 genes) and one tRNA locus argue against one.

## Gate 4 (not started)
Plan: targeted mtRF1a in 18 flagged species with sister controls; 6 of the 18 have a nuclear assembly. Estimated **20 assemblies, 19.7 Gb total: about 5.5 GB download, 20 GB unpacked, peak about 26 GB** including indices. Free disk is ~246 GB. **This also needs Microdevario, eel and Albula mtRF1a, for which only nuclear assemblies for *Albula glossodonta* and some sisters exist.**

## Deviations and corrections to note
- Compound-location bug fix (post-lock, adopted), interim looks, watcher, reordering (all above).
- Earlier statements corrected here: the second "Albula glossodonta" record is a different species (AP002973), not a replicate; Microdevario kubotai and *M. nanus* replicate within genus, while Bufotes does not.

## Downloads and resources
Ledger (`data/manifests/download_ledger.tsv`): 28.32 GB cumulative in 428 rows (Gate 3: 17.0 GB, mostly read runs, deleted after use). 5 Gate 2 rows have no checksum; NCBI API responses (conspecific records, summaries) are not metered. Free disk 250.8 GB (floor 30 GB).
