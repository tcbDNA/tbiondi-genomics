# PREREG (LOCKED by the commit that adds this file; 2026-11-23)

No trnY, UAA or mtRF1a caller ran on a non-control genome before this file was committed (deviations: GATE_2 A4). The commit hash and SHA256 are in `reports/GATE_2_CALIBRATION.md`. Evidence for every number: `reports/GATE_2_CALIBRATION.md` (sections and Addendum A1–A9).

## 0. Disclosed before locking (seen on the control panel or published; not discoveries)
- U33C is not necessary: *Radopholus* and *Pseudotanais* are UAA=Tyr with U33 (Zhao; panel agrees).
- All C33 calls in UAA=Tyr panel genomes are in Anomia, Enigmonia, Placunidae, Plicatulidae and *Cephalodiscus*.
- *Albula vulpes* `NC_082992.1` (UAA=Stop) has a corroborated C33 (T6).
- Li 2019 reports C33 in starfish trnY; Zhao shows U33 in *Asterias* (T4).
- mtRF1a, human `NP_061914.3` coordinates (Zhao position + 25; Zhao's text "R178" is a typo for Q178): S216 in *Anomia chinensis*, *A. ephippium*, *Enigmonia*, *Placuna vitream*, *Cephalodiscus* (2 taxa); **T216 retained in *Plicatula muricata* (loop otherwise remodelled) and *Radopholus*, both UAA=Tyr**; K203 in *Radopholus* and *Rhabdopleura annulata* (UAA=Stop); *Heteranomia* I203; 36 background species Q203/T216.
- Terminal-stop shares and codon O/E for the known reassigned lineages and their sisters (A2).
- **Consequence:** all known reassigned lineages are in the calibration panel. Anything reported for them is descriptive, in a separate section. **Confirmatory inference applies only to non-panel genomes.**

## 1. Data
- Sweep: RefSeq r238 ∪ NCBI organelle assemblies (20,195 species-level representatives) ∪ 12,884 eligible INSDC-only species (title "complete" or circular topology; not UNVERIFIED, not partial). One record per species, RefSeq > circular > longest; open-nomenclature taxa kept separate. Provenance = RefSeq / DToL / other assembly / INSDC_only; every table is stratified by it.
- Excluded from confirmatory tests: the species of the 360 panel genomes. Excluded from UAA calls: `PP708078`, `PP708091` (UAA codons removed, per Zhao), `OZ576062` (apparent concatemer). Excluded everywhere: length outside 8–200 kb or >1% ambiguous bases. All exclusions are written to a TSV with a reason.

## 2. Calls
### 2.1 UAA
Counter: six-frame translation (table 5, TAA/TAG→X), 13 OXPHOS profiles, HMMER filters off, i-E ≤1e-5, first/last 10 codons of each gene masked; conserved-Tyr column = p(Y) ≥0.5. Codetta `863359e`, `-m`, OXPHOS profiles, `-e 1e-5 -r 0.7`.
- **insufficient:** aligned codons <800, or genes <5, or conserved-Tyr columns <15.
- **Tyr (primary, high specificity):** internal TAA ≥5, in ≥3 genes, ≥3 at conserved-Tyr columns, Tyr log-likelihood ≥3, and Codetta TAA=Y.
- **Stop:** internal TAA ≤1, none at conserved-Tyr columns, Codetta not Y.
Thresholds above are set from the 340 negative controls only: each is the negative maximum plus a margin (internal TAA 4+1, genes 1+2, conserved-Tyr TAA 1+2, LLR ≥3). No positive and no non-panel genome enters the thresholds. The margin was chosen with the positives visible and is the author's choice (GATE_2 A11).
- **weak candidate:** not primary, and (internal TAA ≥2 with ≥1 at a conserved-Tyr column, or Codetta Y).
- **ambiguous:** everything else.

**Weak-tier audit (never counted as Tyr in any primary test).** A weak candidate is promoted to "supported weak" only if (i) ≥2 genes have TAA at columns where ≥80% of the closest `Stop` relatives have Tyr — the 10 closest; if fewer than 10 exist, the available ones with a minimum of 5; with fewer than 5, the case is flagged for manual conserved-alignment audit and is not promoted automatically — in an alignment free of frameshifts and internal TAG within 30 codons, **and** (ii) another species in the same family is primary, or is weak with TAA at ≥1 of the same columns. Supported-weak results appear only in a labelled sensitivity analysis.

### 2.2 Measured performance of the UAA call
| Measure | Result |
|---|---|
| In-sample recall | 20/20 genomes = 15 distinct (21-mer clusters), about 5 origins |
| False Tyr on negatives | 0/340 (309 + 25 Stop, 1 ambiguous, 5 insufficient) |
| Planted-null recall, primary | 0 / 24 / 90 / 95 / 97% at 5 / 10 / 25 / 50 / 100% of conserved Tyr codons replaced |
| Planted, weak tier in addition (replicate 0) | 59/65, 48/65, 6/65, 2/65, 2/65 |
| Weak tier on negatives | 1/340 (*Acanthobdella*) |

Leave-one-lineage-out recall (thresholds and Codetta setting set without the lineage):
| Lineage (genomes; distinct) | Fixed rule (in-sample) | Negatives-only | Midpoint | 0.8 × min |
|---|---|---|---|---|
| Anomia (5; 3) | 5 | 5 | 5 | 5 |
| Enigmonia (2; 1) | 2 | 2 | 2 | 2 |
| Placunidae (4; 3) | 4 | 4 | 4 | 4 |
| Plicatulidae (2; 2) | 2 | 2 | 2 | 2 |
| Cephalodiscidae (1; 1) | 1 | 1 | 1 | 1 |
| Radopholus (1; 1) | 1 | 1 | 1 | 1 |
| **Pseudotanais (5; 4)** | 5 | 5 | **3** | **2** |

0 false positives in every fold (fixed rule, negatives-only and the two positive-informed procedures); holding out whole origins gives the same. The fixed rule was hand-set after seeing all positives, and four counter rules and two Codetta settings were compared. **A weakly reassigned lineage can be missed; absence of a Tyr call is not evidence of absence.**

### 2.3 trnY position 33
Infernal 1.1.5 with MITOS2 refseq89m `trnY.cm`; GUA at canonical 34–36, 7-nt loop (only the sequence's own insertions count), stem ≥3/5; ARAGORN 1.2.41 and tRNAscan-SE 2.0.13 as corroborating callers.
- **U33:** CM score ≥20, or ≥10 with ≥1 corroborating caller. **Non-U33:** ≥10 and ≥1 corroborating caller (stem rule below).
- **Locus-identity audit (applied to every call, U33 and non-U33 alike; required before any genome counts as a counterexample or new clade):** verified GUA in a 7-nt loop; **≥4 of 5 anticodon-stem pairs, with G·U wobble counted as paired**; overlap ≤10 nt with any annotated CDS or other tRNA (records without annotation: profile PCG domains and non-GUA tRNA calls); and, if the record annotates a trnY, >50% overlap with it on the same strand. Pass counts are reported at 0, 5 and 10 nt (panel: 188, 254 and 261 of 305). G·U was already counted when the panel numbers were made: of the 18 panel calls failing the stem rule, **0 pass with G·U counted**.
- **Statuses:** *resolved* (passes everything); *stem_audit_review* (passes the gate and every other check, fails only the stem rule; listed explicitly, never silently dropped, and not counted as resolved); *unresolved* (no call, gate failure, or wrong locus). Unresolved and stem_audit_review are never imputed. A rate is always over *resolved* genomes, with the other two counts beside it.
- **Every main table carries the phylum-stratified trnY completeness** (resolved, unresolved, wrong-locus), with the caveat that Nematoda (29% unresolved), Acanthocephala (67%) and Bryozoa (50%) are under-sampled.
- Unresolved genomes carry a reason and are reported by phylum. No looser pass will be run. Every non-U33 call gets a structure audit.
- Failure on panel negatives with an annotated or known gene (non-vertebrate): **26/273 = 9.5% (CI 6.6–13.6)** after a bug fix (10.3% before). Adding accepted calls at a locus other than the annotated one: **33/273 = 12.1% (CI 8.7–16.5)**. Counting all 300, including genomes with no trnY: 51/300 = 17.0%. Treated as a **documented borderline pass**, not a clean one.

| Phylum | Unresolved / known gene (95% CI) | Accepted at wrong locus |
|---|---|---|
| Arthropoda | 2/74, 2.7% (1–9) | 2 |
| Mollusca | 1/52, 1.9% (0–10) | 2 |
| Chordata | 0/40 (0–9) | 0 |
| **Nematoda** | 11/38, **28.9%** (17–45) | 1 |
| Platyhelminthes | 2/25, 8.0% (2–25) | 0 |
| Annelida | 0/25 (0–13) | 0 |
| Porifera | 0/13 (0–23) | 0 |
| Nemertea | 0/8 (0–32) | 0 |
| **Acanthocephala** | 4/6, 66.7% (30–90) | 0 |
| **Bryozoa** | 3/6, 50.0% (19–81) | 2 |
| Brachiopoda, Xenacoelomorpha | 0/4 each (0–49) | 0 |
| Onychophora 1/2; Cnidaria 1/1; Rotifera 1/1; all other phyla 0 of ≤3 | too few to estimate | 0 |

Echinodermata and Anomiida are not in this table by design (blind). Two armless/glocal alternatives were tried and rejected (A1). Conclusions about U33/non-U33 in nematodes (mainly Enoplea), acanthocephalans and bryozoans will be stated as under-sampled.

### 2.4 mtRF1a (targeted)
Taxa: reassigned lineages, matched sisters, background panel, plus any new Tyr clade and ≥2 sisters, plus *Plicatula plicata* (SRR3350473). Sources in order: annotated proteins/OrthoDB, TSA, assemblies (miniprot), raw reads; otherwise "untested". Anchored on `VPxTE`/`GRxHTS`; reported at human 203 and 216. Scored separately: (a) **T216S**; (b) **loop remodelling** (196–222): indel count versus human and identity to the nearest canonical sister. Paralogs (MTRF1, MRPL58, MTRFR) and bacterial RF1 excluded by tree. Reads: E ≤1e-8 spanning 205–220, human-identical loops discarded, ≥3 distinct reads. Guards in `src/mtrf_guard.py`.

## 3. Outcome definitions
- **New reassigned clade:** outside Anomioidea, Plicatulidae, Cephalodiscidae, *Radopholus* and *Pseudotanais*, and either two primary-Tyr species in the same clade (same family or lower), or one primary-Tyr genome plus independent evidence: (a) RNA-seq/TSA reads confirming internal UAA, (b) proteomics, or (c) an independent assembly or specimen from a separate BioProject. Otherwise **"unconfirmed single-genome candidate"**.
- **Sufficiency counterexample:** non-panel species called `Stop` with a corroborated, structure-audited non-U33.
- **Necessity counterexample:** species called `Tyr` with confident U33 (known ones are listed, not claimed).
- **Null result:** no new reassigned clade and no sufficiency counterexample among non-panel genomes, reported with the number examined and the Wilson 95% upper bound.

## 4. Named tests (exact; two-sided α = 0.05; both species and origin counts always shown)
- **T0 necessity (descriptive):** Tyr species and origins by N33 state.
- **T1 specificity:** among non-panel `Stop` species with a confident N33, the rate of non-U33 (Wilson CI, by base, by phylum). **The Fisher exact test (Tyr × non-U33) runs only with ≥3 independent new origins (`gains_min` on confirmatory genomes); otherwise descriptive.** Author choice, fixed before Gate 3. Counts and rates are always reported.
- **T2 predictive power:** precision and recall of non-U33 for `Tyr` in non-panel genomes. **Inferential only with ≥3 independent new origins and ≥5 new Tyr species (so not five species of one genus); otherwise descriptive.** Author choice, fixed before Gate 3. Counts and rates are always reported.
- **T3 routes:** independent origins per route (C33; other non-U33; U33), with `gains_min`–`gains_max` from the NCBI taxonomy.
- **T4 starfish/echinoderms:** N33 by class. *Li-consistent* if ≥2 asteroid species in ≥2 families are confident C33 with UAA=Stop; *Zhao-consistent* if ≥80% of ≥5 confident asteroid calls are U33; *inconclusive* if <2 confident calls; otherwise *mixed*. *Asterias rubens* `LR745847.1` is reported by name.
- **T5 replication** of each new Tyr clade (§3).
- **T6 Albula / Elopomorpha:** N33 in all Elopomorpha in the sweep plus both *A. glossodonta* records (`CM165726.1`, `AP002973.1`) and *Pterothrissus* (`NC_005796.1`). Checks: **NUMT** (circularity, 13 intact PCGs, no second trnY); **sequencing error** (read support if the source BioProject has reads); **conspecifics** (no other *A. vulpes* mitogenome in INSDC; any partial trnY sequence). Readings: *lineage-level* if ≥1 other Albuliformes species or an independent specimen is confident C33; *singleton* if all other confident Albuliformes are U33; *artefact suspected* if a NUMT or read check fails; *inconclusive* otherwise.

## 5. Confirmatory auxiliary checks (non-panel genomes only)
| Prediction | Computation | Status |
|---|---|---|
| A1: U33 in the nearest `Stop` relatives of each Tyr origin | lowest shared taxon with ≥1 confident Stop species | confirmatory; exact binomial over origins |
| B1: higher terminal TAG share in Tyr lineages | "first in-frame stop" rule, identical for every genome, against the nearest `Stop` sister taxon with ≥3 called genes; sign test over origins; genes not pooled across origins | confirmatory; conservative. TAG-first reported as an upper bound only |
| B3: T216S / loop remodelling by origin | §2.4; Fisher over origins only if ≥3 origins with data per group | confirmatory only if that minimum is met, else descriptive |
| A2: TAA share of Tyr codons versus third-position composition | `codon_oe` on aligned codons; Spearman across Tyr species | descriptive |
| B4: `Stop` species with standalone T216S versus T216 sisters | "first" rule | descriptive |

Coverage (genes called) accompanies B1 (41–64% of annotated genes on the panel). *Pododesmus* and *Rhabdopleura* are TAG-rich without reassignment; effects are judged against sisters, never the metazoan mean. **Calibration controls (the 360 panel genomes) are reported in a separate descriptive section and never enter confirmatory counts.**

## 5b. Exploratory section (not confirmatory, no p-values reported as tests)
- Incomplete T/TA stop loss (36% detection on the panel).
- Backup TAG downstream and 3' gene overlap (not tested on the panel).
- Extra tRNAs (UUA/CUA/AUA, duplicate GUA): single-caller noise; two-caller hits listed only.

## 6. Statistical plan
Units: species and independent origins. Origins from `src/origins.py` on the NCBI taxonomy: minimum changes (Hartigan) and no-reversal gains, as a min–max range (polytomies soft/hard). Only `Tyr` and `Stop` species enter; `weak`, `ambiguous`, `insufficient` and `unresolved` are excluded and tabulated. With ≤10 expected origins, power is low, so origin-level results are primarily counts and exact intervals. Number of tests run is reported; no correction is applied across T0–T6 because each is a separately named claim.

## 7. Conduct
- Locked rules are not changed after commit; anything else goes under "Post-Hoc Exploratory Analyses", dated.
- Gate 3 runs only with `GATE3_UNLOCK=1`, which the code accepts only if `PREREG.md` is committed.
- Resource rule of 2026-11-23: ≥30 GB free disk; every download in `data/manifests/download_ledger.tsv`; totals in each gate report; long jobs in the background with a progress log.
- Cold review (Gate 5) before `REPORT.md`, as in the task brief.

## 8. Limitations fixed in advance
1. **The UAA thresholds were set by hand with the positives in view.** They are written as "negative maximum plus a margin" (5 / 3 / 3 / LLR ≥3) and no positive enters the formula, but the margins and the Codetta setting were chosen after seeing the positives (GATE_2 A11). Leave-one-lineage-out recall under procedures that do use positives drops to 3/5 and 2/5 for *Pseudotanais*.
2. **False calls at sweep scale:** 0/340 negatives has a 95% upper bound of about 1.1%, i.e. up to ~370 false Tyr calls among ~33,000 genomes. A Tyr call is a candidate until the clade rule, alignments and locus audit support it.
3. **Sensitivity:** recall is 24% when 10% of conserved Tyr codons are replaced; weakly reassigned lineages can be missed, and absence of a Tyr call is not evidence of absence.
4. **trnY:** borderline criterion (9.5% unresolved, 12.1% with wrong-locus calls); Nematoda, Acanthocephala and Bryozoa under-sampled; wrong-locus is judged against annotations that can themselves be wrong.
5. **Calibration coverage:** every known reassigned lineage is in the panel; confirmatory evidence can come only from new genomes. Positives are 15 distinct genomes, about 5 origins.
6. **Taxonomy-based origins:** NCBI taxonomy stands in for a phylogeny; origin counts are ranges.
7. **Author choices:** 10-nt overlap tolerance, stem ≥4/5, weak-tier cut-offs, test minimums.

## 9. Gate 3 execution
Uniform pipeline (`src/gate3_run.py`): every caller, including Codetta, on every genome; **no pre-filter**. Classification and tests: `src/gate3_classify.py`, `src/gate3_build.py`, committed with this file. Panel equivalence with the Gate 2 outputs is in GATE_2 A13.
