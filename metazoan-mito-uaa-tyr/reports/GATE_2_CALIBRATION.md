# Gate 2: Control benchmark and calibration

> **Read the Addendum (sections A1–A9, at the end) first.** It supersedes §2–§3 trnY numbers, the mtRF1a guard note in §0 and the terminal-stop plans. The text above the Addendum is the state committed as `c8f9f5c`.

Date: 2026-11-23. **Status: partially met — one criterion (trnY extraction failure) is not cleanly met; see §3. Nothing was run on any non-control genome.** `PREREG.md` is drafted but **not committed** (this file will record its hash only once it is locked).

## 0. Control panel (frozen before any caller ran)
`results/gate2_panel.tsv`, sha256 `b169a634…` (360 genomes). The allowlist is enforced in code (`src/panel_guard.py`) for the trnY, UAA and Codetta callers; a test with non-panel OZ576062 was blocked. **The `mtRF1a` scripts (`src/mtrf_*.py`) do not call the guard**: they were restricted only by the species lists written into them. Those lists include control genera/species whose nuclear data are not panel ids (*A. ephippium* nuclear genome, *Cephalodiscus* sp. and *Rhabdopleura annulata* RNA-seq). The guard must be added before Gate 3/4.
- **Positives (20, UAA=Tyr documented):** Lin Figshare: *Anomia* (4 genomes of 3 species: *A. chinensis*, *A. trigonopsis* ×2, *A. ephippium*), *Enigmonia aenigmatica* ×2, *Placuna* ×3, *Plicatula* ×2; Zhao MOESM4 re-assemblies: *A. chinensis*, *P. placenta*, *Pseudotanais* sp. a; GenBank: *Cephalodiscus hodgsoni* MF374802, *Radopholus* FN313571, *Pseudotanais* OZ281197–200.
- **Matched negatives (25):** *Pododesmus* ×3, *Heteranomia* ×5, *Rhabdopleura* ×3, tanaids ×6 (*Arctotanais*, *Sinelobus*, *Neotanais*, *Carpoapseudes*, *Apseudes* ×2), nematodes ×8 (incl. *B. xylophilus*, *Meloidogyne*, *Aphelenchoides*, *C. elegans*).
- **Broad canonical panel (315):** seeded stratified sample of RefSeq (circular, complete, ≥12 CDS, 0 ambiguous bases, tables 2/4/5/9/13), all phyla with quota (70 Arthropoda, 45 Mollusca, 40 Chordata, 30 Nematoda, 25 Platyhelminthes, …). **Excluded on purpose:** Echinodermata, non-*Rhabdopleura* Hemichordata, Anomiida, Tanaidacea, and any record whose annotation uses table 14/33 (so the pre-registered starfish test stays blind).
- **Arctotanais (decision 4):** LC597489 annotation gives 13 CDS, 0 internal TAA/TAG, 6 genes with complete TAA terminators, transl_table 5; Matsushima: TAA=Tyr "not observed in other peracarids". Kept as a negative. (Annotation under table 5 cannot show the alternative; the counter below is an independent test and also calls it Stop.)
- *Plicatula* stem and *A. ephippium* 71-kb genome: OZ576062 (concatemer) was not used; the Figshare *A. ephippium* was.

## 1. UAA caller
`src/uaa_caller.py`: pyhmmer 0.12.3 hmmsearch of the 13 metazoan OXPHOS profiles (Codetta's `examples/mito-db`) against six-frame translation (TAA/TAG→X), i-E≤1e-5, first/last 10 codons masked.
**Method bug found and fixed here:** HMMER's default acceleration filters discarded a 723-bit COX1 domain on the AT-rich *Rhabdopleura* translation, so *Rhabdopleura* and *Radopholus* returned zero genes. Filters are now off (`F1=F2=F3=1`, no bias filter). Any earlier per-genome impressions from the pre-fix code are void.

Results (rule: TAA internal ≥5, in ≥3 genes, ≥3 at conserved-Tyr columns, LLR≥3, AND Codetta Y; Stop = ≤1 internal TAA, 0 conserved-Tyr TAA, Codetta not Y):

| | Tyr | Stop | ambiguous | insufficient |
|---|---|---|---|---|
| Positives (20) | **20** | 0 | 0 | 0 |
| Matched negatives (25) | 0 | **25** | 0 | 0 |
| Broad negatives (315) | **0** | 309 | 1 (*Acanthobdella peledina*, 4 TAAs in 1 gene, Tyr-log-likelihood −14) | 5 |

False-positive Tyr calls on 340 unaltered negatives: **0**. Positives have 20–53 internal TAAs in ≥6 genes; negatives have ≤4 TAAs, in ≤1 gene. The 5 "insufficient" are 2 Ctenophora, 2 Platyhelminthes and 1 Cnidaria with <800 aligned codons or <5 genes.

**Read "20/20" as in-sample, not as a held-out recall.**
- *Distinct genomes:* 21-mer containment ≥0.5 collapses the 20 positives to **15 clusters** (e.g. Lin and Zhao *A. chinensis* 0.95; *Pseudotanais* sp. a and OZ281199 0.99; Lin *Placuna vitream* and Zhao *P. placenta* 0.90, possibly one species under two names) in 7 lineage groups, about 5 independent origins. The 25 matched negatives collapse to 22.
- *Thresholds were chosen after seeing these positives.* The LLR floor (3) sits just under the lowest positives (3.4, the duplicated *Pseudotanais* pair). Internal-TAA floor was lowered from 8 to 5 after seeing the planted series. Four counter rules were compared. Codetta was switched to Lin's thresholds after the default gave 17/20.
- *Leave-one-lineage-out* (thresholds = 0.8 × minimum over the other lineages): Anomia 5/5, Cephalodiscidae 1/1, Enigmonia 2/2, Placunidae 4/4, Plicatulidae 2/2, Radopholus 1/1, **Pseudotanais 2/5**. A new lineage with as weak a signal as *Pseudotanais* could be missed. The 0/340 false-positive figure is less exposed to this (largest negative: 4 TAAs in 1 gene), but the same genomes informed the floor.

### Planted-null sensitivity (1,430 runs on 65 negatives; TAT/TAC→TAA at conserved-Tyr columns; 5 replicates for 5–50%, 1 for 0/100%)
| Fraction | n | Tyr | ambiguous | Stop | insufficient |
|---|---|---|---|---|---|
| 0% | 65 | 0 | 0 | 63 | 2 |
| 5% | 325 | 0 | 306 | 9 | 10 |
| 10% | 325 | **79 (24.3%)** | 235 | 1 | 10 |
| 25% | 325 | **291 (89.5%)** | 24 | 0 | 10 |
| 50% | 325 | **310 (95.4%)** | 5 | 0 | 10 |
| 100% | 65 | **63 (96.9%)** | 0 | 0 | 2 |

Counter-only rule; with the Codetta conjunction (replicate 0 only, 65 genomes per fraction) Tyr calls are 0/65, 16/65, 59/65, 63/65 and 63/65 at 5/10/25/50/100%. Failures at 100% are the two genomes with insufficient aligned coverage. Position-mapping checks: 0 of 139 recorded codons disagreed with the genome; planted sites lost are lost because weak domains fall below the E threshold. **Minimum reliably detectable level: roughly ≥25% of conserved Tyr codons replaced (≈6+ codons).** Below that, the result is "ambiguous", not "Stop".

### Codetta `-m` (custom OXPHOS profile DB)
Default thresholds (1e-10, 0.9999): 17/20 positives Y, 3 '?', 0 of 340 negatives Y. Lin's thresholds (1e-5, 0.7): **20/20 Y, 0 Y among 340 negatives** (1 broad negative called R at *Acanthobdella*). Pfam `-m` (full Pfam-A): 15/20 Y, 5 '?', 25/25 negatives '?'. On the planted series Codetta (Lin thresholds) calls Y at 5% in 55/65, so on its own it is too sensitive to a couple of codons; the counter gate is what rejects stray TAA. KACI was not built (Linux/WSL only).

## 2. trnY extractor (`src/trny_extract.py`)
MITOS2 `trnY.cm` (Zenodo, CC BY 4.0; md5 verified), consensus columns N33–N36 = 31–34, loop 30–36, stem 25–29/37–41. G34–U35–A36 verified, loop exactly 7 nt, stem ≥3/5.
Control results: positives N33 = **C in 14/20, U in 6/20** on the raw calls (Radopholus and all 5 *Pseudotanais*: U; all Anomiidae/Placunidae/Plicatulidae/Cephalodiscus: C). Matched negatives: U in 25/25. Zhao's C31/C34 and Lin's U36 mapping: the extractor gives canonical N33 independent of linear coordinates. Lin's "36" is an alignment/logo column as inferred in Gate 0; not verified here since their alignment is not available.

**Two bugs found and fixed:** (a) a bad option combination (`--nohmm --noF4b`) in the sensitive pass made 60 genomes crash and be labelled "error" (fixed to `--max`); (b) the sensitive pass accepts noise (score 8–10 hits in Cnidaria, which lack trnY), producing five non-U33 calls that would have become false sufficiency counterexamples.

### 3. trnY criterion (Gate 2): "<10% extraction failure across non-vertebrate negative controls" — **not cleanly met**
Final gate (proposed, from split-half validation): U33 accepted if CM score ≥20, or ≥10 with ≥1 independent caller (ARAGORN/tRNAscan-SE) at the same site; any non-U33 needs ≥10, ≥1 independent caller, ≥4/5 stem.
| Denominator | Failure |
|---|---|
| all 300 non-vertebrate negatives (including 27 with no annotated trnY, mainly Cnidaria) | 53/300 = **17.7%** |
| the 273 with annotated or known trnY | 28/273 = **10.3%** (95% CI 7.2–14.4%) |
| as above, excluding Cnidaria (272) | 27/272 = 9.9% |
Split-half check (choose threshold on one half at ≤3% discordance with annotation, score on the other): test failure 9.7% and 11.3%. Thresholds were **not** set to hit 10%.
Failures: Nematoda 11, Acanthocephala 4, Arthropoda 3, Bryozoa 3, Mollusca 2, Platyhelminthes 2 (truncated or divergent trnY, common in nematodes), plus annotated Cnidaria.
Specificity of the gate: accepted but at a different locus than the annotation 6/260 (2.3%); accepted in genomes without an annotated trnY 2/27.
**Non-U33 among 340 negatives after the gate: 1 (*Albula vulpes*, vertebrate, C33, loop CCGUAGC, score 65.4, corroborated).** Five others were rejected as uncorroborated noise.
**Positives:** 16/20 called (13 C33, 3 U33). Not called: *Radopholus* (score 16.5, U33), *Plicatula muricata* (C33, stem 3/5), *Pseudotanais* OZ281198/OZ281200 (scores 14–16, U33). These four keep "no_confident_call" status; raw calls are listed in the data so the known necessity counterexamples are not hidden.
Settings tried for this gate: five (low, high) pairs on the full panel and a 3×5 grid in the split-half; the split halves share the panel used to design the gate, so they are only partly independent. Before that, a plain score cutoff (12 values) and the ungated extractor (6.3% failure, rejected because it accepts noise) were examined. P. placenta and *R. compacta* anticodon stems have 4/5 pairs, the *Pododesmus* stems 3/5.
**Options for the reviewer:** (i) accept the annotated-gene denominator with the stated CI; (ii) accept a looser CM gate for U33 only, with a documented discordance rate; (iii) keep as is and treat trnY as available for ~90% of genomes with the gene. I have not picked one silently; `PREREG.md` §8 lists it.

## 4. mtRF1a (control species only)
Sources: Figshare proteomes (*A. chinensis*, *P. vitream*, *P. muricata*, *H. squamula*), miniprot models from assemblies (*A. ephippium*, *Enigmonia*, *Pododesmus macrochisma*, *Radopholus*), 36 annotated background proteomes from panel species, and reads (below). Paralogs: candidates best-hitting MRPL58/MTRFR and bacterial RF1s removed with a FastTree tree (E. coli RF1/RF2 and human paralogs).
**Contaminant caught:** *Enigmonia* unplaced contigs `JBTMWK…` carry RF1 genes at 0.79–0.89 from *E. coli* RF1 (bacterial); the mitochondrial copy is on `CM152705`.
Human reference positions (Zhao −25 = human 203/216): residues and loop (`EGGVHRVQRVPKTEKQGRVHTSTMTVA`):

| Species | UAA | 203 | 216 | Loop (196–222) | identity to human |
|---|---|---|---|---|---|
| *A. chinensis* | Tyr | Q | **S** | `EGGVHRFQRETRVKQARNFFSCTITIA` | 41% |
| *A. ephippium* | Tyr | Q | **S** | `EGGVHRFQRETRVKKARTFYSCTITIA` | 44% |
| *Enigmonia* | Tyr | Q | **S** | `EGGVQRFQRETRVKVPRVFFSCTITVA` | 41% |
| *P. vitream* | Tyr | Q | **S** | `EGGVQRFQRETRVKSSRTFFSCTITIA` | 37% |
| *Plicatula muricata* | Tyr | Q | T | `EGGVHRVQRKPVTDRNNKRQTSTITVA` | 63% |
| *C. hodgsoni* (3 reads) | Tyr | Q | **S** | `HKVQRIPKTENKNRIHSST` | — |
| *Cephalodiscus* sp. (reads) | Tyr | Q | **S** | `HKVQRIPKTERNNRIHSST` (≥4 independent reads) | — |
| *Radopholus* (tblastn genome, spliced manually) | Tyr | **K** | T | `RLPVTNSKGKGVLHTST` (+2 aa vs human) | — |
| *Rhabdopleura annulata* (5 reads) | Stop | **K** | T | `EGGVHRVKRVPKTERKGRIHTSTSSVA` | — |
| *Heteranomia squamula* (2 copies) | Stop | **I** | T | `EGGVHRVIRVPSTSKADKVHTSTVHVA` | 70% |
| *Pododesmus macrochisma* | Stop | Q | T | `EGGVHRVQRVPVTEQKGRVHTSTITLS` | 78% |
| 36 background species (36 proteomes, 73 copies) | Stop | Q | T | — | 56–96% |

Reads with loops identical to human (`GGVHRVQRVPKTEKQGRVHTSTMTVA`, 3 of 7 in SRR11101524) are human contamination and were excluded. *Rhabdopleura annulata* and *Pseudotanais*: see §4b.
Descriptive summary (no test): among panel UAA=Tyr species with loop data, 6/8 carry S216 and 2 do not (*Plicatula*, *Radopholus*); in all Stop panel species with data, 0 carry S216. **Anomiidae and Placunidae share S216 in a loop whose `PxTE` and `HTST` motifs are both altered, and Plicatulidae (UAA=Tyr) keep T216.** Zhao's "T191S" therefore has support for Anomiidae/Placunidae/Cephalodiscus (and Pseudotanais per Zhao), not for Plicatulidae or Radopholus. These are control-panel facts, not new findings.
### 4b. Read-based taxa
Reads were filtered to diamond blastx E≤1e-8 spanning human 205–220; peptides were translated, not assembled. The first, looser look (E≤1e-3) was mostly spurious (457 reads for SRR11101524) and was not used.
| Run | Taxon | Strict loop reads | Result |
|---|---|---|---|
| SRR1695471 | *C. hodgsoni* | 3 | Q203, S216 |
| SRR11101524 | *Cephalodiscus* sp. | 7, of which 3 human contamination | Q203, S216 (4 reads) |
| SRR11101525 | *Rhabdopleura annulata* (UAA=Stop) | 5 | K203, T216 (Zhao's Q178K); identical in the corrupt and the validated run |
| ERR15194044 | *Pseudotanais* sp. a | 2, neither usable | **untested**: one read is bacterial RF1 (`QGRVHTSACTVAVMPE`), one stops before 216; an unbiased motif scan found one further read (`HRVQRVPVTESSGRIHTS`), also bacterial-like. Zhao needed an assembly of two runs. |

**Integrity problem found and repaired:** the first SRR11101525 file was 52.7 MB too large (two `curl` processes wrote to one file after I swapped jobs). Its result was discarded; the run was repeated with a server size match and `vdb-validate` (rc 0). SRR11101524 and ERR15194044 sizes match the server; they were not vdb-validated.

## 5. Sources and resources in Gate 2
New downloads: Lin Figshare 6 files (mitogenomes.zip 341 KB + 4 proteomes; MD5 match), Zhao MOESM9 (public) and MOESM4 (user-supplied), MITOS2 refdata (19.6 MB, md5 verified), nuclear assemblies 0.98 GB, background proteomes, SRA reads 6.6 GB (SRR11101524, SRR11101525 twice, ERR15194044) plus SRR1695471, pyhmmer/diamond/FastTree/SPAdes. Cumulative retained or processed downloads since Gate 0 ≈10.6 GB (<20 GB; API responses and aborted ENA transfers not metered). No job exceeded 1 hour. Machine-level memory: one of my orphaned workers briefly held ~9 GB swap and was killed; `diamond` is now run with lower block size.
Failures logged: two bugs above; tRNAscan-SE prompt hang (`-o /dev/stdout`; the discarded run produced no output), ENA slow/dropped transfers (switched to the NCBI SRA open-data bucket), and the corrupt SRR11101525 download (§4b).

## 6. Gate 2 acceptance summary
| Criterion | Status |
|---|---|
| Documented code recovered for all positives and negatives | **met in-sample** (20/20 genomes = 15 distinct; 25/25 = 22 distinct; 0 FP of 340); leave-one-lineage-out misses 3 of 5 *Pseudotanais* |
| Planted-null sensitivity curve quantified | **met** (table) |
| trnY G34–U35–A36 anchoring, <10% failure | **not cleanly met** (17.7% / 10.3%) |
| mtRF1a paralog separation | **met** for proteome/assembly controls (tree placement; FastTree support values not evaluated); *Pseudotanais*, *Placuna placenta* and *Pododesmus* other species untested |
| Controls only used | **met**; enforced in code for trnY/UAA/Codetta, by species list only for mtRF1a (§0) |


---
# Addendum, 2026-11-23 (after the Gate 2 decisions)

## A1. trnY diagnostic pass (bounded; 15:44–16:15, 31 min of the 2 h allowed)
| Variant | Known-gene failure | Accepted non-U33 on negatives | Verdict |
|---|---|---|---|
| baseline (checkpoint) | 28/273 = 10.3% | 1 (*Albula*) | — |
| **A: loop-check bug fix** | 26/273 = 9.5% | 1 (*Albula*) | **adopted** |
| B: A + glocal search | 28/273 | 1 (*Albula*) | rejected: loses one nematode, 14 raw non-U33 hits against 6 |
| C: armless nematode models (leave-one-order-out; 34 panel nematodes) | — | 0 at score ≥20 | **rejected**: on the 11 nematodes failing A it recovers 0 at the right locus and 8 at a wrong locus (score ≥20); 18–46 wrong-locus calls on negatives per model |
- The bug: gap columns created by *other* hits' insertions made a valid 7-nt loop look broken. The patched extractor gives results identical to variant A on all 360 genomes.
- C passes the literal rule (no added false non-U33) and was still rejected, because wrong-locus U33 calls are not recoveries. Only one armless construction was tried.
- **A adds a wrong call:** of the two "recovered" genomes, *Lima* is at the annotated locus and `NC_071830` is 7 kb from its annotated gene (score 20.2).
- **The 10.3% accepted as a borderline pass omitted accepted calls at the wrong locus.** Counting failed or wrong-locus: baseline 34/273 = 12.5%; with A 33/273 = **12.1% (95% CI 8.7–16.5)**. Wrong-locus means no overlap with the RefSeq-annotated trnY; some annotations may themselves be wrong.
- Nematode failure is unchanged at 11/38 (29%); 10 of the 11 are Enoplea; their annotated trnY are 51–62 nt.

Per-phylum failure with the adopted extractor (negatives with an annotated/known gene): `results/gate2/trny_phylum_failure.tsv`, reproduced in `PREREG.md`.

## A2. Auxiliary modules (built and tested on the panel)
| Module | Test | Result | Use in PREREG |
|---|---|---|---|
| Origin counter (`src/origins.py`) | 6 hand cases; 400 random trees by exhaustive enumeration | 0 mismatches for minimum changes and hard-polytomy gains; `gains_min` checked on hand cases only | locked |
| Terminal stops (`src/aux_checks.py`) | against annotation, negatives | "first in-frame stop" rule: 567/567 TAG and 1595/1596 TAA right. **TAG-first rule calls 43.3% (409/944) of annotated TAA stops TAG** | "first" locked, conservative; TAG-first is an upper bound only |
| – incomplete T/TA stops | against annotation, broad panel | 161/443 (36%) recognised | not testable; exploratory |
| – coverage | | 41% (matched) and 64% (broad) of annotated genes called; the rest not found or not aligned to the 3' end | reported with every result |
| Codon O/E | *A. chinensis* vs Zhao Table S4 | TAA/TAT/TAC 0.31/0.18/0.51 (83 codons, aligned regions) vs 0.34/0.22/0.44 (whole CDS) | descriptive |
| Extra tRNA (`src/extra_trna.py`) | Zhao's two claims; false positives | finds *Radopholus* UUA and *Pseudotanais* AUA **with ARAGORN only**; same class of hit in 12/315 (UUA) and 38/315 (AUA) broad negatives. Two-caller hits: 0/20 positives, 1/25 matched (*Heteranomia aculeata* M, UUA), 1/315 broad (AUA) | exploratory |

**A circularity was found and removed.** The first terminal-stop version looked for TAG first in Tyr genomes; it agreed with annotation in 69/69 genes only because both assume the code. Both rules used now treat every genome identically.

**Seen on the panel, therefore not confirmatory:** TAG share under the "first" rule — Anomioidea (Tyr) 36/57, Plicatulidae 8/12, Pseudotanais 15/16, Cephalodiscidae 1/4, Radopholus 1/5; canonical sisters *Pododesmus* 9/15, *Rhabdopleura* 5/8, *Heteranomia* 1/17, Tanaidacea 13/35, nematodes 10/39. *Pododesmus* and *Rhabdopleura* are TAG-rich without reassignment.

## A3. mtRF1a firewall
`src/mtrf_guard.py` is now called by every `mtRF1a` script that reads data (`mtrf_bg_fetch.py` uses `panel_guard`); tested: `SRR11101525` and an unknown run are refused, panel sources run. miniprot was run from the command line in Gate 2 and must go through a guarded script from now on.

## A4. Deviation log: data used that is not in the frozen panel
| Data | Taxon | When | Why it is out of panel |
|---|---|---|---|
| SRR11101525 (mtRF1a) | *Rhabdopleura annulata* | Gate 2 | genus is a control; species/run not a panel id |
| SRR11101524 (mtRF1a) | *Cephalodiscus* sp. | Gate 2 | same |
| GCA_990049579.1 (mtRF1a) | *Anomia ephippium* | Gate 2 | panel species, different specimen and assembly from the panel mitogenome |
| OZ576062.1 (trnY only; tRNAscan-SE, ARAGORN) | *Anomia ephippium* | Gate 0 | scanned before the panel was frozen; excluded from it afterwards |
| Figshare *Isomonia*, *Idas*, *Lima*, *Spondylus*, *Ctenoides*, *Planostrea*, *Terua*, *Leptopecten* | bivalves | Gate 2 | record lengths, topology and gene counts listed; no caller run |
The first two were mtRF1a-caller runs outside the frozen panel without the guard. They are closed to reuse in code and their results are descriptive only.

## A5. *Plicatula muricata* mtRF1a (UAA=Tyr with T216 retained)
- **Loop** (human 196–222): `EGGVHRVQRKPVTDRNNKRQTSTITVA`; human `EGGVHRVQRVPKTEKQGRVHTSTMTVA`. No indels, 63% identity. Q203 = `CAG`; **T216 = `ACA`**. Changed: V205K, K207V, E209D, K210R, Q211N, G212N, R213K, V214R, **H215Q**, M219I, so both the `PxTE` and the `HTS` motifs are altered although T216 is kept.
- **Gene model:** `PMU02G001039`, chr02:78305761–78334157 (−), 5 exons. CDS rebuilt from genomic DNA matches the deposited protein exactly (276 codons, no internal stop); the loop lies inside one exon. The model has 276 residues against 380 in human, so the N-terminus looks incomplete. GGQ motif intact (`GAGGQH`).
- **Paralog audit:** best human hit MTRF1L (292 bits) over MTRF1 (235); the proteome's other RF-like genes are MRPL58 (`PMU03G000913`) and MTRFR (`PMU06G002011`). In the FastTree tree it groups with the two *Heteranomia* copies (support 0.917), in a clade of RF1-type proteins (0.612). `tblastn` against the assembly (E≤1e-3): 17 HSPs, all at this locus; one covers the loop; several are partial tandem copies of the exon at residues 129–174.
- **Limits:** one individual, one assembly; no transcript evidence read; *P. australis* has no nuclear data. *P. plicata* RNA-seq (SRR3350473) is out of panel and is a pre-registered Gate 4 target.

## A6. Leave-one-lineage-out (`results/gate2/lolo.tsv`)
Every lineage is fully recovered by every procedure except *Pseudotanais*: fixed rule 5/5, negatives-only thresholds 5/5, midpoint thresholds 3/5, 0.8×min 2/5. Holding out whole origins gives the same. 0 false positives among 340 negatives in every fold. *Cephalodiscidae* and *Radopholus* rest on one genome each. The fixed rule was hand-set after seeing all lineages, so its 5/5 is in-sample.

## A7. Weak tier (panel)
On negatives: 1/340 (*Acanthobdella*). Planted series (65 genomes, replicate 0; primary / weak / neither): 5% 0/59/6; 10% 16/48/1; 25% 59/6/0; 50% 63/2/0; 100% 63/2/0; 0% 0/0/65.

## A8. *Albula vulpes* (calibration observation)
`NC_082992.1` (from OR546140): trnY 5488–5418 (−), loop `CCGUAGC`, N33 = **C**, stem 5/5, CM score 65.4; ARAGORN (G34 at 5456) and tRNAscan-SE (51.8) agree; overlaps the annotated gene. UAA = Stop: 0 internal TAA over 3,376 aligned codons, 13 genes, 72 conserved-Tyr columns; table 2. A candidate sufficiency counterexample seen during calibration, not a discovery. Follow-up is pre-registered (T6). OR546140 is the only *A. vulpes* mitogenome in the INSDC list, so no conspecific mitogenome is available there; two *A. glossodonta* records exist (CM165726.1, AP002973.1).

## A9. Resources
Resource rule of 2026-11-23 applies (≥30 GB free; ledger). **No ledger was kept during Gates 0–2.** `data/manifests/download_ledger.tsv` was reconstructed from files on disk and session logs; 5 SRA/ENA rows have no checksum because the files were deleted, and NCBI/ENA API responses were never metered. Cumulative: **11.33 GB** in 277 rows (5 without checksum), including the *Plicatula* genome and GFF fetched for A5. Earlier totals in this file (≈10.6 GB) were estimates.

## A10. Why 1 of 340 negatives is in the weak tier
*Acanthobdella peledina* `NC_063475.1` (leech; mitogenome "complete", full length). All four in-frame TAA codons lie in repeated, degraded ATP6-like copies downstream of the annotated gene (annotated ATP6 ends at 14,592; copies between about 14,660 and 18,330, 11 partial profile domains in several frames). Three sit at an Arg column (p(Arg) = 0.79) and one at a Tyr column (p(Tyr) = 0.74). The full-length ATP6 (13,893–14,592) has none. The genome is annotated table 5; the caller returns **ambiguous** (weak tier), not Stop. With 4 internal TAAs in one gene it fails the primary rule (≥5, in ≥3 genes), and ambiguous genomes do not enter the primary tests. It is a pseudogene/duplicate-copy artefact, and is the likely reason for Codetta's Arg call in that genome (inferred, not tested). The caller does not distinguish duplicate gene copies; this class is the kind the weak-tier audit must reject.

## A11. Threshold method for the primary UAA caller
Definition: every count threshold is the maximum over the 340 negative controls plus a margin, set without any positive or non-panel genome: internal TAA ≥ 4+1 = 5; genes with TAA ≥ 1+2 = 3; conserved-Tyr-column TAA ≥ 1+2 = 3; Tyr log-likelihood ≥ 3 (no negative exceeds 0). Codetta thresholds are Lin's published ones; the Codetta setting was preferred over the default after seeing 17/20 vs 20/20 positives, so that choice did use positives. The margin values were chosen while the positives were visible and the formula was written afterwards to describe them; it reproduces the fixed rule exactly. The choice of the "most conservative on negatives" principle is the author's.
- Specificity check, leave-one-phylum-out over the negatives (thresholds from the other phyla, both with and without the margin): **0 false Tyr** in every held-out phylum.
- Sensitivity, leave-one-lineage-out: thresholds of this method do not use positives, so recall equals the in-sample 20/20 (every lineage in A6). Two other procedures that do use positives give *Pseudotanais* 3/5 (midpoint) and 2/5 (0.8×min); other lineages full. The in-sample figures are not held-out estimates for a new lineage.
- Because 0/340 has a 95% upper limit of 1.1% (Wilson), a sweep of about 33,000 genomes could still hold up to ~370 false Tyr calls from the error rate alone. The clade rule and the locus/alignment audits, not the caller, carry the specificity at sweep scale.

## A12. Locus-identity audit (`src/locus_audit.py`; panel test)
Rules: verified GUA in a 7-nt loop with stem ≥4/5; overlap ≤10 nt with any annotated CDS or other tRNA (if no annotation, with profile PCG domains and non-GUA tRNA calls); if a trnY is annotated, >50% overlap on the same strand. Of 305 confident calls on the panel, 261 pass. Failures: stem <4/5 18, disagreement with the annotated trnY 6, overlaps the rest. *Albula* passes. The 10-nt tolerance and the stem cut-off are author choices; with tolerance 0 the pass count is 188/305, 5 nt 254, 10 nt 261, 20 nt 264; 73 calls pass only because of the 10-nt tolerance, because neighbouring mitochondrial genes normally overlap by a few nucleotides. Two *Enigmonia* calls fail because ARAGORN reports an overlapping UUG tRNA, a caller artefact rather than a wrong locus.

## A13. Final pre-lock checks (2026-11-23)
- **Stem rule and G·U.** The extractor already counted G·U as paired. Of the 18 confident panel calls failing the ≥4/5 rule, **0 pass with G·U counted** and 0 with Watson–Crick only. (If only Watson–Crick pairs counted, 42 of the 305 confident calls would fail.) 12 of the 18 share the arm `UUAGG…CCUUU` (including the birds *Trogon*, *Apteryx*, *Gavia* and *Rhynochetos*) and three are the *Pododesmus* arm. An independent recount of the contiguous arm disagrees with the recorded pair count in 4 of 305 calls (the model alignment allows insertions there).
- **Gate 3 runner equivalence.** `src/gate3_run.py` (Codetta on every genome, no pre-filter) was run on the 360 panel genomes: trnY, second callers, tRNA lists, UAA counter, Codetta code string, terminal-stop ("first", "tagfirst") and codon O/E are **identical for all 360** (`src/gate3_equiv.py`); 0 module errors; median 12.7 s per genome, maximum 107 s. Codetta pre-filtering was not adopted, so no skip needs justifying.
- **Classifier on the panel** (`src/gate3_classify.py`): UAA 20 Tyr / 25 Stop / 309 Stop, 5 insufficient, 1 weak; locus audit passes 188/254/261 at 0/5/10 nt; pass/fail identical to the Gate 2 audit. 17 genomes are stem_audit_review: one of the 18 stem failures also fails another check and is unresolved.
- **Sweep set:** 33,080 representatives, 32,982 included; excluded: 94 ambiguous bases >1%, 3 length, 1 by PREREG (`OZ576062`). Roles among included: 32,624 confirmatory, 334 panel, 24 panel species. All 13,199 GenBank records requested for Gate 3 were retrieved.

## A14. PREREG lock
- PREREG.md commit: `5aa5774eb2ee906a9838a578bf4ccbe1d012bb92` (2026-11-23T16:51:04-05:00)
- PREREG.md SHA256: `e23bc13f37ca48647b86cdc6b3d70abc5f7e0b69c71968de05e7b32285b4fa0c`
- Gate 3 code is refused unless PREREG.md is committed, unmodified, and its SHA256 appears in this file (`src/panel_guard.py`).

## A15. Post-lock bug fix: compound annotation locations (2026-11-23)
**Found after the lock and after interim look 2** (`results/gate3/interim/LOG.md`), while checking why strong calls were rejected. Origin-spanning and spliced CDS/tRNA features were stored as min..max, so a gene wrapping the origin covered the whole genome and every trnY call in that record was rejected as wrong-locus. The rule in PREREG §2.3 is unchanged; the code did not implement it. Fix: intervals are the parts of the location (`src/gate3_features_fix.py`, `src/panel_features.py`; classifier default `data/sweep/features_fixed.jsonl`). The original `data/sweep/features.jsonl` is kept for the without-fix sensitivity run. PREREG.md is not edited.
- Scope: 675 of 28,524 annotated sweep records had a >5 kb interval; 4 after the fix.
- **A12 recomputed on the panel:** confident calls 305; locus audit passes at 0/5/10 nt **193 / 263 / 271** (were 188 / 254 / 261; 10 calls change, all broad negatives); 17 stem_audit_review unchanged. PREREG §2.3 quotes the older numbers.
- The trnY failure percentages (A1) do not change: 26/273 unresolved, 33/273 (12.1%, CI 8.7–16.5) with wrong-locus calls. Three panel records had whole-genome trnY intervals, which could only have hidden errors; recomputing with parts gives the same 7 wrong-locus calls.
- Effect on the 6,743-genome interim snapshot (no final results were used to choose the fix): 89 unresolved become resolved and 2 become stem_audit_review; the non-U33 Stop candidates are unchanged.
- Final Gate 3 results are reported **with the fix (primary) and without it (sensitivity)**.

