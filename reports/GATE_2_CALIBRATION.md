# Gate 2: Control benchmark and calibration

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
