# Gate 4: targeted nuclear mtRF1a

Date: 2026-11-24. PREREG.md (`5aa5774`) and the committed Gate 3 outputs are unchanged. The plan and source list were printed before downloading (`reports/GATE_4_PLAN.md`). Numbering: human `NP_061914.3` positions 203 and 216, loop 196–222. Raw tables: `results/gate4/` (`loop.tsv`, `labels.tsv`, `tree.nwk`, `s216_reads.txt`, `reads_SRR3350473.txt`). Downloads: 4.45 GB of Gate 4 sources plus the COX1 fragment batches (logged after the fact; see the ledger). Ledger total 88.37 GB in 613 rows. Free disk 236.7 GB (floor 30 GB).

Tier A is confirmatory in scope (PREREG §2.4) except *Rotylenchulus reniformis*, which is blocked by the stem audit and so is labelled **sensitivity**. Tier B and the pre-Gate 4 checks are exploratory.

## 1. Main finding
**_Rotylenchulus reniformis_ mtRF1a has Q203 and S216 (T216S).**
- Two genome loci (`LDKF01008021.1`, `LDKF01000164.1`; assembly `GCA_001026735.1`) and TSA transcripts from five sets (GGVO, GGVP, GGVQ, GGVR, GGVV; five life stages of one project, PRJNA286314, so they count once): all 48 models give Q203 and S216. The loop reads `EGGVHRVQRLPFTEKTGRVHSSTVSVA` (0 insertions or deletions, 74–78% identical to human at these positions).
- Reads from two projects separate from the TSA sets: RNA-seq `ERR744555` (PRJEB8325) has 25 exact reads with codon TCC (serine) and 0 with any Thr codon; genome run `SRR26324223` has 6 with TCC and 0 with Thr.
- The protein groups with nematode RF1 and with *Rotylenchus robustus* (support 0.97) and is 35–38% identical to *E. coli* RF1, so it is not bacterial.
- Tylenchid Stop sisters have **T216 and Q203**: *Pratylenchus vulnus* (Gate 3 Stop, confirmatory), *Meloidogyne hapla* (Stop, exploratory UNVERIFIED set) and *Rotylenchus robustus* (Hoplolaimidae, the same family as *Rotylenchulus*). *R. robustus* has **no mitogenome and so no UAA call**. It is named as a Stop sister on weaker grounds: 6 GenBank COX1 fragments with no internal UAA (the fragment test has low power, since *Radopholus*, a known Tyr, shows UAA in only 1 of 12), plus *R. goodeyi*, the same genus, which is Stop in the exploratory set (1 internal UAA).
- *Hoplolaimus columbus* gave only a bacterial contaminant (67% identical to *E. coli* RF1, placed beside it in the tree), so it is **untested**.

**Comparison with _Radopholus_ (Gate 2): different changes in the two tylenchid Tyr lineages.** *Radopholus* has K203, T216 and a 2-residue insertion. *Rotylenchulus* has Q203, S216 and no insertion.

## 2. mtRF1a by species
| Species | UAA call | 203 | 216 | Notes |
|---|---|---|---|---|
| *Anomia* ×2, *Enigmonia*, *Placuna vitream* (Gate 2) | Tyr | Q | **S** | loop remodelled |
| *Cephalodiscus* ×2 (Gate 2 reads) | Tyr | Q | **S** | |
| *Plicatula muricata* (Gate 2) | Tyr | Q | T | loop remodelled |
| ***Plicatula plicata*** `SRR3350473` (verified: *P. plicata* RNA-seq, PRJNA316843) | Tyr lineage | Q | **T** | 3 clean reads spanning 205–220; loop identical to *P. muricata* |
| *Radopholus* (Gate 2) | Tyr | **K** | T | +2 aa |
| ***Rotylenchulus reniformis*** (sensitivity) | Tyr (blocked by stem audit) | Q | **S** | see §1 |
| *Rhabdopleura annulata* (Gate 2) | Stop | **K** | T | |
| *Heteranomia* (Gate 2) | Stop | **I** | T | |
| *Pododesmus macrochisma*, 36 background species (Gate 2) | Stop | Q | T | |
| *Pratylenchus vulnus*, *Meloidogyne hapla* (sisters) | Stop | Q | T | one residue shorter loop in both |
| *Rotylenchus robustus* (Hoplolaimidae, sister; Stop inferred from 6 fragments only) | no mitogenome call | Q | T | one residue shorter loop |
| *Albula glossodonta* (annotated `GCA_019788955.1` and hap1 `GCA_056783665.1`) | Stop, C33 | Q | **T** | Tier B |
| *Megalops cyprinoides*, *Anguilla anguilla* (sisters) | Stop, U33 | Q | T | Tier B |
| *Lycenchelys verrillii* `GCA_056747025.1` | Stop, C33 | Q | **T** | Tier B |
| *Lycodopsis pacificus*, *Zoarces viviparus* (sisters) | Stop, U33 | Q | T | Tier B |
| *Craniella infrequens*, *C. zetlandica* (hap1) | Stop, G33 | Q | **T** | Tier B |
| *Cinachyrella apion* (sister) | Stop, U33 | Q | T | Tier B |

*Microdevario* and Serrivomeridae are **untested (no nuclear assembly)**; no read search was done. *Pseudotanais* and *Placuna placenta* have been untested since Gate 2. The vertebrate loops are 89–93% identical to human across positions 196–222.

**Tier B answer.** Every confirmed non-U33 Stop species with a nuclear assembly (*Albula*, *Lycenchelys*, both *Craniella*) carries canonical T216 and Q203. **No Stop species with both a non-U33 tRNA-Tyr and S216 was found, so the combined state could not be tested.** This fits a model in which S216 matters, but it is an absence of a combination, and Microdevario and the Serrivomeridae were not tested.

## 3. Pre-registered checks (all descriptive here)
- **B3 (S216 by origin).** Confirmatory comparison is not possible: the panel genomes are excluded by §5, *Rotylenchulus* is blocked, and no origin group has ≥3 new origins with data. Descriptive: UAA=Tyr origins with loop data: S216 in Anomiidae + Placunidae, *Cephalodiscus*, and *Rotylenchulus* (sensitivity); T216 in Plicatulidae (*P. muricata*, *P. plicata*) and *Radopholus*. That is 3 of 5 with *Rotylenchulus*, 2 of 4 without. Stop species with loop data: 0 with S216 (36 background, *Heteranomia*, *Pododesmus*, *Rhabdopleura annulata*, *P. vulnus*, *M. hapla* and the Tier B sisters).
- **B4 (Stop with standalone T216S).** Not testable: no Stop species with S216 was found.
- **Necessity and sufficiency.** S216 is not necessary, because Plicatulidae and *Radopholus* read UAA as Tyr with T216. Whether it is sufficient was not tested. S216 appears in three unrelated Tyr lineages (Anomiidae + Placunidae, *Cephalodiscus*, *Rotylenchulus*), which fits convergence.
- The residue mapping of *R. reniformis* S216 has one caveat: it is at the same human position as in Zhao's T191S (human T216), by alignment.

## 4. Exploratory checks before Gate 4
### 4a. *Pratylenchus goodeyi* (`ON525593.1`)
- **GenBank has no COI or other mitochondrial record for *P. goodeyi* apart from `ON525593`** (the species has 61 records; the rest are rRNA, 18S or ITS). Its 35 internal UAA therefore cannot be tested in independent specimens, and it stays unconfirmed.
- **van de Vossenberg et al. (2023), PhytoFrontiers 10.1094/PHYTOFR-04-22-0041-FI**: the paper's DOI is the one on the record. **I could not verify the trnX/UUA claim:** the publisher's page is behind a bot challenge, which I did not try to bypass; no open copy in Europe PMC, DOAJ or the repositories I checked, and the GenBank record has only a `source` feature. The user has decided not to pursue the paper: the UUA-locus question is settled as **not evidence**, on the base-rate result below.
- **The locus.** ARAGORN reports a UUA-anticodon tRNA-like locus at 15202–15260 (59 nt, 90% AT) in `ON525593`, and tRNAscan-SE reports none. The hit on its own is not evidence, because 29 of 444 non-tylenchid nematode genomes (6.5%) give an ARAGORN call like it. So do 6 of 36 tylenchid genomes: *P. goodeyi*, *Radopholus* (390–452) and four Stop *Meloidogyne*.
- *Pratylenchus vulnus* (Stop) is the only other Pratylenchidae with a mitogenome in Gate 3. Fragments of 1,883 *Pratylenchus* records: 0 with two or more internal UAA, 8 with one at a conserved-Tyr column (scattered, unconfirmed).

### 4b. Which tylenchid families show UAA=Tyr evidence?
Fragment scan of 10,294 of 10,295 GenBank COX1 fragments in Tylenchomorpha (conserved-Tyr columns; internal UAA anywhere in the aligned frame).

| Family | Evidence |
|---|---|
| **Hoplolaimidae: *Rotylenchulus*** | **91 of 92 records have internal UAA (86 with two or more)**, across at least 10 species: *reniformis* 49 of 49, *macrosoma* 16 of 16, *macrodoratus* 8 of 8, *borealis* 8 of 9, others. The two records without are *macrosomoides* and *sacchari* (1 each) |
| Hoplolaimidae: *Helicotylenchus*, *Scutellonema* | mixed: *H. oleae* 2 of 2 (UAA at Tyr columns; same lab); *H. pseudorobustus* 10 of 11 with UAA and UAG; *S. brachyurus* 8 of 20 with UAA and 7 with UAG; *S.* sp. A 4 of 4. Most other species in these genera have none. **Candidate lineages, unconfirmed** (no read checks; could be pseudogenes or sequence problems) |
| Hoplolaimidae: *Rotylenchus*, *Hoplolaimus* | 0 (about 60 and about 220 records) |
| Pratylenchidae | *Radopholus* (known), *Pratylenchus goodeyi* (weak, one genome) |
| Heteroderidae, Meloidogynidae, Criconematidae, Aphelenchoididae, others | 0 records with two or more internal UAA among about 7,500 |

**One origin or several?** On this evidence, **at least two**: *Radopholus* and *Rotylenchulus* are in different families, each with Stop relatives, and they differ in mtRF1a (§1). *Pratylenchus goodeyi* is unconfirmed, and the other Hoplolaimidae signals are unconfirmed. The scan does not show a single origin, and it cannot rule one out.

**Corrections to the 41-of-49 statement (Gate 3B).** `transl_table` 14 is **not** evidence of UAA=Tyr on its own: *Ditylenchus dipsaci* and *D. gigas* fragments and *Radopholus* fragments with no UAA are also filed under table 14, so the choice reflects labs' practice. Within *Rotylenchulus* all but the missing annotations use it. The fragments themselves, with internal UAA in 91 of 92 records, are the evidence.

**Power.** Fragments are short. Two or fewer Tyr columns are typical, and even *Radopholus* (known Tyr) shows UAA in only 1 of 12 records. Absence is weak evidence for any family with few records.

### 4c. *Rotylenchulus* internal TAG sites (`CM003310.1`)
`results/gate3/roty_tag_sites.txt`. Relatives are the 10 closest Stop species (none closer than order).
| Site | p(Tyr) | Relatives' residues | Reads |
|---|---|---|---|
| COX1 column 404 | 0.02 | M 4, L 4, S 1, N 1 | **RNA-seq: 5 TAG**; genome: 434 TAG, 5 TAT |
| COX3 column 239 | 1.00 | Y 8 of 8 readable | genome: 227 TAG, 2 TAA, 1 TAT; no RNA-seq coverage |
| CYTB column 62 | 0.04 | S 9, D 1 | genome: 353 TAG, 1 TAT; no RNA-seq coverage |

**Only COX1 404 appears in RNA-seq.** Only COX3 239 is a tyrosine site. So UAG is read as sense somewhere (or these are not stops), but that it is tyrosine is not shown at two of the three sites.

## 5. Caveats
- The *Rotylenchulus* genome is fragmented (N50 6 kb); the S216 result does not depend on assembly because the TSA and reads agree.
- Assemblies and proteins: *H. columbus* contamination was found and dropped; *R. robustus* has no mitogenome call.
- *Plicatula plicata*: 3 independent clean reads, at the pre-set minimum; the fourth read is truncated and has an odd end.
- In the *Rotylenchulus* result I chose loci after seeing the alignment, and the loop alignment uses the Gate 2 method.
