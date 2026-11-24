# Gate 3B: exploratory UNVERIFIED run and candidate audits

Date: 2026-11-24. **Post-hoc and exploratory. Nothing here changes the committed Gate 3 totals** (`7e0541f`). Gate 4 is on hold. The 5-minute hit watcher ran during the UNVERIFIED run at the user's request, so this set was not analysed blind. The classification and audits use the locked functions unchanged (`src/gate3u_classify.py`, `src/gate3u_weak.py`). My exploratory weak-audit scripts first omitted the Gate 3 Tyr genomes from the replication test. The fix changed only *Pratylenchus goodeyi*. Outputs: `results/gate3_unverified/` and `results/gate3/` (audits).

## 1. UNVERIFIED-flagged species
Of 8,189 excluded species, 7,707 ran (480 excluded for >1% ambiguous bases, 2 excluded by PREREG). Run complete: 7,707 of 7,707, no module errors, 212 species also in Gate 3. `out.jsonl` sha256 prefix `d68cdb68efc65b5d`.

| Result | n |
|---|---|
| Stop | 7,636 |
| ambiguous | 41 |
| insufficient | 18 |
| weak | 12 |
| **Tyr** | **0** |
| trnY 33 resolved / unresolved / stem_audit_review | 6,560 / 1,002 / 145 |
| non-U33 among resolved Stop | 1 of 6,506 = 0.015% (Wilson 0.003–0.087%) |

Coverage by phylum (n, resolved): Arthropoda 5,323 / 4,704; Chordata 1,008 / 870; Bryozoa 556 / 322; Mollusca 349 / 307; Annelida 131 / 123; Platyhelminthes 87 / 82; Porifera 67 / 60; Cnidaria 54 / 0; Nematoda 45 / 33; Bivalvia 150 (all Stop, 137 resolved).

**Weak calls (12). The locked audit rejects 11; 1 is "supported weak".**

| Species | Group | Internal UAA / genes / Tyr sites | Audit |
|---|---|---|---|
| *Pratylenchus goodeyi* `ON525593` | nematode | 35 / 7 / 9 of 28 | **supported weak** (10 sites, 4 genes; replication = *Radopholus similis*, panel) |
| *Ortalia* sp. `PX770128` | ladybird | 2 / 2 / 1 | rejected: 1 supported site |
| *Scolopendra morsitans* `MW810062` | centipede | 1 / 1 / 1 | rejected: 1 supported site |
| *Nasutitermes* sp. `OQ675028` | termite | 1 / 1 / 1 | rejected: frameshift, internal UAG |
| *Joubiniteuthis portieri* `MW255583` | squid | 3 / 1 / 1 | rejected: frameshift |
| *Chrysso nordica* `MW776437` | spider | 1 / 1 / 1 | rejected: 1 supported site |
| *Spilarctia alba* `KX753670` | moth | 9 / 2 / 1 | rejected: frameshift in 5 genes, LLR −29.8 |
| *Metynnis luna* `PV962571` | fish | 3 / 2 / 1 | rejected |
| *Lophiobrycon weitzmani* `PV962563` | fish | 4 / 3 / 1 | rejected: 7 genes in >1 frame |
| *Ammocryptocharax elegans* `PV962435` | fish | 3 / 3 / 1 | rejected: 5 genes in >1 frame |
| *Anser* sp. `KY626009` | goose | 1 / 1 / 1 | rejected |
| *Sotalia guianensis* `KM893424` | dolphin | 7 / 5 / 2 | rejected: 0 supported sites |

Read tests (record codon at each UAA versus TAT/TAC/TAG; exact 31-mers):
- *Metynnis luna*, *Ammocryptocharax elegans*, *Chrysso nordica*: the reads carry a tyrosine codon where the record has UAA at the conserved-Tyr site: TAT in 27 *Metynnis* reads, TAC in 5 *Ammocryptocharax* reads, and TAT 379 + TAC 92 against TAA 1 in *Chrysso*. **So the UAA in these three records is an error.**
- *Lophiobrycon*, *Sotalia*, *Joubiniteuthis*, *Scolopendra*: no reads covering the site or no matching mitochondrial reads in the run, so uninformative. Others have no usable public reads.
- **Suspected bad submission: the `PV962` fish batch.** 216 records with this accession prefix (the two I opened are linear Megahit assemblies that GenBank could not verify) give 20 non-Stop calls (9.3%: 17 ambiguous, 3 weak), against 0 of 304 other bony fish and 51 of 7,187 (0.7%) elsewhere. I treat their calls as one suspected assembly problem and list them separately. Three were tested against reads: two records are contradicted by the reads and the third was not covered. I did not test the other 213.

***Pratylenchus goodeyi*** is the only call of interest. UAA falls at conserved-Tyr columns 19 times more often than elsewhere (9 of 28 against 26 of 1,072; post-hoc Fisher p = 5×10⁻⁸). Against that, LLR is −2.68, 18 of 35 UAAs sit where Tyr probability is below 0.1 and 5 best-fit phenylalanine; no conspecific record, reads or assembly exist in NCBI. GenBank could not verify the sequence. *Pratylenchus vulnus* (Gate 3) is Stop with 0 internal UAA. Supported weak is not Tyr, and appears only in the sensitivity analysis.

**Non-U33 Stop: one, *Desmacella informis*** `OM729615` (sponge, A33, score 10.1, 4/5 stem, 1 caller). It passes the locked rules. I judge it a probable artefact: linear 12.9 kb record with 6 of 13 genes, the 23-base stem sequence is a direct repeat, and ARAGORN calls two overlapping Tyr tRNAs 30 bases apart. No conspecific record, reads or assembly to test it. No second independent A33 exists.

## 2. Candidate audits (non-U33 Stop candidates from Gate 3)

| Candidate | Independent support | Verdict |
|---|---|---|
| *Microdevario* (2 sp.) | earlier checks; see Gate 3 report | holds |
| Albuliformes (4) | independent *A. vulpes* mitogenome; reads 469 to 0; **0 NUMT over N33 in 2 of 3 *A. glossodonta* assemblies** (hap1, BYU) | holds |
| *Serrivomer* / *Stemonidium* (4 sp.) | *S. sector*: 4 independent records, reads 363 to 0. *S. beanii*: 3 others, reads 272 to 0. *Stemonidium*: 4 other specimens, reads 85 to 0. *S. lanceolatoides*: 1 record, reads 419 to 0 from its own specimen | family holds |
| *Lycenchelys verrillii* | 4 records C33, 2 from USNM vouchers; reads 622 to 1 (USNM specimen, not the record's) | holds |
| *Coloconger meadi*, *Ophichthus asakusae* | own-specimen reads only (88 and 200 to 0) | unconfirmed beyond the specimen |
| *Xyrias revulsus*, *Toxocara apodemi* | 1 record (RefSeq copy of the same for *Toxocara*), no reads | unconfirmed |
| Crane (*Antigone rubicunda*) | 1 record and duplicate, stem 3/5, no reads; 33 other *Grus* records U33 | unconfirmed |
| *Craniella* (2 sp., G33) | reads: *C. infrequens* DNA 40 to 0, RNA from a second study 66 to 0 and own project 2 to 0; *C. zetlandica* DNA 55 to 0, RNA 35 to 0. Both records come from the same Sanger pipeline (BioProjects PRJEB107262, PRJEB96775; different specimens) | holds as a *Craniella* feature |
| *Globodera ellingtonae* | reads C33 in two populations (176 and 271 to 0), but chromosome II carries another GUA tRNA-like locus with **U33** (score 11.1; the C33 locus scores 10.0) | **dropped**: the C33 locus may not be tRNA-Tyr |
| *Bufotes pewzowi* | RAD-seq runs hold no mitochondrial reads; no independent conspecific | **dropped** |

Serrivomeridae carries C33 in all four species tested. Two GenBank labels look wrong, and I set both records aside. `OR546184`, labelled *S. sector*, has U33 but matches *Antennarius striatus* (a frogfish) at 90.7% and no eel above 72%. `OR482581`, labelled *Stemonidium*, is *S. sector* (98.5–99.2% identical to four *S. sector* records). Neither was a Gate 3 genome, so the Gate 3 counts do not change.

**Eel origins of C33: 3 or 4.** Counting by NCBI taxonomy gives 4 changes within Anguilliformes: Serrivomeridae, *Ophichthus*, *Xyrias* and *Coloconger*. The gene tree (FastTree, 9 mitochondrial genes, 242 tips, rooted on Elopiformes) joins *Xyrias* with *Ophichthus asakusae* (support 1.0), which gives 3. Support for the deeper nodes is weak, and Anguilliformes itself has support 0.54. *Coloconger cadenati*, the other *Coloconger* in the tree, is U33. Albuliformes does not form a single group on this tree, so I did not use it to count origins in that lineage.

**Craniella details.** Whole-mitogenome identity: the two species 97.1%; to *C. villosa* / *C. arb* (U33) 86.7–88.8%; to *Cinachyrella* (U33) 88.9–90.8%. The trnY is identical in both species (`CGGUAAA`, 4/5 stem); Cinachyrella and *C. arb* / *villosa* have `CUGUAAA`, 5/5. One U→G change plus a stem-base change separates them. One origin is the simplest reading. A second GUA-like hit turns up in these Tetillidae genomes, but in the annotated `NC_010198` it lies inside the ND2 coding sequence, so there is **no second trnY**. Reads from *C. arb* itself show 216 with U33 and none with another base. The second *C. infrequens* run's species label is from its SRA record and was not checked separately.

## 3. NUMT checks
| Species (assembly) | Nuclear hits ≥200 bp, ≥80% | Covering N33 |
|---|---|---|
| *Albula glossodonta* hap1 `GCA_056783665.1` | 0 (assembly holds the mitochondrion as contig `CM165726.1`) | 0 |
| *A. glossodonta* BYU `GCA_019788955.1` | 0 | 0 |
| *Lycenchelys verrillii* `GCA_056747025.1` | 0. An 82 kb contig (`JBQQXK010000041.1`) is ≥99% mitochondrial over most of its length (tandem copies or assembly artefact); excluded | 0 |
| *Craniella infrequens* hap1 | 217 | 4, all G (record + strand base G) |
| *C. zetlandica* hap1 | 322 | 38, all G |
| *Globodera ellingtonae* | 10 | 0 |

The method found a planted 10%-diverged 2 kb copy. The accession you gave, `GCA_056783515.1`, is *Draco maculatus*; the bonefish isolate is **`GCA_056783665.1`**. For *Craniella*, nuclear copies also carry G33, so the assemblies cannot separate a mitochondrial from a nuclear origin of the record. The read counts (no U33 in any read) show no U33 copy at comparable depth. Other haplotype assemblies (hap2) and eels, *Microdevario*, *Xyrias*, *Toxocara* were not searched: no nuclear assembly.

## 4. *Rotylenchulus reniformis* (Gate 3 loose end)
**Strong candidate, blocked by the pre-registered stem audit** (trnY stem 3/5, score 12.8, U33; `CM003310.1`). Evidence: 32 internal UAA in 6 genes, 15 at Tyr columns, LLR 67.4, Codetta Y.
- Independent genome reads (`SRR26324223`, PRJNA1026095): all 32 TAA sites present at depth 70–474, 10,613 reads TAA against 23 TAT/TAC in total.
- RNA-seq (`ERR744555`, PRJEB8325): 16 of 32 TAA sites covered, 55 reads, 0 Tyr-codon reads. Low depth; COX2 and ND1 sites have none.
- **Internal TAG sites (loose end b):** WGS reads: COX1 column 404 TAG 434 (TAT 5), COX3 239 TAG 227 (TAA 2, TAT 1), CYTB 62 TAG 353 (TAT 1). RNA-seq: COX1 404 TAG 5; the other two have no coverage. Model probability of tyrosine: 0.02, 1.00, 0.04. The TAGs are real in the DNA; the COX1 one is also in RNA. COX3 239 is the only one in a conserved-Tyr column, so UAG may be read as Tyr there, or the column may be a site variation; unresolved.
- 49 independent *COX1* fragments: TAA in 47–49 of 49 at three sites, and 41 of 49 (TAG in 8) at a fourth. **41 of 49 fragments are already filed under translation table 14, in which UAA is tyrosine.** No paper stating it was found.
- It does not meet PREREG §3 as a new clade while the locus audit blocks it. Relatives tested (*Rotylenchus goodeyi*, *Hoplolaimus galeatus*, *H. columbus*, *Meloidogyne* ×3, *Heterodera* ×3, *Nacobbus*): UAA = Stop or insufficient, except weak *Pratylenchus goodeyi*.

## 5. Answer: does PREREG require the 8.5 Gb Bufo/Rhinella download?
**No.** PREREG §2.4 lists mtRF1a taxa as reassigned lineages, matched sisters, the background panel, any new Tyr clade with at least two sisters, and *Plicatula plicata* (`SRR3350473`). Non-U33 Stop species, including *Bufotes* (which is also dropped here), are not on the list. My Gate 3 estimate (20 assemblies, 19.7 Gb) used the wider wording of the original brief ("every species flagged"), and the 8.5 Gb was genome length (*Bufo bufo* 5.04 + *Rhinella marina* 3.47), about 2.4 GB compressed. *Bufo bufo* (`GCF_905171765.1`) has annotated proteins, so a protein query would be the first source under §2.4. `GCA_900303285.2` (*Rhinella*) was not checked for annotation.

## 6. Corrections to earlier statements
- **Albula nuclear assemblies.** The committed Gate 3 report (line 67) says no species among the three named candidates has its own nuclear assembly. That is wrong for *Albula*: *A. glossodonta*, the species of record `CM165726.1`, has three (`GCA_056783665.1` hap1, `GCA_056776265.1` hap2, `GCA_019788955.1`). *A. vulpes* has none. The NUMT search was run on hap1 and `GCA_019788955.1` (§3), not on hap2. The Gate 3 report is left as committed; this is the correction.
- The Gate 4 size in the Gate 3 report follows the original brief, which is wider than PREREG §2.4 (see §5).
- *Globodera* and *Bufotes* no longer count as candidate evidence. *Coloconger*, *Ophichthus asakusae*, *Xyrias*, *Toxocara* and the crane stay unconfirmed, because no independent specimen or reads of another individual exist for them.

## 7. Downloads and resources
Ledger: 542 rows, **79.54 GB cumulative** (Gate 3 was 28.32 GB; most of the rest is reads that were deleted after use). Disk free 243.5 GB (floor 30 GB).
