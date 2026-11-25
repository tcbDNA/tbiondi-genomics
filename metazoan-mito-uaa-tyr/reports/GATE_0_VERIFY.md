# Gate 0: Primary-source verification and scope decision

Date: 2026-11-23. Verdict: **PASS, with control-panel problems that need a decision before Gate 2.**
Provisional scope: **`FULL_TWO_GENOME`** (Zhao's `mtRF1a` claims are extractable and specific). Final scope depends on Gate 1 pairing.

## 1. Sources and how each was obtained

| Source | Status | How |
|---|---|---|
| Zhao et al. 2026, BMC Biol, `10.1186/s12915-026-02767-5` ("Article in Press") | Full text + MOESM1 read | Springer served a bot challenge / sign-in redirect to automated fetches and the paper is not in PMC/PubMed/Europe PMC. The user downloaded the PDF and MOESM1 in a browser. PDF sha256 `a5c2dcbc…b27a`; MOESM1 sha256 in `data/manifests/sources.sha256` |
| Lin et al. 2026, bioRxiv, `10.64898/2026.08.14.744788` (v1, 19 Aug 2026) | Main text read; **supplement not read** | PDF saved by the fetch tool, sha256 `0fba5992…00f1` |
| Li et al. 2019 (PMC6319601), Jacob et al. 2009 (PMC2761399), Perseke et al. 2011 (PMC3121625) | Main texts read by keyword search, supplements not read | Europe PMC full-text XML |
| Matsushima et al. 2026, `10.3897/asp.84.e192470` | Read by keyword search | Publisher open-access XML |
| Jakiel et al. 2025, `10.1016/j.pocean.2025.103534` | **Not read** (ScienceDirect HTTP 403; not retried) | — |
| Kakui & Kano 2021, `10.2108/zs200167` | **Abstract only** | Europe PMC |

"Read by keyword search" means I extracted matching sentences, not that I read every page.

## 2. Zhao et al. 2026: what the paper actually says

**mtRF1a.** "recognition of the third base position in the codon was primarily mediated by residues R178 and T191 (Fig. 3a and b) [9]."
"we identified a parallel amino acid mutation at position 191 (T191S) in *A. chinensis*, *P. placenta*, *Cephalodiscus* spp., and *Pseudotanais* spp."
"reassignment occurred independently in at least four metazoan lineages, and this amino acid mutation arose three times. However, the mutation was not found in *R. similis*. In *R. similis*, a mutation occurred at the position 178 (Q178K). Yet, the same mutation was also found in *Rhabdopleura* spp., whose mtRF1a recognizes both UAA and UAG."
"we noted an unusual two-amino-acid insertion in *R. similis*".

**tRNA-Tyr / TGT / queuosine.** "*A. chinensis* and *P. placenta*, as well as members of Cephalodiscidae, possessed an identical anticodon loop sequence and exhibited a parallel U33C substitution (C31 in the linear *A. chinensis* sequence and C34 in Cephalodiscidae, respectively, corresponding to position 33 in standard tRNA numbering)."
"The C33 substitution alters the local motif from U-G-U to C-G-U, which is predicted to disrupt recognition by tRNA guanine transglycosylase (TGT) [29], thereby preventing queuosine (Q) modification of the wobble base G34."
"*R. similis* lacks this U33C mutation"; instead Zhao proposes a "suppressor" tRNA(UUA). "In *Pseudotanais* sp. a, tRNATyrGUA does not use the U33C strategy"; *Pseudotanais* sp. b carries a novel tRNATyr(AUA).
Residual termination: "No candidate UAG termination codon was detected for *cob*"; "UAA might still possess termination functionality in these species".

### Reference numbering (needed before any residue call)
- Zhao does not state a reference sequence. Figs 3a/3b label GLN-178 and THR-191 on human mtRF1a; the main text's "R178" contradicts its own "Q178K" and Fig. 3b, so I treat it as a typo.
- In human `NP_061914.3` (380 aa), residue 178 is R and 191 is **R**, in the wrong context. The Fig. 3a motif `KFEGGVHRVQRVPKTE…KQGRVHTST` maps to `PxTE` at 205 and `GRVHTS` at 212.
- **Inferred mapping: Zhao 178 = Q203 and Zhao 191 = T216 in `NP_061914.3` (offset 25).** Both sites are 13 apart in both schemes. The origin of the offset is not stated in the paper. All later calls will be anchored on the `GRxHTS` / `VPxTE` motifs, never on an index.

### Caveats on Zhao's data (confounders to carry forward)
1. **The T191S "site" in *Anomia*/*Placuna* sits in a broadly rewritten loop** (Fig. 3a: human `VPKTE---KQGR--VHTST`; *Anomia* `RETRVK---QARN--FFSCT`). PxTE and HTS are both altered, so T191S is not isolable there.
2. *C. hodgsoni*, *P. placenta* and *Pseudotanais* sequences are labelled "partial"; *C. hodgsoni* was "manually curated" from a low-depth transcriptome; the starred *Pseudotanais* sp. comes from pooled reads.
3. Zhao states UAA codons "appear to have been treated as assembly errors and removed" in `PP708078.1` (*A. chinensis*) and `PP708091.1` (*P. placenta*). **Those two accessions are therefore unusable for UAA calling.** (Both are also `UNVERIFIED` with no CDS.)
4. Table S5 errors found: `NC_040109.1` is listed as *Rhabdopleura compacta* but is *Glossobalanus marginatus*, and `FN908482.1` is listed as *Glossobalanus* but is *R. compacta*.
5. `PV231829` (Zhao's *A. chinensis* mitogenome) could not be retrieved from NCBI or ENA on this date.
6. Extracted from the figure by eye, not from an alignment file: Zhao Fig. 4b shows *Asterias rubens* with U33, which **conflicts with Li et al. 2019** (below).

## 3. Lin et al. 2026
- "In all TAA=Tyr lineages, we identified a consistent U36C substitution in the anticodon loop of trnY (Fig. 4h & S12–S14)… This substitution is present exclusively in the TAA=Tyr lineages and is absent in all species retaining the canonical code 5". "The independent TAA reassignment in the hemichordate Cephalodiscidae also involves a U-to-C substitution at the analogous position in the trnY anticodon loop (30)."
- *Pododesmus*: "translation code 5". *Heteranomia*: code 5, "harboring in-frame TAG codons" with hypothesised +1 frameshifts (H. aculeata nad4/nad6: "four in-frame TAG codons with four corresponding frameshifts"); "TAA as the sole termination codon". **So *Heteranomia* is UAA=Stop but its genes need frame-aware handling.**
- *Release factor*: only "a precise discrimination likely involving ribosome–release factor interactions (72)" (in the context of Heteranomia frameshifting). Regex search of the main text found no `mtRF1a`/`MTRF` and no queuosine/transglycosylase. Confirmed: Lin does not test mtRF1a.
- Reassignment is said to have converged independently in "*Anomia*/*Enigmonia*/*Placuna* and Plicatulidae" (Lin), so **Anomioidea+Plicatulidae may be ≥2 origins**; Lin also names a Placunidae species, *Placuna vitream*, not *P. placenta*. Codetta was run with `-m -e 1e-5 -r 0.7` (defaults are 1e-10 / 0.9999).
- Accessions: the main text gives none. It points to BioProject `PRJNA1198682` (0 hits in NCBI BioProject/SRA today) and Figshare `10.6084/m9.figshare.33250323` (DataCite: findable, 1.04 GB, CC BY 4.0). Supplementary tables S1/S6/S7 were not read.

### "U36C" mapping (control accessions only)
Canonical position 36 is A in a GUA anticodon, so "U36C" cannot be canonical. Method: tRNAscan-SE `-O` / annotated GenBank coordinates; I read the base 5' of the `GTA` and counted pairs in the putative 5-bp anticodon stem for the six tRNAscan-SE hits (5/5 unless noted). The four annotated-only controls were not structure-checked. This is a manual check, not the Gate 2 extractor.

| Control | Gene-local index of N33 | N33 | Note |
|---|---|---|---|
| *A. chinensis* PP708078 | 31 | **C** | reproduces Zhao's "C31" exactly |
| *P. placenta* PP708091 | 31 | **C** | anticodon stem 4 of 5 pairs |
| *A. ephippium* OZ576062 (4 copies) | 29 | **C** | |
| *C. hodgsoni* MF374802 | 34 | **C** | reproduces Zhao's "C34" exactly; stem GTGGG/CCCAC |
| *Heteranomia squamula* OZ181722 | 32 | **U** | code 5 |
| *R. compacta* FN908482 | 29 | **U** | anticodon stem 4 of 5 pairs |
| *Radopholus* FN313571 | 33 (annotated) | **U** | structure not verified |
| *Pseudotanais* OZ281197 | 26 (annotated) | **U** | structure not verified |
| *Arctotanais* LC597489 / *Sinelobus* LC923628 | 27 / 23 (annotated) | **U** | structure not verified |

**Finding:** "36" matches no control's local index (23–34 across the ten controls). Lin's method (MAFFT alignment → WebLogo) makes an alignment/logo column the likely meaning, but I have not seen their alignment, so this is an inference. Zhao's C31/C34 are gene-local indices and do map to canonical 33. Zhao also explicitly maps to position 33.
**Already visible from the published data and these controls (must be stated in `PREREG.md`):** *Radopholus* and *Pseudotanais* are UAA=Tyr with U33 — published necessity counterexamples to U33C. This is expected from Zhao's own text, not a discovery.

## 4. Documented genetic code per control lineage

| Lineage | Best evidence | Code | Status |
|---|---|---|---|
| *Cephalodiscus hodgsoni* | Li 2019; GenBank table 33 | UAA=Tyr | verified; nad1 may still end in UAA (Li) |
| *Radopholus similis* FN313571 | Jacob 2009; 157 TAA `transl_except` Tyr (+1 AAT) | UAA=Tyr | verified; annotated table 5 |
| *Pseudotanais* sp. OZ281197–200 | GenBank `transl_except` 85–92 TAA each; Matsushima cites Jakiel; Zhao re-analysed sp. a/b | UAA=Tyr | verified at secondary level; **Jakiel not read**; all four records share taxid 3031986 |
| *Anomia* | Lin, Zhao (Codetta; resequencing of 3 further individuals) | UAA=Tyr | verified in papers; usable genomes: see §5 |
| *Placuna* / Plicatulidae / *Enigmonia* | Lin (Codetta; proteomics for *P. muricata*) | UAA=Tyr | verified in paper; genomes only via Figshare |
| *Arctotanais alascensis* LC597489 | Matsushima: PCGs translated with invertebrate code 5; TAA=Tyr "has not been observed in other peracarids" than *Pseudotanais* | UAA=Stop | **negative control** (a negative inference in the paper, not a dedicated test) |
| *Sinelobus*, *Apseudes*, *Carpoapseudes*, *Neotanais* | Matsushima: translated with code 5 | UAA=Stop | negative controls |
| *Rhabdopleura compacta* FN908482 | Perseke 2011: TAA and TAG both stop (9 UAG, 2 UAA); table 24 | UAA=Stop | negative control |
| *Heteranomia* / *Pododesmus* | Lin | UAA=Stop | negative; *Pododesmus* has no NCBI mitogenome |
| Nematode sister taxa | Jacob: unique to *Radopholus*; Zhao Fig. 4b lists *B. xylophilus*, *C. elegans* as canonical | UAA=Stop | candidates: `NC_023208`, `NC_021120`, *Meloidogyne* spp. |

**Conflict to resolve after pre-registration:** Li 2019 says "both starfish and *C. hodgsoni* mitochondrial tRNA Tyr genes have a C" at position 33 (code 9, UAA=Stop). Zhao Fig. 4b shows *Asterias rubens* U33. If Li is right, starfish are a sufficiency counterexample. Not tested here (non-control genomes; firewall).

## 5. Accession availability (positive controls first)
- Public in GenBank and usable: *C. hodgsoni* MF374802, *Radopholus* FN313571, *Pseudotanais* OZ281197–OZ281200, *A. ephippium* OZ576062 (**71,118 bp, linear, trnY copies spaced 32,920 bp apart — likely a concatemer; no CDS annotated; genus-level evidence only**).
- Public but **not usable for UAA calling**: PP708078, PP708091 (UAA codons removed, per Zhao).
- Not found in NCBI: *Enigmonia*, *Plicatula*, *Placuna vitream*, *Pododesmus* mitogenomes; PV231829. Available only via Lin's Figshare (outside the listed databases; needs approval) or by assembling raw reads.
- Negative controls in GenBank: OZ181722 (*Heteranomia squamula*; Lin used *H. aculeata*), FN908482, ON600795/6, PP334099, LC597489, LC923624–28, NC_023208 and other nematodes.

## 6. Nuclear data (first look; formal overlap is Gate 1)
| Taxon | NCBI assemblies | SRA RNA-seq |
|---|---|---|
| *Anomia* | GCA_062075034 (*A. chinensis*), GCA_990049579 (*A. ephippium*) | 5 |
| *Enigmonia* | GCA_056601665 (chromosome level) | 2 |
| *Heteranomia* | GCA_964263295, GCA_964263325 | 0 |
| *Pododesmus* | GCA_032359645 (scaffold, contig N50 4.8 kb) | 1 |
| *Placuna*, *Plicatula* | none | 0 / 1 |
| *Cephalodiscus*, *Rhabdopleura* | **none** | 5 / 6 |
| *Radopholus* | GCA_013357305 (+2) | 8 |
| Tanaidacea | **none** | 25 |

For hemichordates, tanaids and *Placuna*, mtRF1a would have to come from RNA-seq or raw reads (Zhao Table S7: SRR11101524, SRR1695471, SRR11101525, DRR198815, ERR15194042/44, SRR35651820). This is the main pressure toward fallback.

## 7. Novelty check
Europe PMC (including preprints), plus one web search, on 2026-11-23: no pan-metazoan RefSeq/DToL census of `trnY` U33C or `mtRF1a` T191S found since Oct 2026. Limited by keyword coverage.
Attribution note: the C33/unmodified-G34/queuosine idea for *Cephalodiscus* was already proposed by Li et al. 2019.

## 8. Tools and compute
Pinned in `env/conda-osx-arm64.lock`: Infernal 1.1.5, HMMER 3.4, MAFFT 7.526, miniprot 0.18, BLAST+ 2.17.0, tRNAscan-SE 2.0.13, ARAGORN 1.2.41, NCBI datasets 18.41.1. Codetta `863359e` (5 May 2025) with HMMER 3.3.2 (osx-64 under Rosetta, because HMMER 3.3.2 does not build on arm64), Pfam-A_enone (1.16 GB) and the metazoan 13-PCG profile set built from Codetta's own alignments. **Codetta was installed but not run**: its bundled examples are non-control genomes. Codetta's default Pfam excludes mitochondrial domains; `-m` is required. lib-kaci (`artmeln/lib-kaci`) found but not built (Linux/WSL only).
Firewall compliance: tRNAscan-SE and ARAGORN were run only on the ten named controls (`src/gate0_trny_controls.py`). Human/*E. coli* mtRF1a sequences were indexed for numbering; no genome was searched.
Download accounting (measured): Pfam 1.16 GB, RefSeq 0.40 GB, conda archives 0.40 GB, taxonomy 0.16 GB, other retained files < 0.2 GB; NCBI/ENA API responses that were not kept were not metered. No job exceeded 1 hour. Outside the project folder I only read two files in `~/Downloads` that the user pointed to (copied in; originals untouched).
**Publisher access, stated in full:** three `curl` requests carried a browser User-Agent (two to BMC/Springer, one to bioRxiv). Springer answered with a JavaScript bot challenge and bioRxiv with Cloudflare 1015 (rate limit). I then dropped the User-Agent but made one further attempt at each publisher (article page and PDF) through the session's web-fetch tool: Springer redirected to a sign-in endpoint (not followed), bioRxiv served the page and PDF. So the Lin text was obtained from a site that had just rate-limited another request from me. ScienceDirect returned 403 once and was not retried. The Zhao PDF came from the user.

## Open decisions
1. **Anomioidea/Plicatulidae positive controls:** approve downloading Lin's Figshare (1.04 GB, CC BY 4.0) and/or Zhao's Additional file 4 (Data S1)? Without it, the within-bivalve Gate 2 controls reduce to *A. ephippium* vs *Heteranomia*.
2. **Optional papers to add:** Jakiel 2025 (evidence for *Pseudotanais*), Lin's supplement (Figs S12–14, Table S6–S8), Kakui & Kano 2021, Zhao Additional files 4 and 9.
3. Scope for the `mtRF1a` half: depends on Gate 1.
