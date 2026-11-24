# Cold review of REPORT.md (at commit 7ba7789)

**Status of this document.** Received 2026-11-24. The text between the two rules below is the reviewer's, reproduced verbatim as pasted by the user. Nothing in REPORT.md, the figures or the dossiers has been edited in response yet.

**Independence note (from the user).** The review session loaded the project's CLAUDE.md and memory. It is therefore a fresh-session, same-family adversarial review, and not a strictly isolated one. The original brief asked for a *different* model in an isolated session with no prior context; that condition is not met, and this is a deviation to disclose.

**Verification log and resolutions.** To be added after the read-only raw-data verification that the user requested. No critique point has been accepted or rebutted in this file.

---

Cold review

  Scope: No files were attached, so I read PREREG.md, reports/*.md, the dossiers and the figures in the repo at 7ba7789. I didn't open src/ or results/. I used git history for one check (A1).

  Independence: This session loaded the project's CLAUDE.md and memory index, and the attribution settings suggest earlier Claude sessions did the work. A truly independent review needs a clean session.

  Short answers

  - Q1, most likely wrong: "At least two origins … in tylenchid nematodes" (REPORT §5e; see A2). The C33 result probably survives, but with narrower reach than stated (A3, A4).
  - Q6, Rotylenchulus: Not a new clade, and not "sensitivity" either. PREREG §3 already has the label: "unconfirmed single-genome candidate." Reasons are at the end.

  A. Critical

  1. §3.6 misstates when the stem-count problem was found (Q4, Q6).
     - §3.6 says it was found "while building the dossiers" and "after the sweep".
     - GATE_2 A13, "Final pre-lock checks", already records that a contiguous-arm recount disagrees in 4 of 305 calls. Git shows that text entered in the lock commit 5aa5774.
     - The ambiguity was measured, then locked unspecified. REPORT never mentions A13.
  2. §5e upgrades its source (Q1, Q3, Q7).
     - GATE_4 §4b says the scan "cannot rule out" a single origin. REPORT drops this and says two are "implied".
     - The evidence is a blocked genome, NCBI family labels, fragments of admitted low power, and different mtRF1a changes, which could postdate a shared reassignment.
     - The dossiers show 9 shared UAA columns between Rotylenchulus and P. goodeyi, one of them at a Phe site (COX1 324). This is not discussed.
     - From background knowledge: molecular trees often place Radopholus outside Pratylenchidae. Check this.
  3. The C33 result misses the prior claim's domain (Q3, Q5).
     - GATE_0 §3 quotes Lin as saying the change is absent in species "retaining the canonical code 5". REPORT §1 omits this.
     - All 9 supported C33 genomes are vertebrates. Among code-5 taxa the report has only Toxocara (unconfirmed) and Globodera (dropped).
     - Genetic code is never treated as a confounder. GATE_3's concession of 0 of 1,491 molluscs is also dropped.
  4. "9 supported" has no stated criterion and is applied inconsistently (Q3).
     - S. lanceolatoides (own reads only), M. nanus and Pterothrissus (no species-level checks) are "supported".
     - Coloconger and O. asakusae (own reads) are "unconfirmed". So is Xyrias, which the gene tree joins to O. asakusae at support 1.0.
     - Bufotes is "dropped" on roughly the evidence that leaves Xyrias and Toxocara "unconfirmed".
     - Requiring independent same-species evidence leaves about 5 genomes in 3 lineages. The sufficiency claim survives, but the abstract's number does not.
  5. A pre-registered T6 check was not run, and panel evidence stands in for it (Q4).
     - T6 requires read support if the source BioProject has reads. The A. glossodonta dossier says no read test was done on the record, and no document says whether reads exist.
     - "469 to 0" in §5a is from A. vulpes, a panel genome, and the table does not say so. Neither counted Albuliformes genome has read support. §3.5 omits this.
  6. Second tRNA-Tyr candidates are unexplained (Q5, Q7).
     - Dossiers report 2 verified GUA arms in A. glossodonta, Pterothrissus, both Microdevario, Lycenchelys, Coloconger, Bufotes and both Craniella. ARAGORN finds 2 Tyr loci in S. beanii and Xyrias. Only Craniella
       is explained.
     - GATE_3 T6 says "no second trnY locus … in all four", which appears to conflict.
     - Globodera was dropped for a second locus.
  7. The new-clade endpoint was never calibrated (Q2, Q4).
     - PREREG §2.3 requires the trnY audit "before any genome counts as a … new clade". §3.6 quotes that bullet without this clause.
     - Of 20 known positives, 4 have no confident trnY call and both Enigmonia fail the audit (GATE_2 §3, A12). At most 14 of 20 pass; Radopholus and Enigmonia score 0.
     - "20 of 20 recall" (§2) describes the UAA call only. The [C] null cannot see Cnidaria (0 of 892) at all.

  B. Major

  Circular or uninformative tests (Q2)

  8. False-positive bound: "0 of 340" and its 1.1% bound use the negatives that set the thresholds, so zero is guaranteed. A11's leave-one-phylum-out test is the fairer figure, and REPORT omits it.
  9. Fig 3, right panel: "Negatives-only" reproduces the fixed rule exactly (A11) and uses no positives. Holding out a lineage cannot change it, so two of four columns are the same result.
  10. T2 and A1: T2's 0 of 17 is forced, because no resolved Tyr genome exists. A1's 1 of 1 is near-certain at a U33 rate above 99%, and its only origin is the blocked genome.
  11. Read checks:
      - Most use runs from specimens whose assemblies already show C33.
      - Exact 31-mers miss any U33 copy with another nearby difference. S. beanii matched 0 reads for that reason.
      - Globodera passes 447 to 0 and was still dropped.
  12. GATE_3B Fisher p = 5×10⁻⁸: P. goodeyi was selected partly for UAA at Tyr columns, then tested for it. PREREG §5b bars p-values as tests in exploratory work.

  Labels and overclaims (Q3)

  13. §6 "[C]": Most rows are panel data disclosed in PREREG §0. REPORT defines [C] as non-panel only.
  14. Abstract: The post-hoc stem reading sits inside a [C] sentence. "Neither is necessary" is unlabelled and already published.
  15. "Three unrelated lineages" for S216: One of the three is the blocked genome. GATE_0's caveat that T191S is "not isolable" in Anomia/Placuna is dropped, and so is GATE_4 §5's admission that loci were chosen
      after seeing the alignment.
  16. §9 "not contradicted": If both changes are required, §5c contradicts the model. Radopholus has neither.
  17. Scope: "All public mitogenomes" (§1) conflicts with §2. The title's pan-metazoan mtRF1a conflicts with GATE_1 §4.

  Undisclosed deviations (Q4)

  18. Missing from REPORT: T0, T3 and T5 by name, the number of tests run, origin ranges (only "11"), provenance stratification, and unresolved and wrong-locus counts by phylum.
  19. Loop remodelling: PREREG §2.4 says to score identity to the nearest canonical sister. Reports use human. Against its tylenchid sisters, whose loops are one residue shorter (Fig 2), Rotylenchulus is +1, so "no
      insertion" (§5e) depends on the comparator.
  20. PREREG §7: Post-hoc work must be dated under its own heading. The stem recount and the supported/dropped audits sit undated in §3, §5 and the abstract.
  21. mtRF1a skipped "by instruction" (GATE_4_PLAN): Globodera, a tylenchid C33 Stop species with a nuclear assembly, is missing from §3.5. It bears directly on "C33 plus S216 is untested".
  22. B1 "not testable": PREREG applies the ≥3-gene minimum to the sister taxon. The Rotylenchulus dossier has terminal calls for COX3 (TAA), CYTB (TAA) and ND1 (T--), with no TAG, against the prediction. Its
      termination is never discussed.

  Confounders (Q5)

  23. Closest relatives: Rotylenchulus' "10 closest" are Steinernema, Brugia and Caenorhabditis, not P. vulnus, a tylenchid Stop genome in the sweep. This affects the TAG analysis, A1 and the weak audit.
  24. Rotylenchulus NUMTs: No search is reported, although its nuclear assembly was used for mtRF1a.
  25. Labels: OR546184 (mislabelled) is near A. vulpes OR546140, and their submission series' runs SRR22396609 and SRR22396674 are close. I infer one batch. S. lanceolatoides and S. sector have identical trnY, and
      whole-genome identity is not reported.
  26. Non-independence: "All 9" Microdevario calls appear to count two RefSeq/source pairs twice, and four more are one study. C33 in Albula was known before the lock. T1's interval treats related species as
      independent.
  27. Detection of non-U33: Non-U33 always needs a second caller, but U33 at score ≥20 does not. Recall for non-U33 was never measured.
  28. Lycenchelys: "0 nuclear copies" follows exclusion of an 82 kb mitochondrial-like contig whose base 33 is not given.

  Mismatches (Q7)

  29. "91 of 92" (§5e): GATE_4 §4b says two records lack UAA, and also borealis 8 of 9. The text supports 89 or 90.
  30. Fig 1:
      - It labels Rotylenchulus "new", against §1 and §5g.
      - It omits the N33 state of Tyr lineages and the 28,566 U33 Stop genomes.
      - U33 relatives are missing, so 11 origins cannot be read from it.
  31. Fig 2 and §6: P. plicata is shown as Tyr without a mitogenome call, while R. robustus gets that caveat.

  C. Minor

  32. GATE_3 gives stem_audit_review as both 2,020 and 2,033. The confirmatory set is 32,624 in A13 and 32,623 later, unexplained.
  33. "No paper found" (GATE_3, 3B) versus Vovlas, and "no other A. vulpes mitogenome" (PREREG T6) versus OZ556087, are not logged in §3.7.
  34. REPORT states Lin's U36C = position 33 as fact; GATE_0 calls it inference. Jacob is "read" in REPORT but keyword-searched in GATE_0.
  35. Coloconger and Ophichthus dossiers blame a missing assembly for mtRF1a; the plan says skipped.
  36. §3.6's "about 28,575" should be 28,574.

  Q6: why "unconfirmed single-genome candidate"

  - Keep the locked count. PREREG cites GATE_2, which specifies model columns, and A13 shows the alternative was known at lock. Switching now changes the headline with results in view.
  - "Rule as written" is loaded. The text doesn't say "beside the loop". Call it a post-hoc alternative count.
  - Don't offer two readings for a reviewer to choose between. That makes the primary outcome a post-results choice.
  - Explain the asymmetry with A15. That fix was adopted because code contradicted the text and candidates did not change. Here the text is silent and the headline changes.
  - 5 of 5 doesn't validate the tRNA. The call has 9 nt before the anticodon arm and 27 after. The other three nematode dossiers have 24–26 and 14, so its boundaries look wrong.
  - Say the block is procedural. The UAA evidence is strong, and trnY quality says little about UAA decoding. It legitimately blocks the necessity claim; for the clade claim it is a design flaw.
  - Route to a claim: a dated amendment tested on new data, such as a second Rotylenchulus mitogenome.

  Check first against raw data: 6, 23 and 29.

---
