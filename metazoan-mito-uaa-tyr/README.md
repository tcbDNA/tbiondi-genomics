# Mitochondrial UAA and tRNA-Tyr position 33 in animals

In most animal mitochondria the codon UAA is a stop: the release factor mtRF1a reads it and ends translation. In a few
lineages UAA codes for tyrosine. Zhao et al. 2026 (`10.1186/s12915-026-02767-5`) propose that two changes drive the
switch:
- **U33C** in the tyrosine tRNA, tRNA-Tyr(GUA). Position 33 is the base just before the anticodon, and it is normally
  U.
- **T191S** in mtRF1a, which is encoded in the nucleus. The same residue is position 216 in the human protein
  (`NP_061914.3`), so it is written T216 or S216 below.

Lin et al. 2026 (bioRxiv `10.64898/2026.08.14.744788`) report a tRNA-Tyr change they call "U36C"; we infer that it is
position 33. They write that it is "absent in all species retaining the canonical code 5".

This project checks in public data whether the two changes go together with UAA = tyrosine. The codon and the tRNA
were read in 32,982 animal mitochondrial genomes. mtRF1a was read in a subset of species: some from the calibration
panel, some that the sweep flagged, and relatives of those.

The short answer: zero new reassigned clades under the locked rules, and C33 is not sufficient for UAA=Tyr under the
vertebrate mitochondrial code (code 2).

Start with [`reports/REPORT.md`](reports/REPORT.md). The analysis rules, fixed before the sweep, are in
[`PREREG.md`](PREREG.md).

The analysis and the write-up were done by Claude agents and checked by AI. They have not been reviewed by a human
expert, and there was no lab work.

## Terms

- **Calibration panel:** 360 genomes used to set and check the rules, 20 of them known to read UAA as tyrosine and 340
  not. They are never counted in the confirmatory totals.
- **UAA = tyrosine rule:** at least 5 internal UAA codons in at least 3 genes, at least 3 of them where other animals
  have a conserved tyrosine, a log-likelihood ratio of at least 3, and agreement from
  [Codetta](https://github.com/kshulgina/codetta), which infers genetic codes.
- **Stem count:** a genome counts as a new clade or a counterexample only if its tRNA-Tyr has at least 4 of 5 base
  pairs in the anticodon stem, among other checks.
- **Readable:** a genome with a confident tRNA-Tyr call that passes those checks. The call is a screen result, not proof
  that the locus is a working tRNA-Tyr.
- **Codes:** NCBI translation tables. Code 2 is vertebrate mitochondrial, code 4 is assigned to sponges, code 5 is
  invertebrate mitochondrial, and in code 14 UAA is tyrosine.

## Results

The paragraphs below are the Abstract of [`reports/REPORT.md`](reports/REPORT.md), word for word. Labels: **[C]**
confirmatory (rules in `PREREG.md`), **[E]** exploratory, **[P]** post-hoc (defined after the results were seen). Section
numbers (§) are those of REPORT.md.

We asked whether two changes proposed to drive mitochondrial UAA reassignment to tyrosine track the code in public
data: tRNA-Tyr(GUA) U33C (Lin et al. 2026; Zhao et al. 2026) and mtRF1a T191S, human T216 (Zhao et al. 2026). The
scope is 32,982 public mitogenomes under rules fixed in `PREREG.md`, and mtRF1a only in flagged species and sisters.

**[C]** Zero new reassigned clades under the locked rules. One non-panel genome, *Rotylenchulus reniformis*, met the
Tyr rule; it fails only the locked model-column count of anticodon-stem pairs in its tRNA-Tyr (3 of 5), so it is an
"unconfirmed single-genome candidate" in PREREG §3's terms. **[P]** A check whose criteria were written before mapping
finds a real mitochondrial U33 tRNA-Tyr locus there: 545 of 545 reads read U33, the assembled flanks are identical to
the record, and both callers' boundaries meet the structural criterion. The two boundaries differ, so the boundary is
unresolved, and about 17% of random stretches of this mitogenome would meet the acceptor-stem part of that criterion.
The genome stays unconfirmed by the locked rule. NCBI Taxonomy already assigns the genus *Rotylenchulus* code 14 (UAA
= Tyr), and the flatworm-code reading of its COI was reported for *R. borealis* (Vovlas et al. 2026).

**[C]** 17 of 28,583 readable UAA=Stop genomes (0.059% per genome, 95% CI 0.037–0.095%; 0.062% per taxon) have a
non-U33 tRNA-Tyr call: 14 with C33 and 3 with G33. "Stop" rests on the internal-UAA counter: Codetta returned "?" for
28,580 of the 28,583. **[P]** Under one evidence rule, 13 of the 17 have Tier 1–3 support, 3 are unconfirmed and 1 is
dropped. The C33 genomes are 6 with same-species independent support (4 lineages), 3 with reads from their own
specimen only, 2 lineage-backed only, 2 unconfirmed and 1 dropped. The G33 genomes are 1 with independent support and
1 with own-specimen reads (both *Craniella*, 1 genus), and 1 unconfirmed.

**[P]** C33 is not sufficient for UAA=Tyr under the vertebrate mitochondrial code (code 2). This rests on the 6
independently supported C33 genomes, all vertebrates. It needs C33 to be the dominant tRNA-Tyr allele with no second
functional copy, not zero U33 reads; the read and contig searches found no second copy in any of the six (§7.2.10).
The contig search detects copies at roughly half the mitochondrial depth or more. It used pre-repeat-resolution
contigs for *M. kubotai* and *A. glossodonta*, a whole-genome-amplified library from a different fish for *A.
glossodonta*, and, for *M. kubotai*, a library from a fish whose mitogenome differs from the counted record by 2.4%.
G33 is tolerated with UAA=Stop in one sponge genus, *Craniella* (2 species, code 4), with two caveats for *C.
infrequens*: its independent RNA reads match a nuclear copy as well as they match the mitogenome, and the upper bound
on a U33 share in its reads is above 4%. **[C]** Among 16,902 resolved UAA=Stop genomes with code 5, including 1,491
molluscs (489 bivalves), 2 have a C33 call, neither in a mollusc; **[P]** one is unconfirmed (*Toxocara*) and one was
dropped (*Globodera*), so none is confirmed. Lin et al.'s statement is limited to species "retaining the canonical
code 5", and no confirmed genome here contradicts it. **[C]** None of 77 resolved starfish has C33; **[P]** the
figures Li et al. 2019 cite for a starfish C show only hemichordates, and the one starfish in their accession table is
U33 here.

**Calibration [panel, descriptive; prior art]:** U33 occurs in UAA=Tyr lineages, and T216 occurs in UAA=Tyr lineages,
so neither change is necessary. Serine at 216 occurs in two calibration lineages (Anomiidae with Placunidae, and
*Cephalodiscus*) and, **[E]**, in the unconfirmed candidate *Rotylenchulus*; none of the 47 UAA=Stop species with both
a UAA call and a readable mtRF1a position 216 has S216 (34 Gate 2 background species, *Heteranomia squamula*,
*Pododesmus macrochisma* and 11 in Gate 4; two more background species show only the paralogue MTRF1 and are not
counted; §6). The panel pass rate for the molecular gates (UAA rule plus tRNA-Tyr gate and audit) is 14 of 20 known
positives, against 20 of 20 for the UAA rule alone, and it is zero for *Enigmonia* and *Radopholus*. A null on new
clades is therefore weak, and the panel does not test the novelty and replication parts of PREREG §3.

The 17 screen positives under the post-hoc evidence rule:

| Support | C33 (14) | G33 (3) |
|---|---|---|
| Same-species independent support (Tier 1) | 6 | 1 |
| Reads from their own specimen only (Tier 2) | 3 | 1 |
| Lineage-backed only (Tier 3) | 2 | 0 |
| Unconfirmed | 2 | 1 |
| Dropped: a second tRNA-Tyr in the genome has a different base | 1 | 0 |

Other limits are in REPORT §2.1 and §8, and departures from the rules are in §3.

## Files

The work was done in stages called gates, and each has a report.

| Path | Contents |
|---|---|
| `PREREG.md` | the rules, locked in commit `5aa5774` before the sweep |
| `reports/REPORT.md` | the report |
| `reports/GATE_*.md` | gate reports: source checks (0), inventory (1), calibration (2), sweep (3), exploratory runs (3B), mtRF1a plan and results (4) |
| `reports/COLD_REVIEW.md`, `COLD_REVIEW_2.md` | two AI reviews of the report; each point is answered in REPORT §11 and §11b |
| `reports/dossiers/` | one page for each of the 17 genomes, *Rotylenchulus*, and *Pratylenchus goodeyi* (a weaker, exploratory UAA case) |
| `reports/figures/` | figures 1–3 (SVG, PDF, PNG) |
| `results/mitogenome_master.tsv` | the per-genome table from the sweep |
| `results/gate2/`, `gate3/`, `gate4/`, `cold_review/` | outputs of each stage |
| `results/gate3_unverified/` | an exploratory run on records that GenBank flags UNVERIFIED |
| `results/posthoc2/` | post-hoc checks run for the second review (its own README lists them) |
| `results/audit/` | accession-level decisions (excluded, duplicate, fetched) and their reasons |
| `data/manifests/` | download ledger and checksums |
| `src/` | scripts, named by gate or topic |
| `env/*.lock` | conda environments: arm64 macOS, and x86-64 macOS for the HMMER 3.3.2 that Codetta uses |

The lock was made in a local git repository, without an external timestamp. Where a gate report and REPORT.md give
different numbers, REPORT.md is the current one.

## Rerunning

- Downloaded data are not in the repository. The ledger lists sources and checksums, and marks files that were
  downloaded by hand in a browser. The papers' full texts are not redistributed.
- Codetta (commit `863359e`) goes in `env/tools/codetta`.
- No single script reruns everything. Run scripts from the project folder; the shell scripts change to it themselves.
- `src/check_home_paths.py` fails if any file holds a path from someone's home folder.
- The sweep took about 16 hours on an Apple silicon Mac.
