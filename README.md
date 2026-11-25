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
panel, some that the sweep flagged, and relatives of those. The short answer: no new lineage with UAA = tyrosine is
claimed, and C33 also occurs in a few vertebrates where UAA is a stop.

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
- **Readable:** a genome with a confident tRNA-Tyr call that passes those checks.
- **Codes:** NCBI translation tables. Code 2 is vertebrate mitochondrial, code 4 is assigned to sponges, code 5 is
  invertebrate mitochondrial, and in code 14 UAA is tyrosine.

## Results

### Confirmatory (rules in `PREREG.md`)
- One genome outside the panel, the nematode *Rotylenchulus reniformis*, met the UAA = tyrosine rule. Its tRNA-Tyr has
  3 of 5 stem pairs under the locked count, so **no new reassigned clade is claimed** and the genome is an unconfirmed
  candidate (a post-hoc recount is in REPORT §3.6 and §7.2.2). NCBI Taxonomy already assigns its genus code 14.
- 17 of 28,583 readable genomes with UAA as a stop (0.059%, 95% CI 0.037–0.095%) have a tRNA-Tyr without U33: 14 with
  C33 and 3 with G33. A genome is one GenBank or RefSeq record, so a species can appear more than once.

### Post-hoc (defined after the results were seen)
All 17 genomes were sorted by one rule, on how well the base is supported:

| Support | C33 (14) | G33 (3) |
|---|---|---|
| Another individual of the same species (Tier 1) | 6 | 1 |
| Reads from the same specimen only (Tier 2) | 3 | 1 |
| Related species only (Tier 3) | 2 | 0 |
| Unconfirmed | 2 | 1 |
| Dropped: a second tRNA-Tyr in the genome has a different base | 1 | 0 |

- C33 is not sufficient for UAA = tyrosine under code 2. This rests on the 6 Tier 1 C33 genomes, all vertebrates.
- G33 occurs with UAA as a stop in one sponge genus, *Craniella* (code 4). Its two species are the Tier 1 and Tier 2
  genomes.
- Under code 5, 2 of 16,902 genomes with UAA as a stop have a C33 call; one is unconfirmed and one was dropped. No
  confirmed genome contradicts the statement by Lin et al.

### Calibration panel (descriptive, known from the literature)
- Some lineages with UAA = tyrosine have U33, and some have T216, so neither change is necessary.
- S216 occurs in two calibration lineages (the bivalves Anomiidae with Placunidae, and *Cephalodiscus*) and in the
  *Rotylenchulus* candidate. None of the 47 species that have UAA as a stop in this study and a readable position 216
  has S216.

### Limits
The UAA rule finds all 20 known positives, which were in view when it was set. Only 14 of them also pass the tRNA-Tyr checks
that a new clade needs, so a new clade could have been missed. Other limits are in REPORT §2.1 and §8, and departures from the rules are in §3.

## Files

The work was done in stages called gates, and each has a report.

| Path | Contents |
|---|---|
| `PREREG.md` | the rules, locked in commit `5aa5774` before the sweep |
| `reports/REPORT.md` | the report |
| `reports/GATE_*.md` | gate reports: source checks (0), inventory (1), calibration (2), sweep (3), exploratory runs (3B), mtRF1a plan and results (4) |
| `reports/COLD_REVIEW.md` | an AI review of the report; each point is answered in REPORT §11 |
| `reports/dossiers/` | one page for each of the 17 genomes, *Rotylenchulus*, and *Pratylenchus goodeyi* (a weaker, exploratory UAA case) |
| `reports/figures/` | figures 1–3 (SVG, PDF, PNG) |
| `results/mitogenome_master.tsv` | the per-genome table from the sweep |
| `results/gate2/`, `gate3/`, `gate4/`, `cold_review/` | outputs of each stage |
| `results/gate3_unverified/` | an exploratory run on records that GenBank flags UNVERIFIED |
| `results/audit/` | accession-level decisions (excluded, duplicate, fetched) and their reasons |
| `data/manifests/` | download ledger (118.49 GB) and checksums |
| `src/` | scripts, named by gate or topic |
| `env/*.lock` | conda environments: arm64 macOS, and x86-64 macOS for the HMMER 3.3.2 that Codetta uses |

The lock was made in a local git repository, without an external timestamp. Where a gate report and REPORT.md give
different numbers, REPORT.md is the current one.

## Rerunning

- Downloaded data are not in the repository. The ledger lists sources and checksums, and marks files that were
  downloaded by hand in a browser. The papers' full texts are not redistributed.
- Codetta (commit `863359e`) goes in `env/tools/codetta`.
- No single script reruns everything. The nine shell scripts in `src/` start with a `cd` to the original machine's
  path, which has to be changed.
- The sweep took about 16 hours on an Apple silicon Mac.
