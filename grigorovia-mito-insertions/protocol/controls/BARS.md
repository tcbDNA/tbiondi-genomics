# Control bars (re-run after DESIGN.md amendment 4 and the PREREG log entry of 07:10)

## Headline sensitivity (unconditional)

Window rule, window cap = 4 nt (smallest width holding >= 95% of clean planted shifts); alpha_sense = 0.001, alpha_stop = 0.0032; split-hit shifts do not vote.

| planted sites at one codon | sense: genome flagged | stop: genome flagged | sense, all sites counted | stop, all sites counted |
|---|---|---|---|---|
| 2 | 0.317 (n=240) | 0.342 (n=120) | 1.000 (n=76) | 0.976 (n=42) |
| 3 | 0.483 (n=120) | 0.575 (n=120) | 1.000 (n=22) | 1.000 (n=20) |
| 5 | 0.717 (n=120) | 0.767 (n=120) | 1.000 (n=12) | 1.000 (n=8) |

"All sites counted" = every planted site detected as a clean shift voting for the planted codon (the conditional B5 definition). The unconditional columns are the headline.

### Planted shifts: detected in the gene, at the site, and at the right codon (share of planted)

| planted | group | n | in gene (same direction, anywhere) | at site (±9 nt) | clean at site | counted at right codon |
|---|---|---|---|---|---|---|
| random positions (planted null) | animals | 14355 | 0.906 | 0.886 | 0.680 | 0.564 |
| random positions (planted null) | fungi | 6054 | 0.832 | 0.816 | 0.650 | 0.528 |
| random positions (planted null) | protists | 4251 | 0.903 | 0.876 | 0.655 | 0.553 |
| random positions (planted null) | all | 24660 | 0.887 | 0.867 | 0.668 | 0.553 |
| same-codon trial sites | animals | 1452 | 0.894 | 0.830 | 0.568 | 0.470 |
| same-codon trial sites | fungi | 396 | 0.937 | 0.929 | 0.765 | 0.654 |
| same-codon trial sites | protists | 792 | 0.891 | 0.856 | 0.636 | 0.539 |
| same-codon trial sites | all | 2640 | 0.900 | 0.853 | 0.618 | 0.518 |

### Reasons, same-codon trial sites (share of sites)

| group | sites | counted (clean, votes for its codon) | clean_other_codon | clean_wide_window | not_aligned | aligned_no_event | masked_insert_edge | homopolymer | gene_end | gene_start | flank_fail | paired |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| animals | 1452 | 0.470 | 0.082 | 0.017 | 0.070 | 0.100 | 0.000 | 0.140 | 0.000 | 0.000 | 0.113 | 0.009 |
| fungi | 396 | 0.654 | 0.078 | 0.033 | 0.038 | 0.033 | 0.000 | 0.068 | 0.005 | 0.000 | 0.091 | 0.000 |
| protists | 792 | 0.539 | 0.076 | 0.021 | 0.077 | 0.067 | 0.003 | 0.159 | 0.000 | 0.000 | 0.051 | 0.008 |
| all | 2640 | 0.518 | 0.080 | 0.020 | 0.067 | 0.080 | 0.001 | 0.135 | 0.001 | 0.000 | 0.091 | 0.007 |

### Reasons, random planted shifts (planted null) (share of sites)

| group | sites | counted (clean, votes for its codon) | clean_other_codon | clean_wide_window | not_aligned | aligned_no_event | masked_insert_edge | homopolymer | gene_end | gene_start | flank_fail | paired |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| animals | 14355 | 0.564 | 0.099 | 0.018 | 0.070 | 0.044 | 0.000 | 0.156 | 0.000 | 0.000 | 0.046 | 0.004 |
| fungi | 6054 | 0.528 | 0.110 | 0.012 | 0.151 | 0.033 | 0.004 | 0.114 | 0.000 | 0.000 | 0.046 | 0.002 |
| protists | 4251 | 0.553 | 0.091 | 0.012 | 0.085 | 0.038 | 0.001 | 0.177 | 0.000 | 0.000 | 0.040 | 0.004 |
| all | 24660 | 0.553 | 0.100 | 0.015 | 0.092 | 0.041 | 0.001 | 0.149 | 0.000 | 0.000 | 0.045 | 0.003 |

## Scan gates (all must pass for the scan to start)

| Gate | Value | Target | Result |
|---|---|---|---|
| B3/B7 null genomes (120) + homopolymer and random-indel copies: genome flags | 0 of 360 {}  | 0 | PASS |
| B9b intron edges with a false clean shift within 90 nt | 0.0018 (n=17742; via split 0) | <= 0.02 | PASS |
| B6 Aphrocallistes vastus: flags | +1TGG k=2/2 genes=2 p=4.0e-04 | must flag | PASS |
| B6 Perkinsus marinus: flags | +1AGG k=4/7 genes=2 p=2.7e-06 | must flag | PASS |
| B6 Perkinsus olseni: flags | +1AGG k=4/7 genes=2 p=2.7e-06 | must flag | PASS |
| B6 Huiozyma sinensis (CTGY): flags | +1CTG k=2/3 genes=2 p=8.5e-06 | must flag | PASS |
| B6 Monosporozyma servazzii (CTGY): flags | +1CTG k=2/3 genes=2 p=8.5e-06 | must flag | PASS |
| B6 Magnusiomyces capitatus (bypass): bypass-like inserts | 65 | > 0 | PASS |

## Bars (a miss is logged in PREREG.md as a limit)

| Bar | Value | Target | Result |
|---|---|---|---|
| B1 planted single shifts, clean recall in aligned genes (homopolymers excluded) | 0.901 (n=17796; detected 0.963) | >= 0.90 | PASS |
| B5 two sites at one rare sense codon, flagged, when both counted (conditional) | 1.000 (n=76) | >= 0.90 | PASS |
| B2 RefSeq-annotated shifts detected at the site, ±9 nt (nad3) | 0.756 (n=619; clean 0.627; via split 0) | in gene anywhere 0.989 | in the voting window, right codon 0.590 | >= 0.80 (at site) | **MISS** |
| B2 RefSeq-annotated shifts detected at the site, ±9 nt (other) | 0.788 (n=66; clean 0.500; via split 0) | in gene anywhere 0.894 | in the voting window, right codon 0.409 | >= 0.80 (at site) | **MISS** |

## Stated limits (fungal introns)

| Item | Value | Status | |
|---|---|---|---|
| B9a annotated core-gene introns inferred as inserts | 0.619 (n=8871) | limit (was >= 0.80) | — |
| B9c shifts planted 31-90 nt from intron edges: clean recall | 0.279 (n=2101) | limit (was >= 0.80) | — |
| near-edge planted shifts inside an aligned hit | 0.438 | limit | — |
| near-edge outcomes | not_aligned 1180, clean 587, flank_fail 161, homopolymer 90, paired 36, aligned_no_event 31, masked_insert_edge 15, gene_end 1 | limit | — |
| fungal intron genomes: per-genome core-CDS coverage | median 0.795, 10th pct 0.647 (n=663) | limit | — |
| fungal intron genomes: masked coding share | median 0.030, 90th pct 0.061, max 0.288 | limit | — |

## Misplacement and pool sites

- Placement offsets of annotated shifts detected in the gene but not at the site (n=155): median 75 nt, 95th pct 120 nt, max 1531 nt; 153 of 155 within 150 nt (= the 50-residue pool tolerance). Offsets beyond 150 nt are usually a different shift elsewhere in the gene.
- Bird and turtle nad3 clean shifts in the annotated set, distinct pool sites per pool (108 family/order pools): old 2% rule max 6 (Charadriiformes), pools with >= 3 sites 13; new 50-residue rule max 2 (Apodiformes), pools with >= 3 sites 0.
- Exploratory baseline: null sequences with clean shifts in >= 3 different core genes: {'rnd': 3} (of 120 each).

## Notes

- Split-hit probation: split events: {'plant_null': 0, 'null_genomes': 0, 'sites': 0}; planted shifts recovered clean only via a split event: 0 (planted null) + 0 (same-codon sites); null sequences that flag only with split events voting: 0 []. Decision: split-hit shifts do NOT vote (exploratory tally only).
- Pre-written alpha rule applied to this run (record only; alpha stays 0.001 / 0.0032): alpha_sense would be 0.001, alpha_stop 0.0032.
- Expected chance flags across the scan (20343 genomes x 2 directions), estimated from the 120 null genomes' clean-shift counts weighted by phylum (group mean where a phylum is not sampled): 0.00. Clean core shifts per genome-direction in the null originals: {0: 233, 1: 7}. Upper bound (ignores the >= 2 genes condition). Recomputed exactly from the scan's own counts after the scan.
