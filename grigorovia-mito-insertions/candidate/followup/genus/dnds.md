# Branch dN/dS (codeml, PAML 4.10.10; yeast mitochondrial code, icode = 2; F3x4; cleandata = 1)

Taxa: Abovina, Atelluris, Ghumatica, Gjiainica, Gtransvaal, Gyakushima, Mservazzii. CDS from miniprot exon structures (S. cerevisiae proteins); G. jiainica's four extra bases fall in the gaps between miniprot segments and are therefore excluded. Codon alignments: MAFFT L-INS-i on proteins, back-translated. Tree: neighbour joining on concatenated protein p-distances (unrooted), `((Mservazzii,(Abovina,Atelluris)),(Gyakushima,(Gtransvaal,Ghumatica)),Gjiainica);`.

Tests per Grigorovia lineage (its terminal branch as foreground):
- **A: different from background** — two-ratio model (model = 2) against one ratio (M0), 1 df.
- **B: still constrained** — foreground omega free against foreground omega fixed at 1 (neutral), 1 df. Rejecting omega = 1 with omega < 1 means purifying selection on that branch.
- **Power**: expected substitutions on the branch, N·dN and S·dS (two-ratio model). With few synonymous changes, omega on that branch cannot be estimated well and neither test can reject.

## cox1+cob (943 codons before gap removal)

One-ratio M0: omega = 0.033, lnL = -8412.33

| lineage | omega background | omega branch | A: 2ΔlnL, p | B (omega = 1): 2ΔlnL, p | N·dN | S·dS | dS (branch) |
|---|---|---|---|---|---|---|---|
| G. jiainica | 0.033 | 0.031 | 0.04, 0.83 | 27.75, 1.4e-07 | 63.3 | 228.1 | 0.832 |
| G. yakushimaensis | 0.032 | 0.044 | 1.21, 0.27 | 28.94, 7.5e-08 | 60.7 | 156.2 | 0.569 |
| G. humatica | 0.035 | 0.022 | 1.67, 0.2 | 27.92, 1.3e-07 | 31.5 | 159.7 | 0.582 |
| G. transvaalensis | 0.039 | 0.013 | 22.80, 1.8e-06 | 151.67, 7.5e-35 | 50.3 | 421.4 | 1.536 |

## cox2+cox3+atp6+atp8+atp9 (908 codons before gap removal)

One-ratio M0: omega = 0.038, lnL = -8496.82

| lineage | omega background | omega branch | A: 2ΔlnL, p | B (omega = 1): 2ΔlnL, p | N·dN | S·dS | dS (branch) |
|---|---|---|---|---|---|---|---|
| G. jiainica | 0.037 | 0.052 | 1.59, 0.21 | 31.59, 1.9e-08 | 146.3 | 236.3 | 1.136 |
| G. yakushimaensis | 0.037 | 0.046 | 0.42, 0.52 | 23.35, 1.4e-06 | 61.9 | 112.5 | 0.541 |
| G. humatica | 0.044 | 0.007 | 24.23, 8.6e-07 | 65.42, 6.1e-16 | 12.1 | 139.4 | 0.671 |
| G. transvaalensis | 0.044 | 0.015 | 14.65, 0.00013 | 79.10, 5.9e-19 | 30.4 | 165.6 | 0.796 |

## all seven (1851 codons before gap removal)

One-ratio M0: omega = 0.036, lnL = -16967.50

| lineage | omega background | omega branch | A: 2ΔlnL, p | B (omega = 1): 2ΔlnL, p | N·dN | S·dS | dS (branch) |
|---|---|---|---|---|---|---|---|
| G. jiainica | 0.035 | 0.043 | 0.94, 0.33 | 59.78, 1.1e-14 | 205.3 | 467.5 | 0.973 |
| G. yakushimaensis | 0.035 | 0.045 | 1.42, 0.23 | 51.79, 6.2e-13 | 123.9 | 267.4 | 0.556 |
| G. humatica | 0.040 | 0.013 | 20.93, 4.8e-06 | 93.09, 5e-22 | 43.4 | 319.9 | 0.666 |
| G. transvaalensis | 0.042 | 0.014 | 42.82, 6e-11 | 251.00, 1.6e-56 | 80.5 | 557.8 | 1.161 |

## cox1 (558 codons before gap removal)

One-ratio M0: omega = 0.037, lnL = -5048.69

| lineage | omega background | omega branch | A: 2ΔlnL, p | B (omega = 1): 2ΔlnL, p | N·dN | S·dS | dS (branch) |
|---|---|---|---|---|---|---|---|
| G. jiainica | 0.035 | 0.061 | 1.79, 0.18 | 8.73, 0.0031 | 43.7 | 88.7 | 0.512 |
| G. yakushimaensis | 0.035 | 0.053 | 1.45, 0.23 | 17.59, 2.7e-05 | 44.3 | 102.7 | 0.593 |
| G. humatica | 0.038 | 0.027 | 0.58, 0.44 | 14.34, 0.00015 | 19.5 | 89.5 | 0.517 |
| G. transvaalensis | 0.042 | 0.017 | 10.62, 0.0011 | 98.19, 3.8e-23 | 36.8 | 271.0 | 1.565 |

## cob (385 codons before gap removal)

One-ratio M0: omega = 0.028, lnL = -3346.90

| lineage | omega background | omega branch | A: 2ΔlnL, p | B (omega = 1): 2ΔlnL, p | N·dN | S·dS | dS (branch) |
|---|---|---|---|---|---|---|---|
| G. jiainica | 0.031 | 0.011 | 4.51, 0.034 | 22.76, 1.8e-06 | 20.0 | 186.1 | 1.820 |
| G. yakushimaensis | 0.027 | 0.033 | 0.15, 0.7 | 9.93, 0.0016 | 16.6 | 49.6 | 0.484 |
| G. humatica | 0.030 | 0.016 | 1.26, 0.26 | 14.63, 0.00013 | 11.9 | 72.0 | 0.703 |
| G. transvaalensis | 0.035 | 0.009 | 12.74, 0.00036 | 52.85, 3.6e-13 | 13.8 | 152.7 | 1.495 |

## Notes

- Gyakushima cox1: 2 segments, gaps between them [1612] nt (introns, or skipped bases where < 5 nt)
- Gyakushima cob: 2 segments, gaps between them [1491] nt (introns, or skipped bases where < 5 nt)
- Ghumatica cox1: 3 segments, gaps between them [2503, 2183] nt (introns, or skipped bases where < 5 nt)
- Ghumatica cob: 2 segments, gaps between them [1496] nt (introns, or skipped bases where < 5 nt)
- Gtransvaal cox1: 3 segments, gaps between them [1698, 3167] nt (introns, or skipped bases where < 5 nt)
- Mservazzii cox1: 4 segments, gaps between them [475, 1435, 1230] nt (introns, or skipped bases where < 5 nt)
- Mservazzii cob: 2 segments, gaps between them [1676] nt (introns, or skipped bases where < 5 nt)
