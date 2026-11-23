# Test 1: candidate placements of each extra base, scored against OR703946

Window +-30 nt on the coding strand; local alignment, match 2, mismatch -3, gap -5/-2. A site is comparable when the alignment of the window (scan placement removed) is gapless and spans window positions 24-36 around the site.

## cox1 24711

Residues: codons 21-38 of the corrected window, read in the frame of the corrected gene (codon starts at multiples of 3; the codon at 30 is the one right after the scan's extra base). 'stops' = stop codons anywhere in the 60-nt corrected window.

| offset | base removed | score | mismatches | gapless | codons 21-38 | residues | stops |
|---|---|---|---|---|---|---|---|
| +1 | G | 85 | 7 | True | TTA ATT GGA AGC TTT GGT | LIGSFG | 0 |
| +2 | G | 85 | 7 | True | TTA ATT GGA AGC TTT GGT | LIGSFG | 0 |
| +3 | C | 85 | 7 | True | TTA ATT GGA AGG TTT GGT | LIGRFG | 0 |
| -1 | A | 80 | 8 | True | TTA ATT GGA GGC TTT GGT | LIGGFG | 0 |
| +0 | A | 80 | 8 | True | TTA ATT GGA GGC TTT GGT | LIGGFG | 0 |
| -3 | G | 75 | 9 | True | TTA ATT GAA GGC TTT GGT | LIEGFG | 0 |
| -2 | G | 75 | 9 | True | TTA ATT GAA GGC TTT GGT | LIEGFG | 0 |
| congener | | | | | TTA ATT GGT AGT TTT GGT | LIGSFG | |

## cob 26496

Residues: codons 21-38 of the corrected window, read in the frame of the corrected gene (codon starts at multiples of 3; the codon at 30 is the one right after the scan's extra base). 'stops' = stop codons anywhere in the 60-nt corrected window.

| offset | base removed | score | mismatches | gapless | codons 21-38 | residues | stops |
|---|---|---|---|---|---|---|---|
| -1 | T | 41 | 7 | True | TGT GTA TAA GGC CAA ATG | CV*GQM | 1 |
| +0 | A | 41 | 7 | True | TGT GTA TAT GGC CAA ATG | CVYGQM | 0 |
| +3 | C | 41 | 7 | True | TGT GTA TAT AGG CAA ATG | CVYRQM | 0 |
| -2 | A | 36 | 8 | True | TGT GTA TTA GGC CAA ATG | CVLGQM | 0 |
| +1 | G | 36 | 8 | True | TGT GTA TAT AGC CAA ATG | CVYSQM | 0 |
| +2 | G | 36 | 8 | True | TGT GTA TAT AGC CAA ATG | CVYSQM | 0 |
| -3 | T | 34 | 5 | False | TGT GTA ATA GGC CAA ATG | CVMGQM | 0 |
| congener | | | | | TGC GTA TAC GGG CAA ATA | CVYGQM | |

## cob 26477

Residues: codons 21-38 of the corrected window, read in the frame of the corrected gene (codon starts at multiples of 3; the codon at 30 is the one right after the scan's extra base). 'stops' = stop codons anywhere in the 60-nt corrected window.

| offset | base removed | score | mismatches | gapless | codons 21-38 | residues | stops |
|---|---|---|---|---|---|---|---|
| +3 | C | 51 | 1 | True | TCA CAT TGA TGG GCA CTA | SHWWAT | 0 |
| -3 | T | 49 | 1 | True | TCA CAT GAT GGC GCA CTA | SHDGAT | 0 |
| -2 | G | 49 | 1 | True | TCA CAT TAT GGC GCA CTA | SHYGAT | 0 |
| -1 | A | 49 | 1 | True | TCA CAT TGT GGC GCA CTA | SHCGAT | 0 |
| +0 | T | 49 | 1 | True | TCA CAT TGA GGC GCA CTA | SHWGAT | 0 |
| +1 | G | 49 | 1 | True | TCA CAT TGA TGC GCA CTA | SHWCAT | 0 |
| +2 | G | 49 | 1 | True | TCA CAT TGA TGC GCA CTA | SHWCAT | 0 |

## cob 26266

Residues: codons 21-38 of the corrected window, read in the frame of the corrected gene (codon starts at multiples of 3; the codon at 30 is the one right after the scan's extra base). 'stops' = stop codons anywhere in the 60-nt corrected window.

| offset | base removed | score | mismatches | gapless | codons 21-38 | residues | stops |
|---|---|---|---|---|---|---|---|
| -2 | T | 75 | 9 | True | TTA GGA ATA CTA GGC AAT | LGMTGN | 0 |
| -1 | T | 75 | 9 | True | TTA GGA ATA CTA GGC AAT | LGMTGN | 0 |
| +0 | A | 75 | 9 | True | TTA GGA ATT CTA GGC AAT | LGITGN | 0 |
| -3 | A | 70 | 10 | True | TTA GGA TTA CTA GGC AAT | LGLTGN | 0 |
| +1 | C | 70 | 10 | True | TTA GGA ATT ATA GGC AAT | LGIMGN | 0 |
| +3 | A | 70 | 10 | True | TTA GGA ATT ACT GGC AAT | LGITGN | 0 |
| +2 | T | 65 | 11 | True | TTA GGA ATT ACA GGC AAT | LGITGN | 0 |
| congener | | | | | TTA GGT ATC CTT GGT AAT | LGITGN | |

