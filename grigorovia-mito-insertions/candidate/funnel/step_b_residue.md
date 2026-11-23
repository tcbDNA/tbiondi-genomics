# Step b: residue facing the codon read after the skipped base

## cox1 @ 24710: old-frame codons GGA/AGG; codon read after the skip: GGC = G
- reference residues at that position (10 hits): {'G': 8, 'F': 1, 'A': 1}; majority Gly: yes
- G (Saccharomycetaceae, cox1, bits 872); G (Saccharomycodaceae, cox1, bits 853); G (Wickerhamomycetaceae, cox1, bits 837); G (Phaffomycetaceae, cox1, bits 836); G (Saccharomycopsidaceae, cox1, bits 821); F (Pichiaceae, cox1, bits 761); G (Pneumocystaceae, cox1, bits 738); G (Capillidiaceae, cox1, bits 731); A (Debaryomycetaceae, cox1, bits 731); G (Microbotryaceae, cox1, bits 730)

## cob @ 26495: old-frame codons AGG/AGG; codon read after the skip: GGC = G
- reference residues at that position (10 hits): {'G': 10}; majority Gly: yes
- G (Saccharomycetaceae, cob, bits 639); G (Saccharomycodaceae, cob, bits 623); G (Saccharomycopsidaceae, cob, bits 588); G (Phaffomycetaceae, cob, bits 585); G (Debaryomycetaceae, cob, bits 585); G (Wickerhamomycetaceae, cob, bits 582); G (Pichiaceae, cob, bits 572); G (Metschnikowiaceae, cob, bits 570); G (Trichomonascaceae, cob, bits 531); G (Dipodascales incertae sedis, cob, bits 518)

## cob @ 26476: old-frame codons TGG/TGG; codon read after the skip: GGC = G
- reference residues at that position (10 hits): {'G': 10}; majority Gly: yes
- G (Saccharomycetaceae, cob, bits 639); G (Saccharomycodaceae, cob, bits 623); G (Saccharomycopsidaceae, cob, bits 588); G (Phaffomycetaceae, cob, bits 585); G (Debaryomycetaceae, cob, bits 585); G (Wickerhamomycetaceae, cob, bits 582); G (Pichiaceae, cob, bits 572); G (Metschnikowiaceae, cob, bits 570); G (Trichomonascaceae, cob, bits 531); G (Dipodascales incertae sedis, cob, bits 518)

## cob @ 26265: old-frame codons ACT/AGG; codon read after the skip: CTA = T
- reference residues at that position (10 hits): {'T': 9, 'N': 1}; majority Gly: no
- T (Saccharomycetaceae, cob, bits 639); T (Saccharomycodaceae, cob, bits 623); T (Saccharomycopsidaceae, cob, bits 588); N (Phaffomycetaceae, cob, bits 585); T (Debaryomycetaceae, cob, bits 585); T (Wickerhamomycetaceae, cob, bits 582); T (Pichiaceae, cob, bits 572); T (Metschnikowiaceae, cob, bits 570); T (Trichomonascaceae, cob, bits 531); T (Dipodascales incertae sedis, cob, bits 518)

**Sites where most reference proteins have glycine: 3 of 4.**

## cob @ 26265 under the AGG placement

The scan placed this shift at the A of ACT, where the codon read after the skip is CTA = Thr in code 3; the reference
proteins have Thr there (9 of 10). Skipping the A of the next codon, AGG, gives the same protein (Thr-Gly). Under that
AGG placement the codon after the skip is GGC.

Reference residues at the GGC column: **G 6, S 2, N 1** (9 hits). The column before has T 8, N 1.

So, under the AGG placement, glycine faces the GGC at all three AGG sites: 8 of 10, 10 of 10, and 6 of 9 reference
proteins.
