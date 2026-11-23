# Grigorovia jiainica: the four sites, for a lab to test

- Strain: **NRRL Y-48843** (ARS Culture Collection / NRRL, USDA, Peoria IL), BioSample SAMN20341879, BioProject PRJNA736342 (Y1000+), SRA run SRR16974510 (Illumina HiSeq 2500, paired). Other deposited strains of the species (markers only): CBS 10587, ATCC MYA-4454, SF1S05 (the type strain, SF1S05T = CBS 10587T, per Lee et al. 2008), TL1903.
- Mitogenome: GenBank TPA BK073075.1 (Sarhan et al. 2025), circular, 29,214 bp. cox1 and cob are on the minus strand of the record.
- Coordinates are 1-based on BK073075.1. Context is the coding (minus) strand, 5'->3', 30 nt either side, the extra base in brackets; the old-frame codons are split by spaces in the second line.

| gene | extra base (BK073075.1, 1-based) | ref. protein residue | old-frame codons at the site | read after skipping the extra base | reads with / without (duplicates removed) |
|---|---|---|---|---|---|
| cob | 26496 (minus strand) | ~137 of 386 | AGG/AGG | GGC | 1718 / 0 (1343 / 0) |
| cob | 26477 (minus strand) | ~143 of 386 | TGG/TGG | GGC | 1733 / 0 (1340 / 0) |
| cob | 26266 (minus strand) | ~213 of 386 | ACT/AGG | CTA as aligned; GGC if the skipped base is taken as the A of AGG (both placements give the same protein, Thr-Gly, in code 3) | 1785 / 0 (1381 / 0) |
| cox1 | 24711 (minus strand) | ~78 of 534 | GGA/AGG | GGC | 1440 / 0 (1151 / 0) |

**cob, extra base at 26496:**
```
GCAGCATTTTTAGGATACTGTTGTGTATAT[A]GGCCAAATGTCACATTGATGGCGCACTAGT
GCA GCA TTT TTA GGA TAC TGT TGT GTA TAT AGG CCA AAT GTC ACA TTG ATG GCG CAC TAG TTA
```
**cob, extra base at 26477:**
```
GTTGTGTATATAGGCCAAATGTCACATTGA[T]GGCGCACTAGTTATTACTAATTTATTCTCA
GTT GTG TAT ATA GGC CAA ATG TCA CAT TGA TGG CGC ACT AGT TAT TAC TAA TTT ATT CTC AGC
```
**cob, extra base at 26266:**
```
CAACATGGATCATCTAATCCATTAGGAATT[A]CTAGGCAATTTAGATAGATTACCAATGCAT
CAA CAT GGA TCA TCT AAT CCA TTA GGA ATT ACT AGG CAA TTT AGA TAG ATT ACC AAT GCA TAA
```
**cox1, extra base at 24711:**
```
TTCTTTATAGCAATACCAGCATTAATTGGA[A]GGCTTTGGTAATTATTTTTTACCATTGATA
TTC TTT ATA GCA ATA CCA GCA TTA ATT GGA AGG CTT TGG TAA TTA TTT TTT ACC ATT GAT AAT
```

Suggested test: RNA (RT-PCR/RNA-seq) across each site, to show whether the extra base is kept in the mRNA (translational shift) or removed (editing); and mitochondrial proteomics for peptides spanning the sites (Gly at the GGC position would support +1 reading).
