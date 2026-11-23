# Scan flags

Genomes scanned: 20343. Window cap 4 nt; alpha_sense = 0.001, alpha_stop = 0.0032; alpha_pool = 1e-05 (/4).
Genome-directions with >= 2 clean core shifts: 143. Expected chance flags under the rule, from the scan's own counts: 2.5 (upper bound). Observed single-genome flags: 16.
Shuffled-pool flags (seeds 1-3; family, order): [2, 1, 3, 0, 4, 0]. Pool-direction tests: 101.
Split-hit shifts do not vote (tally only). Exploratory: genomes with clean shifts in >= 3 different core genes: 23 (11 in known lineages).
Core-CDS coverage, RefSeq genomes with CDS: median 0.97, 10th pct 0.797.

## Single-genome flags in known lineages: 15
- Graptemys ouachitensis (NC_071766.1; Emydidae, Testudines, Chordata): +1 AGT k=2/2, genes 2, p=8.1e-04, E=1.8, same codon in related NC_031300.1 **STRONG**
- Tachybaptus novaehollandiae (NC_010095.1; Podicipedidae, Podicipediformes, Chordata): +1 AGT k=2/2, genes 2, p=4.0e-04, E=0.76, same codon in related NC_008140.1 **STRONG**
- Podiceps cristatus (NC_008140.1; Podicipedidae, Podicipediformes, Chordata): +1 AGT k=2/2, genes 2, p=4.0e-04, E=0.76, same codon in related NC_010095.1 **STRONG**
- Malaclemys terrapin terrapin (NC_031300.1; Emydidae, Testudines, Chordata): +1 AGT k=2/2, genes 2, p=8.1e-04, E=1.8, same codon in related NC_071766.1 **STRONG**
- Huiozyma sinensis (NC_036748.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=2/3, genes 2, p=8.5e-06, E=0.0066, same codon in related BK073131.1,BK073176.1,BK073185.1,BK073484.1,BK073486.1,NC_004918.1 **STRONG**
- Monosporozyma servazzii (NC_004918.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=2/3, genes 2, p=8.5e-06, E=0.0066, same codon in related BK073131.1,BK073176.1,BK073185.1,BK073484.1,BK073486.1,NC_036748.1 **STRONG**
- Perkinsus marinus (NC_068742.1; Perkinsidae, Perkinsida, Perkinsozoa): +1 AGG k=4/7, genes 2, p=2.7e-06, E=0.0016, same codon in related NC_068743.1 **STRONG**
- Perkinsus olseni (NC_068743.1; Perkinsidae, Perkinsida, Perkinsozoa): +1 AGG k=4/7, genes 2, p=2.7e-06, E=0.0016, same codon in related NC_068742.1 **STRONG**
- Aphrocallistes vastus (NC_010769.1; Aphrocallistidae, Sceptrulophora, Porifera): +1 TGG k=2/2, genes 2, p=4.0e-04, E=0.75 **STRONG**
- Grigorovia jiainica (nom. inval.) (BK073075.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 AGG k=3/4, genes 2, p=3.5e-05, E=0.036 **STRONG**
- Huiozyma naganishii (BK073131.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=3/4, genes 2, p=1.9e-08, E=5e-06, same codon in related BK073176.1,BK073185.1,BK073484.1,BK073486.1,NC_004918.1,NC_036748.1 **STRONG**
- Jamesozyma lodderae (BK073176.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=3/3, genes 3, p=4.8e-09, E=1.1e-06, same codon in related BK073131.1,BK073185.1,BK073484.1,BK073486.1,NC_004918.1,NC_036748.1 **STRONG**
- Kazachstania viticola (BK073185.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=6/7, genes 4, p=1.6e-16, E=1.1e-14, same codon in related BK073131.1,BK073176.1,BK073484.1,BK073486.1,NC_004918.1,NC_036748.1 **STRONG**
- Monosporozyma aerobia (BK073484.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=2/3, genes 2, p=8.5e-06, E=0.0066, same codon in related BK073131.1,BK073176.1,BK073185.1,BK073486.1,NC_004918.1,NC_036748.1 **STRONG**
- Monosporozyma solicola (BK073486.1; Saccharomycetaceae, Saccharomycetales, Ascomycota): +1 CTG k=3/4, genes 2, p=1.9e-08, E=5e-06, same codon in related BK073131.1,BK073176.1,BK073185.1,BK073484.1,NC_004918.1,NC_036748.1 **STRONG**

## Single-genome flags NOT in a known lineage: 1 (1 strong)
- Verticillium dahliae VdLs.17 (NW_009276970.1; Plectosphaerellaceae, Glomerellales, Ascomycota): -1 CAC k=3/29, genes 2, p=3.5e-04, E=0.61 **STRONG**

## Pool flags: 2 (2 in known lineages)
- family_pool Saccharomycetaceae (known): +1 CTG k=6/9 sites, genes 4, members 101, p=1.9e-15, E=1.2e-11 **STRONG**
- order_pool Saccharomycetales (known): +1 CTG k=6/9 sites, genes 4, members 102, p=1.9e-15, E=1.2e-11 **STRONG**

## Stop rule
Strong flags outside known lineages (funnel candidates; funnel not started without asking): 1 genome, 0 pool.
