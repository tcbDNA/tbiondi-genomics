# Grigorovia jiainica (BK073075.1): completeness and other data (exploratory)

## Completeness
- Record: circular, 29,214 bp (TPA_asm "mitochondrion, complete sequence"). The six closest relatives are 22,474-34,643 bp.
- Core genes aligned: atp6, atp9, cob, cox1, cox2, cox3. atp8 (about 48 aa) is not hit; it is also not hit in 3 of the
  6 relatives. No nad genes, as expected for this lineage. Also rps3 and orf146.
- tRNAs: ARAGORN found 20; tRNAscan-SE 2.0.12 (organellar mode, Infernal 1.1.5) found 23, including tRNA-Pro and
  tRNA-Cys that ARAGORN missed. Same isotype set as all six relatives. tRNA-Arg anticodon UCU (reads AGA and, by
  wobble, AGG) as in all relatives. The only tRNA-Gly has anticodon GCC (reads GGC/GGU); the relatives' tRNA-Gly has UCC.
- Conclusion: the assembly looks complete. Completeness gives no reason to expect local assembly errors, though it
  cannot exclude them.

## In-frame AGG usage (owed)
- Grigorovia 0 per 1,000 aligned core codons (0/1,748); relatives 0, 0, 0, 0, 0 and 1.01 per 1,000 (K. bromeliacearum,
  2/1,981). AGA about 19-22 per 1,000 in all seven.
- So "0 in-frame AGG" is normal for this group and does not distinguish Grigorovia (step c uninformative). An extra A
  placed before a GGC glycine codon necessarily creates an AGG in the old frame.

## Other data on the species or genus (NCBI, 2026-11-22)
- Former name: Kazachstania jiainica Lee et Liu 2008 (moved to Grigorovia by Gouliamova & Dimitrov 2020, with three
  other former Kazachstania species: G. transvaalensis, G. humatica, G. yakushimaensis). Also "Kazachstania jiainicus",
  "Kazachstania sp. SF1S05".
- SRA, searched by name for each species (corrected 2026-11-22). A first search by the genus taxid (txid2917336)
  found only G. jiainica, because NCBI files the other three species under a different genus node (txid2862416).
  All entries are DNA; there is no RNA-seq for any Grigorovia species.
  - G. jiainica: SRR16974510 (WGS, Illumina HiSeq 2500 PE, strain NRRL Y-48843, yHAB144, UW-Madison / Y1000+). Its
    nuclear assembly GCA_030571715.1 (2023, contig N50 20.6 kb) comes from the same reads.
  - G. yakushimaensis: SRR6475932 and SRR16974690 (WGS, HiSeq 2500; BioSample SAMN08343276, the source of OR703946).
  - G. transvaalensis: SRR6475934 and SRR16974369 (WGS, HiSeq 2500; SAMN08343274) and ERR9706995 (WGS, NovaSeq 6000).
  - G. humatica: SRR16974597 (WGS, HiSeq 2500; SAMN20341876) and two MiSeq amplicon runs.
- Nucleotide records for G. jiainica: the TPA mitogenome and rRNA/ITS markers only (strains SF1S05, CBS 10587,
  ATCC MYA-4454, TL1903, 25-RM-435-1). No other mitogenome of it.
- None of the other three Grigorovia species is in the scan set.
- So no independent assembly or RNA of G. jiainica exists. The reads behind the assembly are the only available
  test. WGS reads exist for all three sister species; only G. yakushimaensis has a deposited mitogenome (OR703946).
