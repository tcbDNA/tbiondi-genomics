# Literature and database search for *Grigorovia jiainica* (step g; fresh agent, 2026-11-22)

This file condenses the reports of two fresh-agent searches. Each source is listed with how deeply it was read.
Quotes are short and cited.

**Names searched.**
- *Grigorovia jiainica*.
- *Kazachstania jiainica* (Lee et al. 2008) and *Kazachstania jiainicus*.
- *Kazachstania* sp. SF1S05.
- *Pachytichospora jiainica* (Malimas et al. 2023, Jxiv preprint).
- The congeners *G. transvaalensis*, *G. humatica* and *G. yakushimaensis*, under both genus names.

**Terms combined with the names.** frameshift, frame shift, bypassing, jumping, hopping, slippage, RNA editing,
insertion, extra nucleotide, pseudogene, genetic code, reassignment, mitochondrial genome/mtDNA, mitochondri*.

**Sources searched.**
- Europe PMC: a grid of 10 name forms × 14 term groups, 140 queries. 126 returned 0, and every non-zero hit was
  checked. The grid is in `grigorovia_epmc_query_grid.json`.
- PubMed, OpenAlex, NCBI Nucleotide and general web search.

## Verdict

No source found describes +1 shifts at AGG (or TGG) in *G. jiainica*, under any of its names. No published report of
+1 frameshifting at AGG or TGG in a yeast mitochondrion was found.

The nearest case in a fungal mitochondrion is a database annotation with no paper describing it. The *cob* record of
the chytrid *Rhizophydium* sp. 136 (NC_003053.1) carries `/ribosomal_slippage` at an AGG followed by C. The paper the
record cites (Forget et al. 2002) does not mention a frameshift.

## Sources that mention the species or genus in a mitochondrial or translation context

- **Sarhan et al. 2025, *Sci Data*** (keyword search of the full text). The species is not named, "frameshift" does
  not appear, and the genomes were annotated with MFannot. The record BK073075.1 has no gene features.
- **Wolters et al. 2023, *Front Microbiol*** (keyword windows in the full text). No *Grigorovia* species is named. It
  reports frameshifts in "several Kazachstania species" consistent with the +1C mechanism, encoded as single-base
  introns.
  - Its data include the congener *G. yakushimaensis* ("Kazachstania yakushimaensis yHAB166"). The annotation
    files in its figshare deposit (not the GenBank record, which has no CDS features) read cox1 and cob straight
    through. That genome is GenBank OR703946, analysed here as the held-out relative.
- **Liu et al. 2024, *Persoonia*** (keyword search of the PDF text). Uses the NRRL Y-48843 nuclear genome for genus
  delimitation. No frameshift content.
- **Groenewald et al. 2023, *Stud Mycol*** (keyword windows). Taxonomy only ("Grigorovia is nested within
  Kazachstania").
- **Malimas et al. 2023, Jxiv preprint** (keyword search of the full text). Taxonomy only; transfers the species to
  *Pachytichospora*.
- **Older partial mitochondrial sequences of the sister species** (records and abstracts only). No frameshift is
  mentioned.

## +1 frameshifting at AGG or TGG elsewhere

**Fungal mitochondria.**
- *Rhizophydium* sp. 136 cob: a database annotation only (above).
- *Spizellomyces punctatus* cob: a record note, with the site not located.
- *Monosporozyma* (formerly *Kazachstania*) *servazzii* cob (Langkjaer et al. 2003; short excerpts read): the
  authors model its sites as +1 T at CTG-TGG, not as AGG.

**Yeast nuclear genes, not mitochondria.**
- The well-known CUU-AGG-C +1 site in Ty1 and a few nuclear genes (Belcourt & Farabaugh 1990; Farabaugh 2006;
  abstracts).
- Fenton et al. 2025 (*NAR*; keyword windows) relate it to an imbalance between the tRNAs for AGG (Arg) and
  GGC (Gly).

## Respiration (follow-up)

**Lee et al. 2008, species description** (short verbatim excerpts read).
- "Only glucose and galactose are assimilated."
- Glycerol, ethanol, DL-lactate, succinate and citrate are listed among the carbon sources not assimilated. Acetate
  was not tested as a carbon source.

**Opulente et al. 2018, *BMC Biol*, Table S10** (spreadsheet row read). "K. jiainica" is negative for ethanol,
glycerol, DL-lactate, succinate and citrate.

**Sister species.**
- *G. transvaalensis* is petite-positive, and its collection strain "already does not grow on nonfermentable
  substrate" (Fekete et al. 2007, *FEMS Yeast Res*, Table 1 footnote).
- *G. yakushimaensis*: negative for the same carbon sources in Opulente 2018.

**Not tested in any source found.** A petite test or respiration measurement for *G. jiainica*. *G. jiainica* is not
among the strains in the Y1000+ phenotype tables (Opulente et al. 2024).

## Limits of this search

- Several full texts were read only as short excerpts, through a web-fetch tool: Szabóová 2018, Langkjaer 2003,
  Forget 2002 and Lee 2008. Check quotes against the originals before citing them.
- Not read: the supplementary tables of Opulente 2024 beyond the phenotype tables; Gouliamova & Dimitrov 2020 (title
  only); the *Kazachstania* chapter of *The Yeasts* (2011; access refused).

## Strains (added at the second review, 2026-11-22)

- **SF1S05 is the type strain of the species.** Lee et al. (2008) give SF1S05T = CBS 10587T. The referee agent
  confirmed that the ATCC catalogue lists ATCC MYA-4454 as SF1S05, the type strain (catalogue entry, read directly).
- **The sequenced strain is NRRL Y-48843** (BioSample SAMN20341879). No record we found says whether it is a deposit of
  the same isolate as SF1S05.
