# Search log for the prior-art step (2026-11-21, US Eastern)

Everything searched for `PRIOR_ART.md`, by whom, where and with what. Queries are verbatim. Three
groups did searching: me (the main session), a literature-sweep subagent (its own log, verbatim, is
`sweep_agent_log.md`) and a subagent that looked for the Al Arab et al. 2017 full text. Five other
subagents only read saved files and did not search.

## A. Europe PMC REST API (my queries, via `code/lit_fetch.py search`)
Round 1 (≈19:13), to locate the leads named in DESIGN.md:
1. `TITLE:"Accurate annotation of protein-coding genes in mitochondrial genomes"` → Al Arab 2017
2. `AUTH:"Gornik" AND (frameshift OR mitochondrial) AND PUB_YEAR:2022` → Gornik 2022 (+ bioRxiv version)
3. `AUTH:"Andreu-Sanchez" AND (nad3 OR frameshift)` → Andreu-Sánchez 2021
4. `TITLE:"Recoding of translation in turtle mitochondrial genomes"` → Russell & Beckenbach 2008
5. `AUTH:"Haen" AND AUTH:"Lavrov"` → Haen 2014 (and Haen 2007, 2010)
6. `AUTH:"Szaboova" AND mitochondrial` → **target missed** (diacritics); found with query 12
7. `AUTH:"Masuda" AND Perkinsus` → Masuda 2010
8. `Polyrhachis AND frameshift` → **target missed**; found with query 13
9. `oyster AND mitochondrial AND frameshift` → **target missed**; found by web search W1
10. `Temperley AND "hungry codons"` → Temperley 2010 and the 2015 response
11. `(mtRF1 OR MTRF1) AND (AGA OR AGG) AND termination AND PUB_YEAR:[2019 TO 2026]` → Wakigawa 2026, Krüger 2024

Round 2 (≈19:15):
12. `DOI:"10.1093/femsyr/foy028"` → Szabóová 2018
13. `TITLE:"Polyrhachis" AND "cytochrome b"` → Beckenbach 2005
14. `(Crassostrea OR oyster) AND "cytochrome b" AND (frameshift OR "+1")` → nothing relevant
15. `DOI:"10.1186/1471-2164-9-33"` → Rosengarten 2008
16. `"translational frameshift" AND mitochondri* AND (survey OR genomes OR lineages) AND PUB_YEAR:[2010 TO 2026]` → noise
17. `(frameshift OR "frame shift") AND "mitochondrial genome" AND annotation AND (MITOS OR RefSeq OR MFannot)` → noise
18. `"+1 frameshift" AND mitochondrial AND (insect OR beetle OR wasp OR bee OR "true bug" OR ant)` → Towett-Kirui 2022
19. `mitochondrial AND frameshift AND (dinoflagellate OR apicomplexan OR "Chromera") AND (cox3 OR cox1)` → Gornik 2022, Wú 2024 preprint, Macher 2023
20. `"programmed frameshift" AND mitochondri* AND (review OR recoding)` → Russell & Beckenbach 2008 only
21. `Kazachstania AND mitochondrial AND (frameshift OR CTG)` → Wolters 2023

Audit round (19:30–19:55):
22. `(SRC:PPR) AND (mitochondri*) AND (frameshift OR "frame shift" OR "ribosomal slippage" OR quadruplet) AND PUB_YEAR:[2025 TO 2026]` → 15 shown, all human-disease or unrelated
23. `(mitochondri*) AND ("translational frameshift" OR "programmed frameshift" OR "+1 frameshift" OR "frameshift insertion") AND PUB_YEAR:[2024 TO 2026] NOT SRC:PPR` → noise (clinical "frameshift mutation" papers)
24. `(TITLE_ABS:"mitochondrial genome" OR TITLE_ABS:mitogenome OR TITLE_ABS:mitogenomes OR TITLE_ABS:"mtDNA") AND (TITLE_ABS:"translational frameshift" OR TITLE_ABS:"translational frameshifting" OR TITLE_ABS:"programmed frameshift" OR TITLE_ABS:"+1 frameshift" OR TITLE_ABS:"frameshift insertion" OR TITLE_ABS:"ribosomal frameshifting") AND PUB_YEAR:[2023 TO 2026]` → 3 hits: Yuan 2025 (plover *nad3*-174) and two records of the human mitochondrial mRNA structurome
25. `(Spizellomyces OR Rhizophydium OR Hyaloraphidium OR Monoblepharella OR chytrid*) AND (mitochondrial OR mtDNA) AND (genome OR "cytochrome b" OR cob)` → **failed, HTTP 503**; replaced by the PubMed fetch of Forget 2002 (cited by those records)
26. `"Arctica islandica" AND mitochondrial AND genome` → no frameshift paper
27. `TITLE:"Massive programmed translational jumping in mitochondria"` → Lang 2014
28. `(byp OR "translational bypassing" OR "translational jumping" OR "translational hopping") AND (mitochondri*)` → Lang 2014, Nosek & Tomáška 2025 review, O'Loughlin 2020 (phage T4)

Fetches by identifier (`code/lit_fetch.py fetch`, not searches): every paper with a `.json` in `lit/`.

## B. Web search tool (mine)
- W1 (standard): `oyster mitochondrial cytochrome b +1 frameshift Crassostrea programmed frameshift cytb` → Milbury & Gaffney 2005
- W2 (extended): `systematic survey of translational frameshifts in mitochondrial genomes across eukaryotes RefSeq annotation frameshift sites catalogue` → no survey; known case papers only
- W3 (extended): `mitochondrial genome +1 frameshift at recurring codon lineage-wide "frameshift" sites rare codon tRNA missing mitochondria new case 2023 2024 2025` → known cases only
- W4 (extended): `bioRxiv 2025 2026 mitochondrial genome translational frameshift +1 frameshift sites codon mitogenome preprint` → no 2025–2026 preprint on the topic
- Two earlier searches (on TilS literature) belong to another project and are not part of this log.

## C. NCBI E-utilities (mine)
- `code/audit/ncbi_lookup.py`: 11 nuccore queries for accessions of known cases (queries and results in
  `results/audit/ncbi_lookup.txt`).
- PubMed efetch, PMID 11861890 (Forget et al. 2002, the paper cited by the chytrid records) →
  `lit/forget_2002_pubmed.txt`.

## D. Forward citations (papers citing a key paper) — `code/audit/forward_citations.py`
Sources: Europe PMC citation network (`/MED/{pmid}/citations`, journals; no abstracts) and OpenAlex
(`filter=cites:`, includes preprints and theses; abstracts where available). Screen: regex
`frame-?\s?shift|ribosomal slippage|recoding|quadruplet|non-?triplet|\+1 (shift|insertion)|insertion editing|programmed (translational )?frameshift`
on title + abstract. EPMC-only hits were screened on titles alone. Every citing work is listed in
`lit/forward/<seed>_citing.tsv`; screened hits in `lit/forward/screened.tsv`.

| Seed | EPMC | OpenAlex | union | preprints | 2025–26 | screened hits |
|---|---|---|---|---|---|---|
| Al Arab 2017 | 74 | 106 | 178 | 10 | 126 | 1 |
| Gornik 2022 | 10 | 13 | 23 | 5 | 11 | 1 |
| Szabóová 2018 | 2 | 3 | 5 | 0 | 2 | 0 |
| Rosengarten 2008 | 34 | 55 | 89 | 3 | 0 | 18 |
| MitoZoa 2.0 (2012) | 40 | 56 | 95 | 2 | 7 | 0 |
| Lang 2014 (added) | 57 | 70 | 125 | 5 | 7 | 22 |
| Haen 2014 (added) | 23 | 26 | 47 | 2 | 4 | 12 |
| Masuda 2010 (added) | 23 | 32 | 55 | 4 | 1 | 14 |
| Russell & Beckenbach 2008 (added) | 72 | 99 | 170 | 9 | 10 | 50 |

What forward searching added beyond the papers already in hand:
- Atkins et al. 2016 review (NAR) → via its references, **Lang et al. 2014 (translational bypassing,
  *Magnusiomyces*)** and Nosek et al. 2015 — the most important miss of the first round.
- Wang et al. 2020 (*Zool Scr*): programmed frameshifts in *nad4* of three venerid bivalves and *cob* of one.
- Seligmann's papers (2011, 2012 ×2, 2019): computational "frameshifted/overlapping" mitochondrial genes
  and peptides read from "expanded codons" (contested; single species).
- Yuan 2025 (plover), Arborophila 2017 (partridges): further *nad3*-174 reports.
- Poehls et al. 2025 (bioRxiv, SLIPPERRS): frameshift sites across bacteria, not mitochondria.
- 2025–2026 preprints citing any seed: none reports a new mitochondrial frameshift lineage (titles of all
  are in the citing tables; abstracts screened where OpenAlex had them).

## E. Backward following
- Atkins 2016 refs 56, 64, 65 → Nosek 2015, Masuda 2010, Lang 2014.
- GenBank records NC_003052/NC_003053 (chytrids) → Forget 2002 (abstract mentions no *cob* frameshift).
- Wolters 2023 → Lang 2014 (byps).

## F. What was not searched
- Google Scholar, Web of Science, Scopus (none used by anyone).
- bioRxiv/medRxiv site search directly (preprints reached only through Europe PMC and OpenAlex).
- Non-English literature (e.g. CNKI), conference abstracts, patents.
- Theses except where OpenAlex listed them as citing works (titles only).
- Forward citations of Beckenbach 2005, Milbury & Gaffney 2005, Szabóová's later work, Andreu-Sánchez 2021,
  Wolters 2023, Kahn & Geller 2018, Wang 2020.
- All of GenBank's mitochondrial feature tables (only RefSeq 238 and the 2,695 fungal records were counted).
- Land-plant RNA editing; Myxozoa, Placozoa, Ctenophora specifically; plastid frameshifts beyond one sweep query.
