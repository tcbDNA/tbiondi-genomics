# Test 5: the rest of the mitogenome (*G. jiainica* BK073075.1 and the congener *G. yakushimaensis* OR703946.1)

All numbers: `test5_tables.md` (`code/explore/grigorovia_h1h2.py mitogenome`).

**Methods.**
- Scan events and stop columns from the frozen pipeline (DIAMOND, frameshift mode).
- miniprot 0.18 with *S. cerevisiae* proteins, code 3 (`-T 3`); relaxed settings (`--outs=0.3 -p 0.1 -n 1 -l 3`) where the default found nothing.
- Longest open frame for VAR1.
- Local nucleotide alignment of *S. cerevisiae* 15S and 21S rRNA and RPM1 (RNase P RNA), from NC_001224.1.
- tRNAscan-SE 2.0.12 in organellar mode.

## Protein genes

| gene | *G. jiainica* | congener | notes |
|---|---|---|---|
| cox1 | full length (534/534 aa); **1 frameshift** (the 24711 site); no stop; no intron | full length; no frameshift or stop; 1 intron | |
| cob | full length (384/385); **3 frameshifts** (the cob sites); no stop; no intron | full length; none; 1 intron (bI2 type) | |
| cox2 | full length (251/251); no frameshift, stop or extra base | full length; none | |
| cox3 | full length (269/269); none | full length; none | |
| atp6 | full length (259/259); none | full length; none | |
| atp8 | full length (48/48); none (miniprot default; the scan's DIAMOND run missed this short gene) | full length; none (relaxed settings) | |
| atp9 | full length (76/76); none | full length; none | |
| VAR1 (rps3) | open frame of 259 codons, no stop (DIAMOND hits to VAR1/rps3; two in-alignment AT-rich inserts of 24–27 nt, as is common in VAR1) | open frame of 272 codons | VAR1 is too divergent for full-length alignment |
| other | orf146-like fragment at 11,006–11,198 with **one +1 shift (CGC) and one stop column**. This is not a core gene; it is the only other scan event in the genome | — | an intergenic ORF fragment, not interpreted |

**Every scan event in *G. jiainica*:** the four extra bases in cox1 and cob, plus one shift in the orf146-like fragment.
- No stop columns in any core gene.
- The only extra bases in protein genes are the four already known.

## RNA genes

| gene | *G. jiainica* (score, coverage of the *S. cerevisiae* sequence, identity) | congener |
|---|---|---|
| 15S rRNA (rns) | 763, 0.79, 0.82 | 681, 0.83, 0.84 |
| 21S rRNA (rnl) | 1000, 0.79, 0.83 | 1051, 0.72, 0.82 |
| RPM1 (RNase P RNA) | 107, 0.53, 0.82 | 106, 0.53, 0.79 |

All three are present in both genomes, at similar similarity to *S. cerevisiae*. These alignments test presence
only, not function.

## tRNAs (tRNAscan-SE, organellar mode; ARAGORN in `step_d_trnas.md`)

- Both genomes carry 23 tRNA genes with the same isotypes and anticodons. The one difference is **tRNA-Gly: GCC in
  *G. jiainica*, TCC in the congener.**
- No tRNA is marked as a pseudogene in either genome; neither has an intron.
- Infernal scores are similar overall: *G. jiainica* 26.7–64.9, congener 35.1–66.1. The lowest in both is tRNA-Cys:
  26.7 against 35.1.
- tRNA-Trp(UCA) appears as "Sup TCA" in both.
- Both genomes also have an entry "Sup TTA", which is not interpreted here.

## Read against the logged prediction

- **H2 predicted** further disabling changes, especially in cox2/cox3, or a decaying tRNA set.
- **Seen:** no other disabling change.
  - All other protein genes are full length without stops or frameshifts.
  - The rRNAs and RPM1 are present.
  - The tRNA set is complete and matches the congener's, except tRNA-Gly.
- **This favours H1, weakly** (as logged). The extra bases are confined to cox1 and cob.
- Under H2, a genome that had stopped needing complex III and IV would be expected to accumulate damage in cox2 and
  cox3 too, given time. None is seen. So either the loss is very recent, or cox1/cob are still translated.
