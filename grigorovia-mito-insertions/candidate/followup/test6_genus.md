# Test 6: the four *Grigorovia* mitogenomes (exploratory; predictions in PREREG)

All numbers: tables in `genus/` (`code/explore/grigorovia_genus.py`): `mtreads.md`, `pick.md`, `bait.md`, `confirm_*.md`, `genes_in_contigs.md`, `compare.md`, `dnds.md`, `dnds_profile.md`, `refine_cob_junction.md`.

## Data

| species | run | strain | size | used? |
|---|---|---|---|---|
| *G. humatica* | SRR16974597 (WGS, 2×251) | NRRL Y-48839 | 0.19 GB fastq.gz (0.38 GB .sra) | yes, all reads |
| *G. transvaalensis* | SRR16974369, and SRR6475934 (the same BioSample, read count and bases: one dataset listed twice) | NRRL Y-17245 = CBS 2186, type strain | 0.40 GB fastq.gz (0.65 GB .sra) | **no: almost no mtDNA reads** (below) |
| *G. transvaalensis* | ERR9706995 (WGS, 2×151) | CLIB 448 | 0.83 GB fastq.gz (0.58 GB .sra) | yes |

- About 1.6 GB downloaded in total (NCBI's AWS mirror), well under the 20 GB limit.
- *G. yakushimaensis* (OR703946) and *G. jiainica* (BK073075) are deposited assemblies.
- Mitochondrial contigs in the deposited *nuclear* assemblies: none, or single NUMT-like hits only.

**The type strain of *G. transvaalensis* has lost its mtDNA.**
- Of 200,000 reads, 4 hit mitochondrial proteins: 0.002%, likely nuclear copies (DIAMOND blastx, code 3).
- For comparison: *G. humatica* 991 (0.5%), *G. jiainica* 2,933 (1.5%), *G. transvaalensis* CLIB 448 1,326 (0.7%).
- No high-coverage AT-rich contig was assembled.
- This matches Fekete et al. 2007: the collection strain "already does not grow on nonfermentable substrate".
- So that strain's respiration-negative phenotype is a property of the culture, not the species. Strain CLIB 448
  has a mitogenome.

## Assembly and confirmation

**Method.**
- SPAdes 4.3.0 (`--only-assembler`), on all *G. humatica* reads and on the first 2 million read pairs of CLIB 448.
- Mitochondrial contigs: hit to a mitochondrial protein, or aligned to OR703946, or GC < 25%; and coverage ≥ 3× the
  nuclear median.
- The first CLIB 448 assembly was fragmented (27 contigs, 33 kb). Its mitochondrial reads were then baited from all
  6 million pairs (minimap2, ≥ 90% identity over ≥ 80 nt, both mates) and reassembled de novo (SPAdes with error
  correction), in two rounds.
- Confirmation: all reads mapped back with minimap2; depth at MAPQ ≥ 20 (`genus/confirm_*.md`).

| species | mitochondrial contigs | total | median depth | positions ≥ 10× | core genes complete |
|---|---|---|---|---|---|
| *G. humatica* | 5 (largest 53,413 nt, carrying all seven core genes) | 62,112 nt | 161 (main contig) | 100% (main contig) | cox1, cob, cox2, cox3, atp6, atp8, atp9 |
| *G. transvaalensis* CLIB 448 | 7 (largest 15,361 nt) | 37,049 nt | 279–548 | 100% | all seven, each complete within one contig |

Neither assembly is circular, as AT-rich spacers and GC clusters usually break yeast mitogenome assemblies. The CLIB
448 total (37 kb) is short of a full mitogenome, so tRNA and intergenic counts for it are incomplete.

## Results

**1. Frozen scan pipeline, unchanged.**

| species | clean core shifts | +1 at AGG | flags | stop columns in core genes |
|---|---|---|---|---|
| *G. jiainica* | 4 | 3 | 1 | 0 |
| *G. yakushimaensis* | 0 | 0 | 0 | 0 |
| *G. humatica* | 0 | 0 | 0 | 0 |
| *G. transvaalensis* | 0 | 0 | 0 | 0 |

- None of the three congeners carries an extra base, or any other disabling mutation, in a core gene.
- At the four *G. jiainica* sites, **all three lack the extra base: 12 of 12 comparisons** (`genus/compare.md`).
  - Method: *G. jiainica*'s ±60-nt window is aligned to each species' spliced CDS (introns removed). The call is the
    net length difference across 15 nt either side, which does not depend on where the aligner puts a gap in a repeat.

**2. cox1 and cob introns**, mapped onto *G. jiainica* coordinates; distance from the nearest site:

| species | cox1 | cob |
|---|---|---|
| *G. jiainica* | none | none |
| *G. yakushimaensis* | 1 (742 nt from the cox1 site) | 1 (bI2 type): **2–3 nt from cob 26477** (3 nt under the group I splice rule), 20 nt from 26496 |
| *G. humatica* | 2 (695 nt away; the other far) | 1: the same intron, **2–3 nt from 26477** |
| *G. transvaalensis* | 2 (448 nt away; the other far) | **none** |

**3. tRNA-Gly anticodon.**

| species | tRNA-Gly | tRNAs found (tRNAscan-SE) |
|---|---|---|
| *G. jiainica* | **GCC** | 23 |
| *G. yakushimaensis* | TCC | 23 |
| *G. humatica* | TCC | 23 |
| *G. transvaalensis* | TCC | 16 (assembly incomplete) |

**4. Tree.**
- Neighbor joining on the concatenated proteins of the seven genes gives
  ((*M. servazzii*, (*A. bovina*, *A. telluris*)), (*G. yakushimaensis*, (*G. transvaalensis*, *G. humatica*)),
  *G. jiainica*).
- *G. jiainica* is the earliest-diverging of the four, as the nuclear protein identities also suggest (test 3).

**5. Branch dN/dS.**
- Method: codeml (PAML 4.10.10), yeast mitochondrial code, F3x4. Codon alignments by MAFFT L-INS-i on proteins.
- Taxa: the four *Grigorovia* species plus *M. servazzii*, *A. bovina* and *A. telluris*.
- *G. jiainica*'s four extra bases are removed. They are the frameshift steps in miniprot's alignment.
- Full tables in `genus/dnds.md`.

| gene set | *G. jiainica* ω (background) | ω = 1 rejected? | substitutions on the branch, N·dN / S·dS | 95% upper bound on ω (profile) | largest neutral share of the branch |
|---|---|---|---|---|---|
| cox1 + cob | **0.031** (0.033) | yes, p = 1.4e-7 | 63 / 228 | 0.058 | **≈ 2.5%** |
| cox1 | 0.061 (0.035) | yes, p = 0.003 | 44 / 89 | 0.19 | ≈ 16% |
| cob | 0.011 (0.031) | yes, p = 1.8e-6 | 20 / 186 | 0.027 | ≈ 0% (bound below background) |
| cox2, cox3, atp6, atp8, atp9 | 0.052 (0.037) | yes, p = 1.9e-8 | 146 / 236 | 0.091 | ≈ 5.6% |
| all seven | 0.043 (0.035) | yes, p = 1.1e-14 | 205 / 468 | 0.064 | ≈ 2.9% |

How to read the table:
- **Neutral share:** if part of the branch evolved neutrally (ω = 1) and the rest at the background ω, the branch ω is
  roughly share × 1 + (1 − share) × background. The column solves that for the 95% upper bound on ω.
- **Comparison with congeners:** their branches give ω 0.01–0.05 for every gene set (`genus/dnds.md`). *G. jiainica*
  differs from the background in no set (p ≥ 0.18), except cob, where it is *more* constrained (ω 0.011, p = 0.034).
- **Power:** the branch carries hundreds of synonymous substitutions, so ω is well estimated. The test can detect a
  neutral period longer than a few percent of the branch. It cannot see a loss of function more recent than that.

## Read against the logged predictions

- **H1 predicted** that AGG skips at GGC occur together with a GCC tRNA-Gly, and that cox1/cob are under purifying
  selection on the *G. jiainica* branch similar to the congeners'.
  - The co-occurrence holds, but with n = 1. Only *G. jiainica* has both the AGG shifts and tRNA-Gly(GCC); the three
    congeners have neither.
  - Purifying selection is strong: ω 0.031 for cox1 + cob, indistinguishable from the background.
- **H2 predicted** extra bases only in *G. jiainica* (seen), no tie to the tRNA (not seen: the one GCC tRNA is in
  the one genome with the shifts), and relaxed selection on the *G. jiainica* branch (**not seen**).
  - If cox1 and cob became dispensable, that happened within the last ~2.5% of the branch, judged by cox1 + cob. Since
    then no disabling change has appeared anywhere else in the mitogenome (test 5) or in the nuclear genes (test 3).
- **On the intron-scar idea:**
  - The cob intron 2–3 nt from 26477 is present in *G. yakushimaensis* and *G. humatica*, and absent in
    *G. transvaalensis*, which carries no extra base there.
  - The outgroups (*Arxiozyma*, *Monosporozyma*) also lack an intron at that point, and *G. jiainica* is the earliest
    branch. So *G. jiainica*'s lineage may never have had it.
  - Group I introns can move, which makes their history uncertain here. Either way, losing this intron left no extra
    base in *G. transvaalensis*.

## Added after the referee review

**Codon models and outgroups** (`genus/dnds_sensitivity.md`).
- For *cox1* + *cob* on the *G. jiainica* branch, ω = 1 is rejected under F1x4, F3x4 and F61 and without the
  outgroups.
- The largest neutral share is 1.2–5.3% across these settings. For *cox1* alone it is 5–29%, for *cob* 0–3%.
- The *cob*-only contrasts with the background vary with the model.

**GGC** (`genus/ggc_sites.md`).
- *G. jiainica*'s seven core genes have 4 GGC codons in 1,815, all at the extra-base sites.
- Every other taxon has GGT, GGA or GGG there. Gly is ancestral at all four, and the *cox1* Ser of
  *G. yakushimaensis* is derived.

**Reads** (`reads_strand.md`). The extra bases are on both strands in similar numbers, and 1,516 reads span both close
*cob* sites with both bases.
