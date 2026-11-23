# Test 2: introns (exploratory; predictions logged in PREREG before running)

Coordinates:
- The extra base sits at a junction between two *G. jiainica* bases.
- Distance = the number of bases between that junction and an intron junction, after aligning the two genes.
  0 means the intron sits exactly where the extra base is.
- Distances in *S. cerevisiae* CDS coordinates (2a) and in *G. jiainica*'s own coordinates (2b) can differ by a few nt
  where the proteins have indels.

## 2a. Null: how often would random sites fall near an intron position?

Intron positions:
- from the RefSeq annotation of cox1 and cob in budding yeasts (Saccharomycotina; 126 species with introns);
- mapped onto *S. cerevisiae* by protein alignment;
- full table in `test2_introns.md`.

The 2,695 TPA genomes have no annotation. Their intron ends from the scan's alignments are unreliable wherever an
intron ORF continues the upstream exon (seen in the congener: 456 nt off), so they are not used here.

| intron map | distinct positions | share of coding positions within 12 nt | P(≥ 3 of 4 random sites within 12 nt) | observed nearest distances (cob 26496, cob 26477, cox1 24711, cob 26266) | observed within 12 nt | P(≥ observed) |
|---|---|---|---|---|---|---|
| *S. cerevisiae* alone | 12 (7 cox1, 5 cob) | 0.105 | **0.0042** | 7, 3, 12, 117 | 3 of 4 | 0.0042 |
| all RefSeq budding yeasts | 100 (73, 27) | 0.440 | **0.23** | 7, 3, 6, 3 | 4 of 4 | 0.038 |

Within 2 nt (the scar window):
- shares 0.022 (*S. cerevisiae*) and 0.134 (all yeasts);
- observed 0 of 4 under both maps.

**Reading.** Clustering near intron positions is clear against the *S. cerevisiae* map. It is weak against the
all-yeast map, where 44% of coding positions lie within 12 nt of some intron position. P(4 of 4) = 0.038 there, but
"3 or more", the logged test, gives 0.23. No site lies within 2 nt of any annotated yeast intron position.

## 2b. Per site

***G. jiainica* itself.** It has no cox1 or cob introns: the scan found no intron inserts in either gene, and both are
read as one continuous alignment.

**Intron-ORF remnants.**
- Method: the whole mitogenome was searched against 5,183 intron-ORF-like reference proteins (LAGLIDADG and GIY-YIG
  endonucleases, maturases, *S. cerevisiae* AI/BI ORFs, unnamed ORFs). DIAMOND --more-sensitive, frameshift-aware,
  range culling, E ≤ 1.
- Every hit inside cox1/cob is to the exon part of exon–intron fusion proteins, as expected for intact exons.
- The only hit outside core genes and rps3 is an orf146-like fragment at 11,006–11,198 (E 3e-13), about 12 kb from
  the nearest site.
- **No remnant near any site.**

| site | congener (*G. yakushimaensis*) | the six closest relatives | other Saccharomycetaceae (annotated) | other budding yeasts | *G. jiainica* lacks an intron a relative has here? |
|---|---|---|---|---|---|
| cox1 24711 | one cox1 intron, 742 nt from the site (missed by the scan's insert list; found by miniprot) | ***M. servazzii* (rank 2) and *M. solicola* (rank 4): 9 nt downstream** (equal-score splits 9–11 nt). *A. bovina*, *A. telluris*, *A. heterogenica*: no cox1/cob introns. *K. bromeliacearum*: two cox1 introns, after codons 240 and 380 of 544, far from the site | 9 nt (13 species incl. *M. servazzii*), 11 nt (1), 12 nt (6, incl. *S. cerevisiae*) | 6 nt (27, *Metschnikowia*), 9 nt (19) | yes, 9 nt away (*M. servazzii*, *M. solicola* and 12 others) |
| cob 26496 | **20 nt** downstream, the same intron as the next row | none near (*M. servazzii* and *M. solicola* have a cob intron 204 nt away) | 7 nt (6 *Saccharomyces*, the bI1 position); 12 nt (*L. kluyveri*) | 12 nt (1), 15 nt (49) | yes, 7 nt (*Saccharomyces* bI1) and 20 nt (congener) |
| cob 26477 | **2 or 3 nt** downstream (tied by sequence; 3 nt under the group I splice rule): the junction between the Gly codon (GGC in *G. jiainica*) and the Ala codon. By sequence alone, the split ties at 2 or 3 nt. 3 nt under the group I splice rule (exon ends in T, intron ends in G: `…CAT TGG GGT \| AT…ATCTTG \| GCA CTA…`); 2 nt would need the exon to end in G and the intron in T | none | 3 nt (7: *Saccharomyces* and *T. quercuum*, the bI2 position); 11 nt (6) | 3 nt (47) | **yes: the congener's intron (1,482 nt; its ORF matches the *S. cerevisiae* bI2 ORF, NP_009318, bitscore 570)** |
| cob 26266 | none | none | none within 30 nt | 6 nt (1) | no |

**Rfam check on the intron type.**
- The congener intron gave no significant hit to the Rfam group I (RF00028) or group II (RF00029) models: best E 0.32
  and 0.24.
- Divergent mitochondrial introns are often missed by these models, so this does not decide the type.
- The homology of its ORF to the *S. cerevisiae* bI2 ORF, and its position, match bI2, which is described as group I.

***M. servazzii*'s own cob +1 sites (CTG·T·GG).** They map 165–204 nt upstream of the *G. jiainica* cob sites, so they
are not shared positions.

## Read against the logged predictions

- **2a.** H2 (scar) predicted P < 0.05.
  - Met under the *S. cerevisiae* map (0.004).
  - Not met under the all-yeast map for "≥ 3 of 4" (0.23); 4 of 4 gives 0.038.
- **2b.** H2 (scar) predicted sites at or within 2 nt of junctions of introns *G. jiainica* lost.
  - **Not met as written.** The nearest case is cob 26477, 3 nt from the congener's bI2-type intron (2 nt only if
    the group I rule is set aside).
  - The other sites are 7–20 nt from such introns (cob 26496, cox1) or near none (cob 26266).
  - No intron-ORF remnant lies near any site.
- **Overall.** Three of four sites lie near intron positions that the genus or its relatives use. *G. jiainica* lacks
  all cox1/cob introns, including the one its congener has 3 nt from a site. But the sites are not at the junctions,
  and there are no remnants. This weakly supports a connection to intron loss, without the precise scar H2 predicted.

## Correction after a cross-check (same session)

**What changed.**
- Exon structures of cox1/cob in the congener and the six closest relatives were rebuilt with miniprot 0.18
  (*S. cerevisiae* proteins, code 3; records rotated where a gene straddles their ends). Details:
  `test2b_relatives_miniprot.md`.
- The scan's insert list had missed two introns: the congener's cox1 intron (742 nt from the site), and the first
  cox1 intron of *M. solicola*, which lies **9 nt** from the cox1 site, as *M. servazzii*'s does.
- The table above was corrected to match.

**What did not change.**
- No relative has an intron within 2 nt of a site.
- At cob 26477, miniprot's exon ends for the congener agree with the spliced alignment: 2–3 nt, tied.
