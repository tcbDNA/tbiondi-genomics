# GGC codons in G. jiainica's core genes, and the residues at the four sites across the genus (post hoc)

- G. jiainica, seven core genes (cox1, cob, cox2, cox3, atp6, atp8, atp9; extra bases removed): 1,815 codons, **4 GGC**; GC at third codon positions 5.4%.
- Gly codons: GGT 67, GGA 38, GGG 0, GGC 4.

## The codon right after each extra base (scan placement) and the next codon, in every taxon (codon alignments, `codon_alignments/`)

At cob 26266 the scan places the extra base as the A of ACT, so GGC is the second codon; placing it as the A of AGG instead makes GGC the first, with the same protein. Gaps at cob 26477 in G. yakushimaensis and G. humatica are codons lost at their intron junction in the CDS extraction (removed by codeml's cleandata).

| site | Gjiainica | Gyakushima | Ghumatica | Gtransvaal | Mservazzii | Abovina | Atelluris |
|---|---|---|---|---|---|---|---|
| cox1 24711 | GGC (G), then TTT (F) | AGT (S), then TTT (F) | GGT (G), then TTT (F) | GGT (G), then TTT (F) | GGT (G), then TTT (F) | GGA (G), then TTT (F) | GGT (G), then TTT (F) |
| cob 26266 | CTA (T), then GGC (G) | CTT (T), then GGT (G) | CTT (T), then GGT (G) | CTA (T), then GGT (G) | ACA (T), then GGA (G) | ACA (T), then GGT (G) | ACA (T), then GGT (G) |
| cob 26477 | GGC (G), then GCA (A) | --- (gap), then --- (gap) | --- (gap), then --- (gap) | GGT (G), then GCT (A) | GGA (G), then GCA (A) | GGA (G), then GCA (A) | GGA (G), then GCA (A) |
| cob 26496 | GGC (G), then CAA (Q) | GGG (G), then CAA (Q) | GGT (G), then ATA (M) | GGA (G), then CAA (Q) | GGG (G), then CAA (Q) | GGA (G), then CAA (Q) | GGA (G), then CAA (Q) |

GGC codons in G. jiainica's seven core genes: 4; of them at the extra-base sites (first or second codon after the base): 4.
