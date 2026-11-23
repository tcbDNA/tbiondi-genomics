# Gate 1: Data inventory and nuclear–mitochondrial overlap

Date: 2026-11-23. Status: **inventory complete; scope decision needed.** Firewall: only metadata, sequences and annotation fields were read. No `trnY` extractor, UAA caller or `mtRF1a` caller was run on any genome in this gate (the Gate 0 checks used the named controls only).

> **Corrections made 2026-11-23 after the first version was presented.** (a) The INSDC-only query required "complete genome" in the title and excluded chordates. It misses 3 of 6 known INSDC-only controls checked (`FN313571`, `OZ281197`, `FN908482`), so **7,629 is an undercount** and the query is replaced in Gate 2 prep. (b) The pairing percentages quoted in the summary message (48/59/50%) count RNA-seq runs; **assembly-only pairing is Mollusca 385/1,078 (36%), Nematoda 164/292 (56%), Arthropoda 3,909/9,189 (43%)**.

## 1. Mitogenome inventory
Sources: RefSeq mitochondrion release 238 (18,220 records; 16,075 Metazoa, 2,145 excluded as non-Metazoa, logged in `results/audit/gate1_refseq_excluded.tsv`) plus 6,467 NCBI "genome assembly" organelle records (2,640 DToL, 353 ATLASea, 194 other Sanger/ToL-tagged, 3,280 other). 57 RefSeq records have no inline sequence; their FASTA was fetched separately (6,524 FASTA downloaded, all checksummed).
Organelle sequences that duplicate a RefSeq source accession: 75. Species already represented by another record: 2,272 (logged, not deleted). Open-nomenclature taxa ("sp.", "cf.", "aff.") are kept separate, because taxid 3031986 "Pseudotanais sp." covers at least two species (Zhao: sp. a and sp. b). **20,195 species-level representatives remain (11,925 non-chordate).** The unique-taxid rule is a simplification; a species-level set is one sequence per species, so within-species variation is ignored.

| Phylum | Mitogenome species | circular | linear | nuclear assembly (species) | nuclear assembly (genus) | annotated assembly | chromosome-level | OrthoDB | RNA-seq (ENA) | any nuclear source (species) | genus-only |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Arthropoda | 9189 | 9058 | 131 | 3909 | 5532 | 317 | 3515 | 735 | 2519 | 4572 | 1360 |
| Chordata | 8270 | 8020 | 250 | 3187 | 5810 | 1062 | 1641 | 840 | not queried | ≥3187 | 2622 |
| Mollusca | 1078 | 954 | 124 | 385 | 612 | 48 | 308 | 34 | 394 | 517 | 154 |
| Cnidaria | 438 | 386 | 52 | 194 | 275 | 25 | 103 | 15 | 193 | 240 | 64 |
| Nematoda | 292 | 282 | 10 | 164 | 215 | 58 | 113 | 50 | 107 | 171 | 49 |
| Annelida | 253 | 214 | 39 | 113 | 157 | 4 | 101 | 2 | 119 | 152 | 31 |
| Platyhelminthes | 205 | 193 | 12 | 60 | 104 | 37 | 35 | 35 | 62 | 70 | 43 |
| Echinodermata | 176 | 156 | 20 | 77 | 107 | 12 | 61 | 6 | 74 | 93 | 25 |
| Porifera | 159 | 140 | 19 | 109 | 135 | 7 | 86 | 3 | 106 | 121 | 18 |
| Nemertea | 25 | 24 | 1 | 7 | 10 | 2 | 5 | 1 | 8 | 9 | 2 |
| Bryozoa | 18 | 16 | 2 | 9 | 12 | 2 | 9 | 1 | 12 | 13 | 2 |
| Acanthocephala | 15 | 15 | 0 | 2 | 5 | 1 | 0 | 1 | 3 | 4 | 3 |
| Xenacoelomorpha | 14 | 12 | 2 | 5 | 10 | 2 | 4 | 0 | 9 | 9 | 4 |
| Ctenophora | 8 | 8 | 0 | 7 | 8 | 1 | 3 | 0 | 7 | 7 | 1 |
| Brachiopoda | 7 | 6 | 1 | 2 | 3 | 1 | 2 | 1 | 3 | 4 | 1 |
| Nematomorpha | 7 | 7 | 0 | 5 | 6 | 1 | 5 | 0 | 3 | 6 | 1 |
| Hemichordata | 7 | 6 | 1 | 2 | 4 | 2 | 2 | 1 | 4 | 4 | 2 |
| Chaetognatha | 5 | 5 | 0 | 2 | 2 | 0 | 2 | 0 | 3 | 3 | 0 |
| (unmapped) | 5 | 5 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 3 |
| Onychophora | 4 | 4 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 1 | 1 |
| Kinorhyncha | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 |
| Placozoa | 3 | 3 | 0 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 0 |
| Tardigrada | 3 | 3 | 0 | 2 | 2 | 1 | 0 | 1 | 2 | 2 | 0 |
| Priapulida | 2 | 1 | 1 | 2 | 2 | 1 | 2 | 1 | 2 | 2 | 0 |
| Gnathostomulida | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 |
| Rotifera | 2 | 2 | 0 | 1 | 2 | 1 | 0 | 1 | 2 | 2 | 0 |
| Entoprocta | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Phoronida | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 |
| Orthonectida | 1 | 1 | 0 | 1 | 1 | 1 | 0 | 1 | 1 | 1 | 0 |
| Gastrotricha | 1 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 1 | 0 |

Total species-level representatives: 20195

Notes on the table:
- "any nuclear source" = assembly, OrthoDB, or ≥1 RNA-seq run for the species. RNA-seq was queried for all invertebrates (ENA, `library_strategy="RNA-Seq"`, 427,466 runs) but not for chordates; Chordata RNA-seq is therefore "not queried", not zero.
- Completeness flags on species representatives: RefSeq 15,936 (15,676 titled "complete genome"; 15,867 with ≥12 CDS annotated; 48 with 0 CDS). All organelle representatives (ATLASea 278, DToL 2,238, Sanger-other 97, other 1,646) have **no CDS annotation in NCBI**, apart from 302 of the 1,646 "other". They must be annotated from sequence at Gate 2 (homology search), which is the main compute cost.
- Linear ≠ incomplete: Medusozoan mtDNA is biologically linear (Cnidaria: 52 linear representatives in this set are expected to be partly biological). 182 representatives are >30 kb and 7 are <10 kb; both get QC flags.

## 2. Control coverage in the sweep set
| Control | In the RefSeq + organelle set? |
|---|---|
| *Anomia ephippium* OZ576062, *Heteranomia squamula* OZ181722, *B. xylophilus* NC_023208 | yes |
| *Cephalodiscus* MF374802, *Radopholus* FN313571, *Pseudotanais* OZ281197–OZ281200, *Rhabdopleura* FN908482, *Arctotanais*, *Sinelobus* | **no — INSDC only** |

Positive controls will be added by accession at Gate 2, but the **sweep set as defined (RefSeq + DToL) excludes the known reassigned lineages themselves** and roughly 7,629 invertebrate species names (counted by NCBI species name; may include open-nomenclature names) that have a complete mitogenome only in INSDC (29,812 non-chordate records; `results/insdc_invert_summaries.tsv`). Example: Arthropoda 5,714 new species, Mollusca 590, Cnidaria 410, Annelida 277, Platyhelminthes 239, Echinodermata 140, Nematoda 103.
**Recommendation:** extend the sweep to INSDC-only species (cheapest form: one record per species not in RefSeq). Titles say "complete genome" and are unreviewed, so quality flags must be strict. Until approved, this is documented, not done.

## 3. Nuclear pairing for `mtRF1a`
Species-level pairing counting any source, RNA-seq included (assembly-only figures are in the correction note above) across invertebrate phyla: Mollusca 517/1,078 (48%), Nematoda 171/292 (59%), Arthropoda 4,572/9,189 (50%), Platyhelminthes 70/205, Echinodermata 93/176, Cnidaria 240/438, Annelida 152/253; **Hemichordata 4/7, Tanaidacea 0 nuclear assemblies, Pterobranchia 0 nuclear assemblies.**
For the reassigned lineages specifically (Gate 0 §6): *Anomia* (2 assemblies), *Enigmonia* (1), *Radopholus* (3) are usable from assemblies. *Cephalodiscus*, *Rhabdopleura*, *Placuna*, *Plicatula*, *Pseudotanais*, other tanaidaceans and *Pododesmus* depend on RNA-seq or raw reads, i.e. assembly from reads (Zhao did this too).
Implication: a pan-metazoan `mtRF1a` genotyping is feasible for molluscs, nematodes and arthropods, but **the hypothesis cannot be tested genome-wide in the clades where the claim matters most** (Hemichordata, tanaids) without read-based assembly. The planned Gate 4 design — flagged lineages plus ≥2 sister controls plus a background panel — fits the available data; a full census does not.

## 4. Scope recommendation (for your decision)
**`FULL_TWO_GENOME`, with the `mtRF1a` half restricted to (a) reassigned lineages and matched sisters and (b) a background panel of species with annotated or chromosome-level assemblies.** Pan-metazoan statements about `mtRF1a` would not be made. If you prefer, I can fall back to `MITO_ONLY_FALLBACK` (then the `mtRF1a` half is reported as untested).

## 5. Resources used in Gates 0–1
About 8.7 GB on disk (`data/` 5.3 GB, mostly unpacked Pfam and taxonomy; `env/` 3.4 GB, incl. a package cache that can be deleted); measured downloads are listed in Gate 0 §8. Nothing exceeded 20 GB or 1 hour. No large download is pending approval.

## 6. Open questions
1. Scope (section 4).
2. Approve Lin's Figshare (1.04 GB) and/or Zhao Additional file 4? Needed for *Enigmonia*, *Plicatula*, *Placuna vitream*, *Pododesmus* and an uncorrected *A. chinensis*/*P. placenta* control set (Gate 0 §5).
3. Approve the INSDC extension of the sweep (section 2)?
4. Optional papers (Gate 0): Jakiel 2025, Lin supplement, Kakui & Kano 2021, Zhao Additional files 4 and 9.
