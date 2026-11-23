# Scan report (by the stop rule)

Written 2026-11-22 ~10:00 EST.
- Rules: `PREREG.md` section 1, frozen 07:51, commit 8140435. Post-freeze descriptive additions are logged in
  the PREREG log.
- Scan: 07:51–09:52. Analysis: `code/analyze.py`. Context: `code/flag_context.py`.
- Detail files: `results/scan/analysis/`.

## Completeness

**20,343 of 20,343 genomes scored** (17,648 RefSeq non-land-plant + 2,695 fungal TPA), 0 not scored. 12 were
scored with no core gene aligned: 11 in RefSeq code 4, 1 in code 1 (`completeness.md`).

## Outcome by the stop rule

- **Strong flags outside known lineages: 1.** *Verticillium dahliae* VdLs.17 (NW_009276970), −1 at CAC, 3 of
  29 clean shifts, p = 3.5e-4, E = 0.61.
  - **Killed at the second-assembly step (pre-approved): assembly error.**
  - At all 38 clean-shift positions on the three VdLs.17 scaffolds (one assembly, one unit), the independent
    assemblies carry the base VdLs.17 lacks: 38 comparable, 38 carried.
  - Those assemblies are *V. dahliae* NC_008248 (Pantou et al. 2006) and *V. nonalfalfae* NC_029238, compared at
    121-nt windows (`second_assembly_VdLs17.tsv`).
  - Context: all 29 shifts are −1. NC_008248 and NC_029238 have 0 clean shifts. The 2011 Broad scaffolds total
    72 kb against a 27-kb genome.
- **Surviving strong flags outside known lineages: 0.**
- **Pool flags outside known lineages: 0.**

## Known lineages that flagged (15 genomes, all strong)

| lineage | genomes | flag | site sources |
|---|---|---|---|
| glass sponge | *Aphrocallistes vastus* | +1 TGG, 2/2 | both sites in Rosengarten 2008 and the record |
| *Perkinsus* | *P. marinus*, *P. olseni* | +1 AGG, 4/7 | Masuda 2010, Gornik 2022; all sites annotated |
| CTG yeasts (RefSeq) | *Huiozyma sinensis*, *Monosporozyma servazzii* | +1 CTG | Langkjaer 2003, Szabóová 2018, Wolters 2023; mostly annotated |
| **CTG yeasts (fungal TPA, unannotated)** | *Huiozyma naganishii* 3/4, *Jamesozyma lodderae* 3/3, ***Kazachstania viticola* 6/7**, *Monosporozyma aerobia* 2/3, *M. solicola* 3/4 | +1 CTG | clade documented (Wolters 2023: "several *Kazachstania*"); these records annotate nothing |
| grebes | *Podiceps cristatus*, *Tachybaptus novaehollandiae* | +1 AGT: nad3-174 + **nad1 res. 56** | nad3: Mindell 1998 et al.; **nad1: no paper in PRIOR_ART; annotated only in the *T. novaehollandiae* record, not in *P. cristatus*** |
| emydid turtles | *Graptemys ouachitensis*, *Malaclemys terrapin* | +1 AGT: nad3-174 + nad4l res. 78 | Russell & Beckenbach 2008 (CTN-AGN motif); both annotated |
| **Saccharomycetaceae, but not CTG** | ***Grigorovia jiainica*** (nom. inval.; BK073075; fungal TPA) | **+1 AGG**, 3 of 4 clean shifts (cox1, cob ×2), p = 3.5e-5, E = 0.036 | none in PRIOR_ART; record unannotated; Illumina HiSeq 2500 PE; no other assembly found |

**Caveat on *Grigorovia*:**
- It is "known" only because its family, Saccharomycetaceae, is on the known list, which was put there for the
  CTG yeasts.
- Its codon is AGG, the *Perkinsus* codon, not CTG, and nothing in our prior art describes AGG shifts in yeasts.
- By the frozen rule it does not trigger the funnel. Whether that is right is the user's decision; nothing has
  been started.
- Its own chance of flagging, given 4 random clean shifts, is 0.011. No related genome has a clean +1 shift at
  AGG.

**Pool flags:** Saccharomycetaceae and Saccharomycetales, +1 CTG, 6 of 9 sites (known).
- The shuffled-pool null gave 2, 3 and 4 family flags (seeds 1–3), above the design's original "≤ 1" target.
- Every one of those comes from a single CTG-yeast genome with ≥ 3 CTG sites (*K. viticola*, *H. naganishii*,
  *J. lodderae*, *M. solicola*) landing in a random pool.
- So in this scan the pool test only re-detects single-genome signals and adds no independent evidence.

## Chance and sensitivity

- **Chance.** 143 genome-directions have ≥ 2 clean core shifts. The expected number of chance flags under the
  rule, from the scan's own counts, is 2.5 (upper bound). Observed: 16 flags: 15 in known lineages, 1 assembly
  error.
- **Planted sensitivity (headline, unconditional; PREREG 1.12).** A genome with 2, 3 or 5 planted same-codon
  shifts flags 0.32, 0.48 or 0.72 of the time.
- **Real per-site rate**, from documented sites counted at their own codon:
  - **headline, nad3-174 excluded: 0.29** (mean of 8 lineage rates; bootstrap 95% 0.09–0.54; pooled 0.32 over 84
    sites);
  - nad3-174 alone: 0.61 (598 sites);
  - at 0.29, a genome with 2, 3 or 5 documented same-codon sites flags about **0.08, 0.20 or 0.44** of the time
    (ranges 0.01–0.29, 0.02–0.56, 0.07–0.86);
  - **out of sample** (sites that set no rule): 3 of 7 detected as clean shifts. Too few for a rate; all are
    listed in `lineage_sensitivity.md`.
- **Real-data misses** (`lineage_sensitivity.md`):
  - *Perkinsus*: 14 of 38 sites removed by the pair filter (+1 and +2 shifts close together) and 8 misplaced.
  - Birds: 130 of 527 misplaced, 47 homopolymer.
  - CTG yeasts: 4 of 11 documented sites counted.
- **Write-up limits**:
  - a negative rules out lineages only at these real-data rates;
  - the pair filter removes real sites in +1/+2 lineages;
  - the chance figures assume errors land at random codons, but assembly errors cluster, as VdLs.17 shows: strong
    without being real.

## Predictions (scored by criteria fixed before results; `predictions.md`)

| | result | count | does sensitivity explain a failure? |
|---|---|---|---|
| P1 named positives flag | PASS | 5/5 | |
| P2 another glass sponge flags at +1 TGG | **FAIL** | 0/3 | **No, it is absence or unknown, not filters.** See the sponge list below |
| P3 another CTG yeast flags at +1 CTG | PASS | 4 of 16 in the four genera; *Jamesozyma lodderae* (outside the fixed genus list) also flags | |
| P4 no bird/turtle flag involves nad3 | **FAIL** | 4 genome flags, each nad3-174 + a documented second AGT site | **No.** The prediction overlooked documented second sites |
| P5 non-known flags ≤ chance + 5 | PASS | 1 ≤ 2.5 + 5 | |
| P6 *Magnusiomyces* in the bypass tally | PASS | 65 inserts | |

**Glass sponges** (`removed_breaks_known_groups.md`): no break at TGG, and no TGG-then-A signature, in any of the
three.

| sponge | breaks written | status |
|---|---|---|
| *Oopsacas minuta* | 4: homopolymer, paired, flank | **documented absence**: its source paper (Jourda 2015, full text) observed no +1 frameshifting |
| *Caulophacus iocasicus* | 2, one cancelling pair | **absence**: no paper, nothing documented |
| *Vazella pourtalesii* | 8: paired, flank, gene end | **unknown**: its source paper is Haen et al. 2014, readable only as an abstract that names no species or positions |

*Polyrhachis dives*: 7 breaks, none in cob, so its documented-site status is unknown.

## Exploratory (descriptive, never flags)

- **Many-shift tally** (clean shifts in ≥ 3 core genes): 23 genomes.
  - Top three: two known editors (*Eukaryota* sp. BB2, Yang 2017; *Badhamia*), then VdLs.17.
  - Three tilapia (*Oreochromis*) share an identical cluster of +1 shifts in nad2.
  - 14 are marked "possible editing or noisy assembly". As worded, the mark also catches the known CTG
    frameshifters, which are all +1.
- **Paired-break tally:** 37 genome-codon rows.
  - Both *Perkinsus* appear at AGG, so the safety net works.
  - The rest are mostly AT-rich codons (ATT, TTA, AAT) in divergent arthropod genomes (mean alignment identity
    0.31–0.77). At most 3 pairs, in the gall wasp *Synergus nanlingensis* (ATT, identity 0.58).
- Other tallies: bypass-like inserts in 25 Dipodascaceae and 99 Trichomonascaceae genomes; stops read in frame;
  detours; homopolymer kills (files `tally_*.tsv`).
- **Pipeline note:** one duplicated shift event in the whole scan (*K. viticola* cob@44562, counted twice: k = 6
  is really 5 distinct sites, still strong). The planted null has 15 duplicates in 16,671, which is negligible.

## Waiting for the user (nothing started)

- whether *Grigorovia jiainica* (+1 AGG in a yeast) is treated as known by the frozen family rule or looked at;
- the grebe nad1 look (literature + grebe and flamingo comparison);
- base checks for tilapia;
- held-out known positives (*Heteranomia* first);
- the short negative write-up, which by the stop rule follows if no surviving strong flag outside known lineages
  remains.

---
**Corrections added at release review** (this report is kept as written at the time; see protocol/corrections.md):
- "mean alignment identity 0.31–0.77" for the paired-break rows is wrong. The table (paired_break_tally.md) runs from
  0.31 to 0.93.
- "CTG yeasts 4 of 11": 11 documented sites, of which 10 can be checked at their codon (4 counted). The 11th is a
  gene-level site (lineage_sensitivity.md).
- "only G. jiainica has clean AGG shifts" (tRNA paragraph, later notes) holds for core genes. *M. aerobia* has one clean
  AGG shift in a non-core ORF (fam_trna_gly.json).
