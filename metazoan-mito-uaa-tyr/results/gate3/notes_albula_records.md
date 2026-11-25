# Albula records: are the two "A. glossodonta" entries independent? (checked 2026-11-24, at the user's request)

| Record | Specimen / source | Notes |
|---|---|---|
| AP002973.1 (named_extra in the sweep) | voucher CBM-ZF 10501 (Chiba Museum, Japan); Inoue et al. 2004; submitted 2000; Sanger-era | named "A. glossodonta" in the record |
| NC_005800.1 | RefSeq; COMMENT: "The reference sequence was derived from AP002973"; 100.00% identical over 16,255 nt | **not in the sweep** (absent from RefSeq release 238 file; accession search returns nothing) |
| CM165726.1 (species representative in the sweep) | isolate fAlbGlo1, Kaunakakai, Hawaii, 2025-07-14, blood; BioSample SAMN55933885; Vertebrate Genomes Project PRJNA1427880; PacBio HiFi + Hi-C | named "A. glossodonta" |
| NC_082992.1 (= OR546140.1) | A. vulpes, voucher USNM:FISH:413592, Curacao; SRR22396674; panel genome | panel role |

Conclusions
- **Not duplicates.** AP002973 and CM165726 are different specimens, labs, years and technologies. NC_005800 is AP002973 (derived and identical), so it would be a duplicate, but it is not in the sweep.
- **They are not the same species by DNA.** Full-mitogenome identity AP002973 vs CM165726 is 87.0% (about 13% divergence); CM165726 vs A. vulpes NC_082992 is 92.9%. One of the two "A. glossodonta" labels is therefore doubtful, and AP002973 is the more divergent one.
- AP002973 counts as **one record of an Albula lineage with an uncertain species label** (the user cites Pfeiler et al. 2006: possibly A. argentea/forsteri). I could not retrieve that paper; a web search found only general statements that bonefish misidentification is known (Pfeiler 2008, citing Pfeiler et al. 2006). The sequence agrees with "not A. glossodonta" but does not identify the species.
- All three Albula records (A. vulpes NC_082992 [panel], CM165726, AP002973) and Pterothrissus gissu have C33 with loop CCGUAGC, UAA = Stop. The trnY of CM165726 and A. vulpes are identical (71 nt); AP002973's differs at a few positions (loop unchanged).
- Earlier statements to correct: this is **not** within-species replication of A. glossodonta. It is C33 in several Albuliformes lineages, from independent specimens. The T6 reading ("lineage-level") does not depend on AP002973.
