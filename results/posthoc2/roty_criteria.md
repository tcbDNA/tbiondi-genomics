# Check 6: criteria written before mapping or assembly (post-hoc, cold review 2, closing question 3)

Reads: SRR26324223 (paired WGS). References: CM003310.1 (mitochondrion) and GCA_001026735.1 (nuclear assembly; any contig that is the mitochondrion is removed from the nuclear set and listed). What was already known when this was written: the locked call (22,645-ish arm, 3/5 pairs at model columns, score 12.8), ARAGORN's 22643-22704 (62 nt), the 5/5 contiguous count, and exact-k-mer counts at internal UAA sites in this run. The run's reads had not been mapped to the trnY locus or assembled.

| # | Question | Measure | Pass | Fail |
|---|---|---|---|---|
| C1 | Is the locus mitochondrial? | primary alignments covering N33 in competitive mapping; nuclear places with a GUA arm or >=80% identity to the tRNA over >=50 nt | >=95% of them on CM003310.1 with MAPQ >=20, and no nuclear copy with depth >=10% of the mitochondrial depth | otherwise ("nuclear copy" or "ambiguous") |
| C2 | Is the allele U33? | base counts at N33 and at every position of the 62-nt ARAGORN interval, base quality >=20 | >=100 reads, U >=99% at N33, and no other position with a second base >=5% | a second base >=5% at N33 = "conflicting allele" |
| C3 | Does an independent local assembly recover it? | SPAdes contig from the recruited read pairs; also from all reads that map to CM003310.1 | one contig holds the whole locus with >=200 bp flanks on both sides that match CM003310.1 at >=99% | no such contig, or a different tRNA sequence (differences listed) |
| C4 | Does either boundary map make a tRNA? | on the assembled contig, for the Infernal interval and the ARAGORN interval separately: acceptor stem pairs (first 7 against the 7 before the discriminator base, G.U allowed), anticodon arm counted beside the loop, loop of 7 nt with GUA at 34-36, overlap with the neighbouring genes (annotated or profile) <=10 nt (PREREG 2.3) | "supported" if acceptor >=5 of 7 and anticodon stem >=4 of 5 and overlap <=10 nt | otherwise "not supported" |
| C5 | What does the registered rule say on the assembled locus? | locked extractor `trny_extract.call` and `gate3_classify.n33` logic (score, callers, stem at model columns) run on the contig | reported as resolved / stem_audit_review / unresolved, whatever it is | - |

Reading, fixed in advance:
- C1 or C2 fails, or C3 gives a different tRNA: the proposed U33 genotype is overturned.
- C1-C3 pass, exactly one map supported (C4): that map is the better boundary; say whether the registered rule (C5) still blocks it. If it does, the zero-new-clade result depends on the locked alignment implementation. The locked result is not changed.
- C1-C3 pass, both maps supported or neither: the boundary stays unresolved, and the candidate stays "unconfirmed" under the registered rules.
- No outcome changes PREREG, Gate 3 or the status "unconfirmed single-genome candidate".
