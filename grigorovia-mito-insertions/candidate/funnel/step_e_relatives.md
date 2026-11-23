# Step e: do the six closest relatives carry the extra base? (+-60 nt windows, local alignment; comparable = >= 90 aligned nt at >= 70% identity; automatic call, check gaps placed in repeats by eye)

| site (1-based) | Kazachstania bromeliacearum | Monosporozyma servazzii | Arxiozyma bovina | Monosporozyma solicola | Arxiozyma telluris | Arxiozyma heterogenica |
|---|---|---|---|---|---|---|
| 24711 | lacks (0.85) | n.c. (0.87) | lacks (0.83) | n.c. (0.86) | lacks (0.80) | lacks (0.83) |
| 26496 | lacks (0.89) | carries (0.90) | lacks (0.91) | lacks (0.91) | lacks (0.92) | lacks (0.88) |
| 26477 | lacks (0.89) | lacks (0.88) | lacks (0.91) | lacks (0.88) | lacks (0.91) | lacks (0.88) |
| 26266 | lacks (0.88) | lacks (0.85) | lacks (0.88) | lacks (0.85) | lacks (0.88) | lacks (0.87) |

## Checked by eye (coding strand)

The automatic call looks for a gap in the relative within 3 nt of the extra base. In a TA repeat the aligner can put
the gap 4 nt away, which the test then reads as "carries". The one such call in this table was re-aligned on the
coding strand and inspected:

- **M. servazzii at 26496 (AGG site): lacks the extra base.** *G. jiainica* `GTATAT A GGCCAAATGTCACATTG A TGG`
  aligns to *M. servazzii* `GTATAT - GGGCAAATGTCACATTG A - GG`. M. servazzii lacks both the AGG-site A and the TGG-site
  T.
- M. servazzii's own documented CTG sites in cob (record positions 8588 and 8610) lie about 150 nt from this region.
  They are not the same sites.

Corrected summary: none of the six relatives carries the extra base at any of the four positions. Of 24 comparisons,
22 lack it and 2 are not comparable.
