# Test 4: congener phenotypes (exploratory; literature search by a fresh agent, 2026-11-22)

Growth or assimilation on non-fermentable carbon. Quotes are short; "depth" says how each source was read.

| species | ethanol | glycerol | lactate | succinate | citrate | other | source | quote or table cell (length) | depth |
|---|---|---|---|---|---|---|---|---|---|
| *G. yakushimaensis* | − | − | − | − | − | methanol −, gluconate − | Mikata, Ueda-Nishimura & Hisatomi 2001, *IJSEM* 51:2189 (original description, as *Saccharomyces yakushimaensis*) | "Glucose, galactose and trehalose (slow) are assimilated, but not … ethanol, methanol, glycerol … lactic acid, succinic acid, citric acid" (≈ 22 words, with ellipses); Table 3: Eth −, Gly −, LA −, SA −, GA − | full-text PDF; Table 3 viewed as an image |
| | − | − | − | − | − | erythritol, mannitol, inositol − | Opulente et al. 2018, *BMC Biol*, Table S10 (curated from *The Yeasts*, 2011) | each cell "-" (row label misspelt "K. yakushimaenisis") | spreadsheet cells |
| *G. humatica* | − | − | − | − | − | methanol −, gluconate − | Mikata et al. 2001 (original description, as *S. humaticus*) | "Glucose and galactose are assimilated, but not … ethanol, methanol, glycerol … lactic acid, succinic acid, citric acid" (≈ 20 words); Table 3 all − | full text and table |
| | − | − | − | − | − | | Opulente 2018 Table S10 | each cell "-" (row label "K. humantica") | spreadsheet cells |
| *G. transvaalensis* | v | v | − | − | | gluconate + | Mikata et al. 2001, Table 3 (their tests on type strain IFO 1625) | "Eth v, Gly v, LA −, SA −, GA +" (table row) | table viewed as an image |
| | − | − | − | − | − | methanol − | Opulente 2018 Table S10 | each cell "-" | spreadsheet cells |
| | | **+** | 0 | | 0 | | Y1000+ phenotypes (Opulente et al. 2024), type strain CBS 2186 | glycerol growth rate 0.025 (glucose 0.059; no carbon 0), fitted in 3 replicates plus a repeat | data files |
| | type culture does not respire | | | | | petite test | Fekete et al. 2007, *FEMS Yeast Res* 7:1237, Table 1 (NRRL Y-17245 / CBS 2186) | "Original collection strain already does not grow on nonfermentable substrate" (11 words) | full text via a web-fetch tool; quote not checked against the PDF |
| *G. jiainica* (for comparison) | − | − | − | − | − | | Lee et al. 2008; Opulente 2018 Table S10 | see `search_log/grigorovia_literature.md` | |

**Not found or not read.**
- No petite or respiration data for *G. yakushimaensis* or *G. humatica*.
- The *Kazachstania* chapter of *The Yeasts* (2011): access refused. Table S10 says it is curated from it.
- Van der Walt's 1956 description of *S. transvaalensis*: title only.

**Reading against the logged prediction.**
- *G. yakushimaensis* and *G. humatica* assimilate no non-fermentable carbon either, just like *G. jiainica*.
- Yet *G. yakushimaensis* has intact cox1 and cob without the extra bases (OR703946; held-out run, 0 AGG shifts). So
  non-growth on these carbon sources in standard tests does not by itself go with mitochondrial gene decay in this
  genus. **This weakens the inference "it cannot grow on glycerol, so its cox1/cob are dispensable", which H2 relied
  on.** It does not rule H2 out: the congener could keep the genes for a reason other than respiration, or its loss
  could be more recent.
- *G. transvaalensis* is uncertain: glycerol variable or positive in two sources, negative in the curated table, and
  its type culture is respiration-deficient in one stock.
- No congener is documented to respire while *G. jiainica* does not. So the result that would have favoured a
  *G. jiainica*-specific loss was not seen.

## The two controls used in test 3 (same literature search; quotes or cells as reported)

| species | carbon source or test | value | source | depth |
|---|---|---|---|---|
| *M. servazzii* | glycerol | + | Mikata et al. 2001, Table 3 | full text |
| | ethanol | − | Mikata et al. 2001, Table 3 | full text |
| | glycerol / ethanol / succinate | + / v / v | Opulente 2018 Table S10 | spreadsheet cells |
| | glycerol (liquid growth rate) | 0 | Y1000+ (Opulente et al. 2024), CBS 4311 | data files |
| | petite test | "+" with no asterisk for the type strain CBS 4311, which implies growth on the glycerol/ethanol test medium (our inference) | Fekete et al. 2007, Table 1 | full text via a web-fetch tool |
| *K. africana* | DL-lactate / glycerol | + / v | Opulente 2018 Table S10 | spreadsheet cells |
| | petite test | petite-negative at 28 °C (depends on respiration) | Fekete et al. 2007 | full text via a web-fetch tool |

So each control grows on a non-fermentable carbon source in at least one record. *M. servazzii*'s records disagree
(liquid glycerol growth 0 in Y1000+).
