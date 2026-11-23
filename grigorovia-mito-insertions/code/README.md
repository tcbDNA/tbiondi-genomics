# code/

Run with `uv run python code/<script>.py`. `../scripts/run_all.sh` lists every command in order. Inputs come from
`../scripts/setup_data.sh`, and tools from `../scripts/setup_tools.sh`.

| file | what it does |
|---|---|
| `prep.py` | builds `data/work/` from the raw records: genome table (set, taxid, code, family and order pools), CDS and intron tables, per-code query FASTA, reference proteins |
| `refs.py` | final reference set and DIAMOND database (NC_079576 excluded; `--length-filter` was tested and not kept) |
| `fscore.py` | the core shared by real, planted and null genomes: DIAMOND call, BTOP walk, frameshift events, the placement window (`trigger`), homopolymer test, inferred inserts |
| `scan.py` | per-genome processing (`process_genome`) and the scan CLI |
| `flags.py` | the frozen single-genome test (window rule, planted-null rates, alpha), expected chance |
| `controls.py` | every control stage (planted null, null genomes, named positives, annotated shifts, introns, same-codon trials) and `bars` |
| `analyze.py` | genome and pool flags, the strong rule, shuffled pools, exploratory tallies, coverage |
| `units.py` | assembly/strain units for the strong rule |
| `flag_context.py` | post-freeze descriptive context: per-flag context, completeness, lineage sensitivity, tallies, predictions |
| `second_assembly_check.py` | the second-assembly step (local alignment of windows around each shift) |
| `pick_control_sets.py` | fixes the control genome sets (seeds 7 and 11) |
| `plots/site_figure.py` | the static site figure |
| `explore/grigorovia.py` | *Grigorovia* funnel steps a, b, c, f |
| `explore/grigorovia_followup.py` | relatives, codon use, tRNAs, tRNA-Gly, family-wide tRNA search, relatives' sites, held-out relative, figure data |
| `explore/grigorovia_h1h2.py` | follow-up tests 1, 2 and 5: congener codons, introns, intron ORFs, the rest of the mitogenome, reads by strand |
| `explore/nuclear_genes.py` | follow-up test 3: nuclear genes, completeness, read checks, heterozygosity, orthology (uses `nuclear_calls.py`, `read_haplotypes.py`) |
| `explore/grigorovia_genus.py` | follow-up test 6: congener mitogenome assembly and confirmation, comparison, dN/dS and its sensitivity, GGC pattern, downloads |
| `release/export_prereg.py` | writes `protocol/frozen_rules.md` from the working repository's PREREG.md and git log |
| `audit/` | prior-art audit: annotation census, known RefSeq frameshift junctions, NCBI lookups, fungal TPA download, forward citations |
| `lit_fetch.py` | Europe PMC search helper used for the prior-art audit |

Not in code: the blinded re-check was run by a separate agent in an isolated folder with miniprot. Its task and its
commands are in `../candidate/funnel/blind/`.
