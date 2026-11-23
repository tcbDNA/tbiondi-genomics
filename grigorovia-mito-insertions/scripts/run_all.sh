#!/usr/bin/env bash
# Rerun everything, in the order it was run. Outputs go to data/work/ and results/; the curated copies in protocol/
# and candidate/ are what this run produced. Times are from a 12-core laptop.
#
#   bash scripts/setup_tools.sh && bash scripts/setup_data.sh --reads --genus && bash scripts/run_all.sh
#
# prep.py checks its taxonomy assignments against protocol/taxonomy_assignments.tsv and stops on any difference;
# add --use-archived-taxonomy to prep.py to take the published assignments instead.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$PWD/tools/local/bin:$PATH"

# 1. Inputs: genome table, CDS and intron tables, per-code query files, reference proteins (one RefSeq genome per family)
uv run python code/prep.py
uv run python code/refs.py                       # without --length-filter: the filter was tested and not kept

# 2. Prior-art audit tables used by the controls (annotated RefSeq frameshifts)
uv run python code/audit/annotation_census.py data/raw/refseq_mito/mitochondrion.1.genomic.gbff.gz refseq238
uv run python code/audit/known_cases.py data/raw/refseq_mito/mitochondrion.1.genomic.gbff.gz refseq238

# 3. Controls (about 35 min). Writes results/controls/BARS.md and thresholds.json; exit 2 = a scan gate failed
uv run python code/pick_control_sets.py
uv run python code/controls.py all || [[ $? -eq 1 ]]   # exit 1 = only non-gate bars missed (logged as limits)

# 4. The frozen scan (about 2 h) and its analysis
uv run python code/scan.py --genomes all-main --out results/scan/main
uv run python code/scan.py --genomes all-fungal --out results/scan/fungal
uv run python code/analyze.py --scan results/scan/main results/scan/fungal --out results/scan/analysis
uv run python code/flag_context.py --scan results/scan/main results/scan/fungal --analysis results/scan/analysis

# 5. The one strong flag outside known lineages: second-assembly step (pre-approved, conditional)
uv run python code/second_assembly_check.py --scan results/scan/main \
    --flagged NW_009276968.1 NW_009276969.1 NW_009276970.1 --others NC_008248.1 NC_029238.1 \
    --out results/scan/analysis/second_assembly_VdLs17.tsv

# 6. Exploratory: Grigorovia jiainica (sent to the funnel after the result; see protocol/frozen_rules.md)
uv run python code/explore/grigorovia.py a       # sites
uv run python code/explore/grigorovia.py b       # residue in relatives
uv run python code/explore/grigorovia.py c       # in-frame AGG
uv run python code/explore/grigorovia_followup.py relatives
uv run python code/explore/grigorovia_followup.py usage
uv run python code/explore/grigorovia_followup.py trnas
uv run python code/explore/grigorovia_followup.py relsites
uv run python code/explore/grigorovia.py f       # reads (needs setup_data.sh --reads)
uv run python code/explore/grigorovia_followup.py gly
uv run python code/explore/grigorovia_followup.py family
uv run python code/explore/grigorovia_followup.py heldout
# The blinded check was done by a separate agent with miniprot in an isolated folder; its task and verdict are in
# candidate/funnel/blind/ (the agent's own commands are listed in VERDICT.md).

# 7. Exploratory: six follow-up tests, bypassed or broken (needs setup_data.sh --reads --genus and setup_tools.sh)
for c in sister introns introns_gj introns_rel intron_orfs servazzii_skips mitogenome reads_strand; do
    uv run python code/explore/grigorovia_h1h2.py $c
done
for c in fetch completeness pass1 close congener reads het rates rbh; do
    uv run python code/explore/nuclear_genes.py $c
done
for c in assemble mtreads pick bait confirm genes compare dnds profile refine downloads dnds_sensitivity ggc; do
    uv run python code/explore/grigorovia_genus.py $c
done

# 8. Figure
uv run python code/explore/grigorovia_followup.py figure
uv run python code/plots/site_figure.py --data results/explore/grigorovia/figure_site.json --out plots

# 9. Added after an outside review (2026-11-23; needs setup_data.sh --reads --genus): tRNA-Gly reads, read counts with
#    duplicates removed, codon positions, tRNA-Gly in relatives, and every mitochondrial protein gene
uv run python code/explore/trna_gly_tcc.py
uv run python code/explore/site_reads_dedup.py    # uses the duplicate flags written by trna_gly_tcc.py
uv run python code/explore/codon_positions.py
uv run python code/explore/quick_checks.py trna
uv run python code/explore/quick_checks.py genes
