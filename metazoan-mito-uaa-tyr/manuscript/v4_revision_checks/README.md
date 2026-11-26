# Checks on the v4 revision (2026-11-25)

These files check the numbers in `../metazoan_mitochondrial_uaa_to_tyr_preprint_v4_revision.tex` and in its redrawn figures (`../figures_v4/`) against the project outputs and against the v3 manuscript. The v3 files are unchanged.

| File | Contents |
|---|---|
| `number_trace_v4.tsv` | The number trace from `src/preprint_trace/` run on the revision and on the labels drawn in `../figures_v4/`. It has 2,406 rows: 1,743 matched to a source, 661 labels (parts of names such as T216) and 2 untraced. The two untraced rows are the "GPT-6.2" model names in the AI disclosure, as in v3. No row is a mismatch. |
| `spelled_out_counts.tsv` | 75 values the trace's hand-set patterns do not reach, because those patterns were written for v3 sentences. Most are counts written as words ("two fish", "seven genes"). Each was checked by hand against the source named in the row: 67 are counts and all are correct, 7 are wording such as "two ways", and 1 is a figure label ("base 33"). |
| `moved_or_omitted.tsv` | v3 sentences whose numbers no longer appear together in the revision, and two values that were removed. Each row is marked `kept` (the sentence was only split), `moved` (to a table, the Methods or a figure) or `omitted` (with where the value can still be found). |

Three values in the revision do not appear in v3:
- "203", the residue number of Q203, in the heading of the Figure 3 caption.
- "97.6%" and "16,852 bp": the *M. kubotai* read fish (`SRR39401262`) has an assembled mitochondrial contig of 16,852 bp that is 97.6% identical to the screened record (`NC_037360.1`, 18,168 bp). Sources: `reports/REPORT.md` line 246 and `results/posthoc2/conspecific_specimen_notes.txt` line 84; that line also gives 99.76% identity to `PV357792`. The record length is in `results/posthoc2/secondcopy/highcopy_summary.txt` line 31.

As a control, the same scripts were run on the v3 file. That run reproduced `results/posthoc2/preprint_number_trace.tsv` byte for byte.

To reproduce, run from the project folder:

```
env/prefix/bin/python src/figures.py manuscript/figures_v4 --plain
env/prefix/bin/python src/preprint_trace/extract.py manuscript/metazoan_mitochondrial_uaa_to_tyr_preprint_v4_revision.tex
env/prefix/bin/python src/preprint_trace/run.py manuscript/v4_revision_checks/number_trace_v4.tsv
```

For the figure rows, `run.py` was pointed at `manuscript/figures_v4/` instead of `reports/figures/`; in the committed script that path is on the line that opens `reports/figures/{f}.svg`. `extract.py` overwrites `results/work/preprint_trace/doc.json`, which is not tracked. The reference list is unchanged from v3, so `refs.py` was not run again.

`src/figures.py` now takes an output folder and `--plain`. Run without them it draws the v3 figures into `reports/figures/`: a run into a scratch folder matched the committed SVGs apart from the random element IDs that matplotlib writes on every run.
