# Metazoan mitochondrial UAA -> Tyr check
Gated, pre-registered study (see `reports/`). Gates 0 and 1 done; `PREREG.md` not yet written, so no caller has been run on any non-control genome.
Environment: `env/conda-osx-arm64.lock` (arm64), `env/conda-hmmer332-osx-64.lock` (Rosetta; Codetta's HMMER 3.3.2). Codetta commit `863359e` in `env/tools/codetta`.
Data: `data/` (gitignored; checksums in `data/manifests/`). Full texts in `data/sources/fulltext/` are gitignored (copyright).
Scripts: `src/gate0_*.py` (Gate 0), `src/gate1_*.py` (Gate 1). All download steps are resumable and log to `results/`.
