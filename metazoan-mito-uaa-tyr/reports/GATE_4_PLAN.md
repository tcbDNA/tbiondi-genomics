# Gate 4 plan (printed before downloads, 2026-11-24)

Scope set by the user: Tier A confirmatory per PREREG §2.4, plus *Rotylenchulus reniformis* and tylenchid Stop sisters (sensitivity label, because its trnY is blocked by the stem audit); Tier B exploratory. Source order per §2.4: annotated proteins or OrthoDB, TSA, assemblies (miniprot), raw reads; otherwise untested. Guard: `src/mtrf_guard.py`, `GATE3_UNLOCK=1`. Numbering: human NP_061914.3 positions 203 and 216. The pre-Gate-4 checks are in `reports/GATE_3B_EXPLORATORY.md` and `results/gate3/`.

| Tier | Taxon | Source | Download (est.) |
|---|---|---|---|
| A | reassigned lineages, sisters, background of Gate 2 | Gate 2 files on disk | 0 |
| A | *Plicatula plicata* | SRR3350473 (RNA-seq, PRJNA316843) | 3.95 GB |
| A (sensitivity) | *Rotylenchulus reniformis* | TSA GGVO, GGVQ, GGVP, GGVR, GGVV; then GCA_001026735.1; then ERR744555 if needed | 157 MB + 88 MB (+2.1 GB) |
| A (sensitivity) | *Pratylenchus vulnus* GCA_054131405.1 (Gate 3 Stop); *Meloidogyne hapla* GCA_051171035.1 proteins (Stop, exploratory set); *Hoplolaimus columbus* GCA_016433075.1; *Rotylenchus robustus* GCA_034700465.1 | proteins or assemblies | ~90 MB |
| B | *Albula glossodonta* proteins GCA_019788955.1 and assembly GCA_056783665.1 (on disk); sisters *Megalops cyprinoides* GCF_013368585.1, *Anguilla anguilla* GCF_013347855.1 | proteins | ~40 MB |
| B | *Lycenchelys verrillii* GCA_056747025.1 (on disk); sisters *Lycodopsis pacificus* GCA_028022725.1, *Zoarces viviparus* GCA_040110945.1 | assembly; proteins | ~25 MB |
| B | *Craniella zetlandica* GCA_978397455.1, *C. infrequens* GCA_974939545.1 (on disk); sister *Cinachyrella apion* GCA_982936665.1 | assemblies | 82 MB |
| B | *Microdevario*, Serrivomeridae | no nuclear assembly | untested |

Skipped by instruction: *Bufotes*, *Globodera*, *Coloconger*, *Ophichthus*, *Xyrias*, *Toxocara*, the crane. Estimated total 4.5 GB (6.6 GB with the RNA-seq fallback), peak about 15 GB, floor 30 GB free.
