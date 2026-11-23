# NCBI taxonomy rows used by the analysis

Cut by `code/release/make_taxonomy_extract.py` from the NCBI taxonomy dump downloaded on 2026-11-21 (`taxdump.tar.gz`, md5 6d18d6c9d60bbc20759869f6b3592eff): 39029 nodes (the lineages of 20924 query taxa: the scanned genomes and the BUSCO saccharomycetes_odb12 reference species), 7 merged-id rows and 39029 scientific names, in the dump's own format. `prep.py`, `flag_context.py`, `explore/nuclear_genes.py` and `explore/grigorovia_h1h2.py` read these files when they are present, so a rerun reproduces the published taxonomy assignments.
