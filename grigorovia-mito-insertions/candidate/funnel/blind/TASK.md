# Task

This folder holds three mitochondrial genome sequences (genome_1.fna, genome_2.fna, genome_3.fna) from budding
yeasts, and reference mitochondrial proteins from two other yeasts (reference_proteins.faa; header = reference_gene).
All three genomes use NCBI genetic code table 3 (yeast mitochondrial).

Find and describe frameshifts in the protein-coding genes of each genome.

- Use miniprot (/opt/homebrew/bin/miniprot; it supports frameshift-aware protein-to-genome alignment and genetic code
  option -T) and any standard command-line tools. Do not look anything up online, and do not read or use any file
  outside this folder.
- For every frameshift you find, report: genome, gene, genome coordinate, direction (+1 = an extra base, -1 = a
  missing base), the codons around it in the original reading frame, and the amino acid the reference protein(s) have
  at the aligned position immediately after the frameshift.
- Then say, for each genome, whether its frameshifts look systematic (recurring in some pattern) or scattered, and why.
- Write everything to VERDICT.md in this folder, with the exact commands you ran.
