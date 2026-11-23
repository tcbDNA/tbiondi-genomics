#!/bin/bash
# Full pipeline for VERDICT.md. Run from anywhere: bash work/reproduce.sh
set -euo pipefail
D=<blind folder>
cd "$D"
MP=/opt/homebrew/bin/miniprot
S=work/scripts
O=work/run
mkdir -p $O
glen() { case $1 in 1) echo 32558;; 2) echo 29214;; 3) echo 31230;; esac; }   # genome lengths (composition.txt)

# 0. inputs
$MP --version > $O/miniprot_version.txt
python3 -I $S/composition.py genome_1.fna genome_2.fna genome_3.fna > $O/composition.txt
python3 -I $S/prep_refs.py reference_proteins.faa > $O/refs.faa      # refA_cob x3 -> refA_cob_v1/v2/v3

for g in 1 2 3; do
  # 1. linear genomes, default settings, yeast mito code
  $MP -T 3 -t 1 -u genome_$g.fna $O/refs.faa > $O/g${g}_linear.paf 2> $O/g${g}_linear.log
  # 2. circular handling: append first 8 kb (genes crossing the sequence end become contiguous)
  python3 -I $S/circularize.py genome_$g.fna 8000 > $O/genome_${g}_circ.fna
  $MP -T 3 -t 1 -u $O/genome_${g}_circ.fna $O/refs.faa > $O/g${g}_circ.paf 2> $O/g${g}_circ.log
  $MP -T 3 -t 1 --gff --aln $O/genome_${g}_circ.fna $O/refs.faa > $O/g${g}_circ.gff_aln 2>/dev/null
  # 3. independent check of the end-crossing genes: genome rotated to start at its midpoint
  python3 -I $S/rotate.py genome_$g.fna > $O/genome_${g}_rot.fna
  $MP -T 3 -t 1 $O/genome_${g}_rot.fna $O/refs.faa > $O/g${g}_rot.paf 2>/dev/null
  # 4. sensitivity sweep of the frameshift penalty (default -F 23)
  for F in 6 12 15 35; do
    $MP -T 3 -t 1 -u -F $F $O/genome_${g}_circ.fna $O/refs.faa > $O/g${g}_circ_F$F.paf 2>/dev/null
  done
done

for g in 1 2 3; do
  L=$(glen $g)
  echo "### genome_$g" 
  # 5. frameshift events (cg ops F/G) from the circular default run
  python3 -I $S/fs_events.py $O/genome_${g}_circ.fna $O/refs.faa $O/g${g}_circ.paf $L
  # 6. placement windows (all equally good positions for each +1 event)
  python3 -I $S/fs_window.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf
  # 7. in-frame codon usage per genome, start/stop check
  python3 -I $S/codons.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf
  python3 -I $S/termini.py $O/genome_${g}_circ.fna $O/refs.faa $O/g${g}_circ.paf | sort -u -k2,2
  # 8. where GGC and CTG fall relative to the reading frame
  for m in GGC CTG; do python3 -I $S/motif_scan.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf $m; done
  python3 -I $S/motif_null.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf AGGC,TGGC,GGC,CTGT,CTGC,CTG x
done > $O/analysis.txt

# 9. residue-level alignment windows around every frameshift
for g in 2 3; do python3 -I $S/aln_windows.py $O/g${g}_circ.gff_aln 40; done > $O/aln_windows.txt

# 10. codons of all three genomes at each site (reference-anchored)
SP="$O/g1_circ.paf:$O/genome_1_circ.fna:32558 $O/g2_circ.paf:$O/genome_2_circ.fna:29214 $O/g3_circ.paf:$O/genome_3_circ.fna:31230"
for x in "refA_cox1 70 86" "refA_cob_v1 130 150" "refA_cob_v1 206 220" "refA_cob_v1 74 88" "refA_cox2 218 232" "refA_cox3 116 130" "refA_cox3 244 258" "refB_atp6 184 198"; do
  echo "### $x"; python3 -I $S/compare_sites.py $O/refs.faa $x $SP
done > $O/compare_sites.txt

# 11. codons around each site (canonical placement), original-frame run-on to a stop, local GC
{
  echo "## genome_2"; python3 -I $S/site_codons.py genome_2.fna 24711 - 26496 - 26477 - 26263 - 26266 -
  python3 -I $S/stop_distance.py genome_2.fna 24711 - 26496 - 26477 - 26263 -
  python3 -I $S/local_context.py genome_2.fna 24711 - 26496 - 26477 - 26263 -
  echo "## genome_3"; python3 -I $S/site_codons.py genome_3.fna 121 + 27197 + 9815 + 10197 + 10200 + 19803 + 19800 +
  python3 -I $S/stop_distance.py genome_3.fna 121 + 27197 + 9815 + 10197 + 19803 +
  python3 -I $S/local_context.py genome_3.fna 121 + 27197 + 9815 + 10197 + 19803 +
} > $O/sites.txt

# 12. refA_atp6 did not map to genome_3: retry with more sensitive settings
awk '/^>refA_atp6$/{p=1;print;next} /^>/{p=0} p' $O/refs.faa > $O/refA_atp6.faa
for opt in "" "-M 0" "-k 5" "-M 0 -k 5 -n 1" "-S"; do
  echo "opts=[$opt]"; $MP -T 3 -t 1 -u $opt $O/genome_3_circ.fna $O/refA_atp6.faa 2>/dev/null | cut -f1-12
done > $O/refA_atp6_retry.txt

# 13. cob intron-ORF parts of refA_cob_v2/v3 mapped on their own (no splicing)
python3 -I - "$O/refs.faa" > $O/intron_orfs.faa <<'PY'
import sys
seqs = {}; n = None
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith('>'): n = l[1:]; seqs[n] = ''
    else: seqs[n] += l
print('>refA_cob_v2_intronORF_172-517\n' + seqs['refA_cob_v2'][171:])
print('>refA_cob_v3_intronORF_143-423\n' + seqs['refA_cob_v3'][142:])
PY
for g in 1 2 3; do $MP -T 3 -t 1 -u -S $O/genome_${g}_circ.fna $O/intron_orfs.faa 2>/dev/null | cut -f1-12,14-19; done > $O/intron_orfs.paf
echo done
