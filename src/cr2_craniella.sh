#!/bin/sh
# Cold review 2, check 2 (post-hoc): map each Craniella run to its own assembly (mitogenome + nuclear, as deposited) and keep the alignments; HiFi with map-hifi, RNA-seq with sr.
cd "$(dirname "$0")/.."; W=results/work/r2/cran; mkdir -p $W; M=env/tools/mapenv/bin; B=env/prefix/bin
map(){ # LABEL RUN ASM PRESET
 [ -s $W/$1.bam.bai ] && return
 $B/fastq-dump -Z --fasta 0 --split-spot results/work/$2.sra 2>/dev/null | $M/minimap2 -t 8 -ax $4 --secondary=no results/work/numt/$3.fna - 2> $W/$1.log | $M/samtools view -b -F 4 - | $M/samtools sort -@ 3 -m 1G -o $W/$1.bam - && $M/samtools index $W/$1.bam; echo "$(date +%T) $1 done"; }
map infr_hifi ERR14835550 GCA_974939545.1 map-hifi
map zetl_hifi ERR14842678 GCA_978397455.1 map-hifi
map infr_rna SRR33975184 GCA_974939545.1 sr
map zetl_rna ERR14857142 GCA_978397455.1 sr
echo CRAN_DONE
