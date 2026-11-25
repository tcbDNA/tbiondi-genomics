#!/bin/sh
# Cold review 2, check 6 (post-hoc). Competitive mapping of SRR26324223 to CM003310.1 + GCA_001026735.1 (mitochondrial contigs removed), then local assembly. Criteria: results/posthoc2/roty_criteria.md (written first).
cd "$(dirname "$0")/.."; W=results/work/r2/roty; M=env/tools/mapenv/bin; B=env/prefix/bin
MT=$W/mt.fa; cp data/sweep/fa/$(env/prefix/bin/python -c "print(str(sum(map(ord,'CM003310.1'))%200).zfill(3))")/CM003310.1.fna $MT
# nuclear contigs that are the mitochondrion: blastn >=98% identity over >=50% of the contig
$B/makeblastdb -in data/mtrf/g4/Rotylenchulus_reniformis.fna -dbtype nucl -out $W/nuc >/dev/null
$B/blastn -task megablast -query $MT -db $W/nuc -evalue 1e-20 -dust no -max_target_seqs 5000 -outfmt "6 sseqid pident length qstart qend sstart send slen" > $W/mt_vs_nuc.tsv
env/prefix/bin/python -I - <<'P'
import collections
W="results/work/r2/roty/"; c=collections.defaultdict(list); sl={}
for l in open(W+"mt_vs_nuc.tsv"):
    f=l.split("\t"); c[f[0]].append((float(f[1]),int(f[2]),int(f[3]),int(f[4]),int(f[5]),int(f[6]))); sl[f[0]]=int(f[7])
def cov(x,n):
    iv=sorted((min(a,b),max(a,b)) for _,_,_,_,a,b in x); m=e=0
    for a,b in iv:
        if b>e: m+=b-max(a,e+1)+1; e=b
    return m/n
mt=[k for k,x in c.items() if cov([y for y in x if y[0]>=98],sl[k])>=0.5]; print("contigs with mt similarity",len(c),"| treated as mitochondrial and removed:",[(k,sl[k],round(cov(c[k],sl[k]),2)) for k in mt]); print("others:",sorted(((k,sl[k],max(y[1] for y in c[k]),max(y[0] for y in c[k])) for k in c if k not in mt),key=lambda z:-z[2])[:12])
keep=True; o=open(W+"ref.fa","w"); o.write(open(W+"mt.fa").read())
for l in open("data/mtrf/g4/Rotylenchulus_reniformis.fna"):
    if l[0]==">": keep=l[1:].split()[0] not in mt
    if keep: o.write(l)
P
$B/fastq-dump --split-3 --skip-technical -O $W results/work/SRR26324223.sra >/dev/null 2>&1; ls -la $W/*.fastq
$M/minimap2 -t 6 -ax sr -N 10 --secondary=yes $W/ref.fa $W/SRR26324223_1.fastq $W/SRR26324223_2.fastq 2> $W/mm.log | $M/samtools sort -@ 3 -m 1G -o $W/all.bam - && $M/samtools index $W/all.bam && echo MAPPED
