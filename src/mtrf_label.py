"""Label candidates by best human paralog (blastp vs 4 human refs); output results/work/mtrf_labels.tsv and per-class FASTA."""
import subprocess,csv,collections,sys
B="env/prefix/bin/"
FA=sys.argv[1]; OUT=sys.argv[2]
subprocess.run([B+"makeblastdb","-in","data/mtrf/human_refs.faa","-dbtype","prot","-out","results/work/db_hs"],capture_output=True)
r=subprocess.run([B+"blastp","-query",FA,"-db","results/work/db_hs","-evalue","1e-3","-outfmt","6 qseqid sseqid pident length evalue bitscore qstart qend","-num_threads","4"],capture_output=True,text=True)
D=collections.defaultdict(dict)
for l in r.stdout.splitlines():
    q,s,p,ln,e,b,*_=l.split("\t")
    if s not in D[q] or float(b)>D[q][s][0]: D[q][s]=(float(b),float(e),float(p))
L={}
for q,h in D.items():
    s=sorted(h.items(),key=lambda x:-x[1][0]); L[q]=(s[0][0],s[0][1][0],s[1][0] if len(s)>1 else "",s[1][1][0] if len(s)>1 else 0)
w=csv.writer(open(OUT,"w"),delimiter="\t"); w.writerow("protein best_human best_bits second_human second_bits".split())
for q,v in L.items(): w.writerow([q,v[0],v[1],v[2],v[3]])
print(collections.Counter(v[0] for v in L.values()))
