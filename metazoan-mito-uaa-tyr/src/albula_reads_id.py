"""Cold review follow-up 1: species of the PRJNA668352 fish (isolate HI-2016) from mitochondrial reads of SRR13036352. Reads with any 21-mer of CYTB or COX1 (both A. glossodonta and A. virgata seeds, so neither is favoured) are aligned to
CM165726 CYTB and COX1; a per-site majority consensus is built and compared with every Albula record."""
import sys,os,subprocess,collections,json
B="env/prefix/bin/"; W="results/work/albula/"; rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
run=sys.argv[1]; R={r[0]:r for r in json.load(open("results/work/albula_mt.json"))}
def rd(p): return "".join(l.strip() for l in open(p) if l[0]!=">").upper()
gl={g:rd(W+f"CM165726_{g}.fa") for g in("CYTB","COX1")}; vg=R["KJ910045.1"][5]
K=21; ks=set()
for s in (*gl.values(),vg): ks|={s[i:i+K] for i in range(len(s)-K+1) if set(s[i:i+K])<=set("ACGT")}
open(W+"seeds.txt","w").write("\n".join(sorted(ks|{rc(k) for k in ks}))+"\n")
if not os.path.exists(W+"hi2016.reads.fa"): subprocess.run(["bash","-c",f"{B}fastq-dump -Z --fasta 0 --split-spot {run} 2>/dev/null | env/tools/kreads/kreads {W}seeds.txt > {W}hi2016.reads.fa"])
n=sum(l[0]==">" for l in open(W+"hi2016.reads.fa")); print("reads recruited",n)
subprocess.run([B+"makeblastdb","-in",W+"hi2016.reads.fa","-dbtype","nucl","-out",W+"hi2016"],capture_output=True)
for g,s in gl.items():
    open(W+"q.fa","w").write(f">{g}\n{s}\n"); r=subprocess.run([B+"blastn","-task","blastn","-query",W+"q.fa","-db",W+"hi2016","-evalue","1e-5","-word_size","11","-max_target_seqs","5000000","-num_threads","6","-outfmt","6 sseqid pident length qstart qend qseq sseq"],capture_output=True,text=True).stdout
    C=[collections.Counter() for _ in s]; used=set()
    for l in r.splitlines():
        a,pi,ln,qs,qe,q,h=l.split("\t")
        if int(ln)<60 or float(pi)<80: continue
        if a in used: continue
        used.add(a); p=int(qs)-1
        for x,y in zip(q,h):
            if x!="-":
                if y in "ACGT": C[p][y]+=1
                p+=1
    cov=[sum(c.values()) for c in C]; con="".join(c.most_common(1)[0][0] if sum(c.values())>=10 else "N" for c in C); print(g,"aligned reads",len(used),"median depth",sorted(cov)[len(cov)//2],"sites with depth>=10:",sum(x>=10 for x in cov),"of",len(s)); d=[(i,s[i],con[i]) for i in range(len(s)) if con[i]!="N" and con[i]!=s[i]]; print("   consensus vs CM165726:",sum(con[i]!="N" for i in range(len(s))),"sites compared,",len(d),"differences",d[:6]); open(W+f"hi2016_{g}.fa","w").write(f">hi2016_{g}\n{con.replace('N','')}\n"); open(W+f"hi2016_{g}_masked.fa","w").write(f">hi2016_{g}\n{con}\n")
