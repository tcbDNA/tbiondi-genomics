"""Short reads: extract reads with any 21-mer near the site (not overlapping N33), align to a 400-nt window with blastn-short rules, tabulate the base at N33.
usage: t6_short.py RUN.sra ACC G34 STRAND LABEL"""
import sys,os,subprocess,collections
run,acc,g34,strand,label=sys.argv[1:6]; g34=int(g34); B="env/prefix/bin/"; W="results/work/t6/"; os.makedirs(W,exist_ok=True)
RC=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(RC)[::-1]
def load(a):
    p="data/sweep/fa/"+str(sum(map(ord,a))%200).zfill(3)+f"/{a}.fna"; return "".join(l.strip() for l in open(p) if l[0]!=">").upper()
mt=load(acc); pos=g34+1 if strand=="-" else g34-1; K=21
ks=set()
for i in list(range(pos-1-75,pos-1-K-2))+list(range(pos+2,pos+75)):
    if 0<=i and i+K<=len(mt) and not(i<=pos-1<i+K): ks.add(mt[i:i+K])
open(W+label+".pat","w").write("\n".join(sorted(ks|{rc(k) for k in ks}))+"\n"); out=W+label+".reads.fa"
if not os.path.exists(out) or os.path.getsize(out)==0: subprocess.run(["bash","-c",f"{B}fastq-dump -Z --fasta 0 --split-spot {run} 2>/dev/null | env/tools/kreads/kreads {W}{label}.pat > {out}"])
n=sum(1 for l in open(out) if l[0]==">"); print(label,"reads with a near-site 21-mer:",n)
w0=max(1,pos-200); open(W+label+".win.fa","w").write(f">w\n{mt[w0-1:pos+200]}\n"); site=pos-w0+1
subprocess.run([B+"makeblastdb","-in",out,"-dbtype","nucl","-out",W+label],capture_output=True)
r=subprocess.run([B+"blastn","-task","blastn","-query",W+label+".win.fa","-db",W+label,"-evalue","1e-3","-word_size","11","-max_target_seqs","500000","-num_threads","6","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout
for mid,ml in((98,60),(95,50),(90,45)):
    c=collections.Counter(); seen=set()
    for l in sorted((x.split("\t") for x in r.splitlines()),key=lambda x:-int(x[2])):
        s,pi,ln,qs,qe=l[0],float(l[1]),int(l[2]),int(l[3]),int(l[4])
        if ln<ml or pi<mid or s in seen or not(qs<=site-8 and qe>=site+8): continue
        seen.add(s); p=qs; b="-"
        for x,y in zip(l[7],l[8]):
            if x!="-":
                if p==site: b=y; break
                p+=1
        c[{"A":"U","C":"G","G":"C","T":"A"}.get(b,b) if strand=="-" else b.replace("T","U")]+=1
    print(f"  identity>={mid}% length>={ml}: reads {sum(c.values())}, tRNA base {dict(c)}")
