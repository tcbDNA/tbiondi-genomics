"""Cold review item 5 / PREREG T6: base at tRNA-Tyr N33 in long or short reads, by alignment. Reads containing any of 30 spread 25-mers of the mitogenome (either strand) are extracted, aligned with blastn
to a window around trnY (record orientation), and the base aligned to N33 is tabulated. Requires an alignment >=MINLEN bp at >=MINID% that covers N33 with 10 nt each side.
usage: t6_align.py RUN.sra ACC G34(1-based + strand) STRAND MINLEN MINID LABEL"""
import sys,os,subprocess,re,collections
run,acc,g34,strand,minlen,minid,label=sys.argv[1:8]; g34=int(g34); minlen=int(minlen); minid=float(minid); B="env/prefix/bin/"; W="results/work/t6/"; os.makedirs(W,exist_ok=True)
RC=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(RC)[::-1]
mt="".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna") if l[0]!=">").upper()
K=25; S=mt+mt[:K]; step=len(mt)//30; ks=[]
for i in range(0,len(mt)-K,step):
    for j in range(i,i+step-K):
        k=mt[j:j+K]
        if len(set(k))==4 and (mt+mt[:K]).count(k)==1: ks.append(k); break
pat=W+label+".pat"; open(pat,"w").write("\n".join(ks+[rc(k) for k in ks])+"\n")
out=W+label+".reads.fa"
if not os.path.exists(out) or os.path.getsize(out)==0:
    cmd=f"{B}fastq-dump -Z --fasta 0 {run} 2>/dev/null | env/tools/kreads/kreads {pat} > {out}"
    subprocess.run(["bash","-c",cmd])
n=sum(1 for l in open(out) if l[0]==">") if os.path.exists(out) else 0
if n==0: print(label,"NO READS MATCHED (check the filter before reading this as a result)"); print(label,"reads containing a seed k-mer:",n,"size MB",round(os.path.getsize(out)/1e6,1),flush=True)
# N33 in + strand: minus-strand gene -> complement base at g34+1; plus strand -> base at g34-1 (1-based)
pos=g34+1 if strand=="-" else g34-1
win0=max(1,pos-2500); win1=min(len(mt),pos+2500); open(W+label+".win.fa","w").write(f">w\n{mt[win0-1:win1]}\n"); site=pos-win0+1
subprocess.run([B+"makeblastdb","-in",out,"-dbtype","nucl","-out",W+label],capture_output=True)
r=subprocess.run([B+"blastn","-task","blastn","-query",W+label+".win.fa","-db",W+label,"-evalue","1e-20","-max_target_seqs","500000","-num_threads","6","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout
cnt=collections.Counter(); ids=collections.defaultdict(list); seen=set()
for l in r.splitlines():
    s,pi,ln,qs,qe,ss,se,q,h=l.split("\t"); qs,qe,ln,pi=int(qs),int(qe),int(ln),float(pi)
    if ln<minlen or pi<minid or not(qs<=site-10 and qe>=site+10) or s in seen: continue
    seen.add(s); p=qs; b="-"
    for x,y in zip(q,h):
        if x!="-":
            if p==site: b=y; break
            p+=1
    if strand=="-": b={"A":"U","C":"G","G":"C","T":"A"}.get(b,b)     # tRNA base
    else: b={"T":"U"}.get(b,b)
    cnt[b]+=1
print(label,"window site",site,"record tRNA N33:",{"A":"U","C":"G","G":"C","T":"A"}[mt[pos-1]] if strand=="-" else mt[pos-1].replace("T","U"),"| reads aligned across N33:",sum(cnt.values()),"| tRNA base counts:",dict(cnt))
