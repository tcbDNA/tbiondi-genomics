#!/bin/sh
# Cold review 2, check 2 part 3 (post-hoc): the mitogenome record is linear in the assembly, so reads over its origin go to nuclear copies that span it. Remap every read that touched the mitogenome or a nuclear copy to (nuclear contigs + the mitogenome written twice).
cd "$(dirname "$0")/.."; W=results/work/r2/cran; M=env/tools/mapenv/bin; B=env/prefix/bin
one(){ # LAB ACC ASM
 $B/blastn -task blastn -query $W/$2.q -db results/work/numt/$3 -evalue 1e-10 -dust no -max_target_seqs 2000 -num_threads 6 -outfmt "6 sseqid sstart send length" | awk '$4>=200{a=$2<$3?$2:$3;b=$2<$3?$3:$2;print $1"\t"(a-1)"\t"b}' | sort -k1,1 -k2,2n > $W/$1.bed
 $M/samtools view -b -F 2304 -L $W/$1.bed $W/$1.bam | $M/samtools fasta - 2>/dev/null > $W/$1.sub.fa; echo "$1 reads on the mitogenome or a nuclear copy: $(grep -c '>' $W/$1.sub.fa)"
 env/prefix/bin/python -I - $2 results/work/numt/$3.fna $W/$1.ref.fa <<'P'
import sys
acc,fa,out=sys.argv[1:4]; keep=True; o=open(out,"w"); s=[]
for l in open(fa):
    if l[0]==">": mt=l[1:].split()[0]==acc; 
    if mt: s+=[] if l[0]==">" else [l.strip()]
    else: o.write(l)
s="".join(s); o.write(f">{acc}_x2\n{s}{s}\n")
P
 $M/minimap2 -t 8 -ax map-hifi --secondary=no $W/$1.ref.fa $W/$1.sub.fa 2>/dev/null | $M/samtools view -b -F 4 - | $M/samtools sort -o $W/$1.x2.bam - && $M/samtools index $W/$1.x2.bam; rm -f $W/$1.ref.fa; }
one infr_hifi OZ318138.1 GCA_974939545.1
one zetl_hifi OZ407450.1 GCA_978397455.1
env/prefix/bin/python - <<'P'
import subprocess,re,collections,statistics as st,json,sys
sys.path.insert(0,"src"); M="env/tools/mapenv/bin/samtools"; W="results/work/r2/cran/"; o=open("results/posthoc2/craniella_depth3.txt","w"); Pr=lambda *x:(print(*x),print(*x,file=o)); run=lambda *a: subprocess.run([M]+list(a),capture_output=True,text=True).stdout
def base_at(f,pos):
    p=int(f[3]); i=0
    for n,k in re.findall(r"(\d+)([MIDNSHP=X])",f[5]):
        n=int(n)
        if k in"M=X":
            if p<=pos<p+n: return f[9][i+pos-p]
            p+=n;i+=n
        elif k in"IS": i+=n
        elif k in"DN":
            if p<=pos<p+n: return "-"
            p+=n
Pr("POST-HOC check 2, part 3: mitogenome written twice, competing with all nuclear contigs")
for lab,acc,L,n33 in (("infr_hifi","OZ318138.1",18285,14849),("zetl_hifi","OZ407450.1",18452,11015)):
    bam=W+lab+".x2.bam"; c=acc+"_x2"; d=[0]*(2*L+1); d20=[0]*(2*L+1)
    for q,x in (("0",d),("20",d20)):
        for l in run("depth","-aa","-G","2304","-Q",q,"-r",c,bam).splitlines(): f=l.split("\t"); x[int(f[1])]=int(f[2])
    D=[d[i]+d[i+L] for i in range(1,L+1)]; D20=[d20[i]+d20[i+L] for i in range(1,L+1)]
    Pr(f"\n{lab} {acc}: mitochondrial depth (both copies added) median",st.median(D),"10th-90th percentile",sorted(D)[L//10],sorted(D)[9*L//10],"| at N33",D[n33-1],"MAPQ>=20*",D20[n33-1])
    T=[json.loads(l) for l in open("results/gate3/run/out.jsonl") if f'"acc": "{acc}"' in l[:60]][0]; a,b=sorted((T["trny"]["start"],T["trny"]["end"])); G=sorted([(x[0],x[1],x[3]) for x in T["pcg"]],key=lambda x:min(abs(x[0]-b),abs(a-x[1])))[:4]
    Pr("  trnY",st.median(D[a-1:b]),"| neighbouring genes",[(x[2],st.median(D[x[0]-1:x[1]])) for x in sorted(G)])
    R={}
    for p in (n33,n33+L):
        for l in run("view","-F","2304",bam,f"{c}:{p}-{p}").splitlines():
            f=l.split("\t"); R[f[0]]=(base_at(f,p),int(f[4]),[int(x[5:]) for x in f[11:] if x.startswith("NM:i:")][0],sum(int(n) for n,k in re.findall(r"(\d+)([M=XD])",f[5])),sum(int(n) for n,k in re.findall(r"(\d+)([SH])",f[5])))
    Pr("  reads over N33 placed on the mitogenome:",len(R),"bases",dict(collections.Counter(v[0].replace('T','U') for v in R.values())),"| end to end (<=50 nt clipped)",sum(v[4]<=50 for v in R.values()),"| divergence from the record: median",round(st.median(v[2]/v[3] for v in R.values()),5),"reads <=0.2%:",sum(v[2]/v[3]<=.002 for v in R.values()))
    cov=[l.split("\t") for l in run("view","-F","2304","-c","-L",W+lab+".bed",bam).splitlines()]; nuc=run("view","-F","2304","-L",W+lab+".bed",bam).splitlines(); Pr("  reads still placed on nuclear copies (any part):",sum(l.split('\t')[2]!=c for l in nuc))
    # how many of those cover N33 on a nuclear copy
    import os
    q=subprocess.run(["env/prefix/bin/blastn","-task","blastn","-query",W+acc+".q","-db","results/work/numt/"+("GCA_974939545.1" if lab[0]=="i" else "GCA_978397455.1"),"-evalue","1e-10","-dust","no","-max_target_seqs","2000","-num_threads","6","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout; S={}
    for l in q.splitlines():
        f=l.split("\t")
        if f[0]==acc or int(f[2])<200: continue
        ps,k=int(f[5]),int(f[3]); stp=1 if int(f[6])>ps else -1
        for x,y in zip(f[7],f[8]):
            if k==n33 and x!="-" and y!="-": S[(f[0],ps)]=stp
            if x!="-": k+=1
            if y!="-": ps+=stp
    n=collections.Counter(); used=0
    for (ct,ps),stp in S.items():
        for l in run("view","-F","2304",bam,f"{ct}:{ps}-{ps}").splitlines():
            f=l.split("\t"); z=base_at(f,ps)
            if z and z!="-": z=z if stp>0 else {"A":"T","C":"G","G":"C","T":"A"}[z]; n[z.replace("T","U")]+=1
    Pr("  reads over N33 placed on nuclear copies:",sum(n.values()),dict(n),"across",len(S),"copies"); Pr(f"  => mitochondrial share of N33 reads {len(R)}/{len(R)+sum(n.values())} = {len(R)/(len(R)+sum(n.values())):.2f}")
P
