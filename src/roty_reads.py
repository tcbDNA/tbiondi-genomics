"""R. reniformis: exact 31-mer read counts at the internal TAA sites (and internal TAG sites) of CM003310, from runs in separate BioProjects.
For every site (codon start npos, strand) the 31-mer (codon at positions 15-17, forward-strand orientation) must be unique in the genome (both strands). Counts: genome allele, TAT, TAC, TAG/TAA alternatives.
usage: roty_reads.py RUN LIMIT_SPOTS"""
import sys,os,subprocess,collections,json; sys.path.insert(0,"src"); import uaa_caller as U
os.environ["GATE3_UNLOCK"]="1"; run=sys.argv[1]; cap=int(sys.argv[2])
fa="data/sweep/fa/"+str(sum(map(ord,"CM003310.1"))%200).zfill(3)+"/CM003310.1.fna"; n,L0,kept,recs,me,idx,H=U.analyse(fa,True); s=U.load(fa)[1]; rc=U.rc
S=s+s[:40]; sites=[]
for r in recs:
    if r[1] in("TAA","TAG") and not r[5]:
        a=r[6]; w=s[(a-14)%L0:][:0] if a<14 else s[a-14:a+17]
        if len(w)==31: sites.append((r[0],r[2],r[1],a,r[7],w,r[3]))
cnt31=collections.Counter(); import itertools
for i in range(len(S)-30): k=S[i:i+31]; cnt31[k]+=1
uniq=[x for x in sites if cnt31[x[5]]==1 and cnt31[rc(x[5])]==0]; print("sites",len(sites),"unique 31-mers",len(uniq))
pat={}
for j,x in enumerate(uniq):
    w=x[5]; codon=w[14:17]
    # alternatives: third-position / TAA<->TAG on the coding strand; for minus-strand genes the genomic codon is rc
    coding=codon if x[4]=="+" else rc(codon)
    for lab,c in (("genome",coding),("TAT","TAT"),("TAC","TAC"),("TAA","TAA"),("TAG","TAG")):
        g=c if x[4]=="+" else rc(c)
        k=w[:14]+g+w[17:]
        for kk in {k,rc(k)}: pat.setdefault(kk,[]).append((j,lab))
tmp=f"results/work/{run}.sra"; subprocess.run(["sh","src/dl.sh",f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}",tmp,"3"])
C=collections.Counter(); reads=0
p=subprocess.Popen(["env/prefix/bin/fastq-dump","-Z","--fasta","0","-X",str(cap),tmp],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
K=31; keys=set(pat); import re
for l in p.stdout:
    if l[0]==">": continue
    q=l.strip().upper(); reads+=1
    for i in range(len(q)-30):
        kk=q[i:i+31]
        if kk in keys:
            for (j,lab) in pat[kk]: C[(j,lab)]+=1
p.wait(); os.remove(tmp)
print("run",run,"reads scanned",reads)
for j,x in enumerate(uniq):
    print(x[0],"col",x[1],"genome codon",x[2],"pY",round(x[6],2),{lab:C[(j,lab)] for lab in("genome","TAT","TAC","TAA","TAG")})
