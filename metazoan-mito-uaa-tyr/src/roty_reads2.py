"""R. reniformis internal TAA/TAG sites: exact 31-mer read counts (codon at positions 15-17 of the k-mer), both strands, from independent projects. Uses src/kmer_support.py.
usage: roty_reads2.py RUN"""
import sys,os,collections; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import uaa_caller as U, kmer_support as KS
run=sys.argv[1]; fa="data/sweep/fa/"+str(sum(map(ord,"CM003310.1"))%200).zfill(3)+"/CM003310.1.fna"; n,L0,kept,recs,me,idx,H=U.analyse(fa,True); s=U.load(fa)[1]; rc=U.rc
S=s+s[:30]; c=collections.Counter(S[i:i+31] for i in range(len(S)-30)); sites=[]; lab={}
for r in recs:
    if r[1] in("TAA","TAG") and not r[5]:
        a=r[6]
        if a<14 or a+17>L0: continue
        w=s[a-14:a+17]
        if c[w]!=1 or c[rc(w)]!=0: continue
        coding=r[1]; tag=f"{r[0]}:{r[2]}:{r[6]}:{coding}:{'pY%.2f'%r[3]}"
        for nm,cd in (("genome",coding),("TAT","TAT"),("TAC","TAC"),("TAA","TAA"),("TAG","TAG")):
            g=cd if r[7]=="+" else rc(cd); k=w[:14]+g+w[17:]; lab.setdefault((tag,nm),[]).append(k)
cnt,err=KS.support(run,{f"{t}|{nm}":ks for (t,nm),ks in lab.items()},maxgb=20)
if cnt is None: print("run",run,err); sys.exit()
tags=sorted({t for t,_ in lab}); print("run",run,"sites with unique k-mers",len(tags))
for t in tags: print(t,{nm:cnt[f"{t}|{nm}"] for nm in("genome","TAA","TAT","TAC","TAG")})
