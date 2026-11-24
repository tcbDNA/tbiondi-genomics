"""EXPLORATORY: read test of every internal TAA of an UNVERIFIED hit (unique 31-mers with the codon at 15-17, exact, both strands), record codon versus TAT/TAC/TAG at the same site. usage: gate3u_reads.py ACC RUN..."""
import sys,os,collections,csv; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import uaa_caller as U, kmer_support as KS
acc,runs=sys.argv[1],sys.argv[2:]; p=[r["path"] for r in csv.DictReader(open("results/gate3_unverified/manifest.tsv"),delimiter="\t") if r["acc"]==acc][0]
n,L0,kept,recs,me,idx,H=U.analyse(p,True); s=U.load(p)[1]; rc=U.rc; S=s+s[:30]; c=collections.Counter(S[i:i+31] for i in range(len(S)-30)); lab={}
for r in recs:
    if r[1] in("TAA","TAG") and not r[5]:
        a=r[6]
        if a<14 or a+17>L0: continue
        w=s[a-14:a+17]
        if c[w]!=1 or c[rc(w)]!=0: continue
        tag=f"{r[0]}:{r[2]}:{r[1]}:pY{r[3]:.2f}"
        for nm,cd in(("genome",r[1]),("TAT","TAT"),("TAC","TAC"),("TAA","TAA"),("TAG","TAG")): lab[(tag,nm)]=[w[:14]+(cd if r[7]=="+" else rc(cd))+w[17:]]
ctl=[]; step=max(1,len(s)//6)
for i in range(step,len(s)-40,step):
    if c[s[i:i+31]]==1: lab[("control@%d"%i,"genome")]=[s[i:i+31]]
print(acc,"internal TAA/TAG sites with unique 31-mers",len({t for t,k in lab if not t.startswith("control")}),flush=True)
for run in runs:
    cnt,err=KS.support(run,{f"{t}|{k}":v for (t,k),v in lab.items()},maxgb=20)
    if cnt is None: print("  run",run,err); continue
    print("  run",run,"controls",[cnt[f"{t}|genome"] for t,k in lab if t.startswith("control")])
    for t in sorted({t for t,k in lab if not t.startswith("control")}): print("    ",t,{k:cnt[f"{t}|{k}"] for k in("genome","TAT","TAC","TAA","TAG")},flush=True)
