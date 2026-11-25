"""[P] 25-mer patterns (both strands) around the 3' ends of the 13 protein genes of A. glossodonta CM165726.1 and L. verrillii CM163399.1, for the RNA read filter (src/kreads.c).
Gene 3' ends: profile alignment of the sweep (uaa_caller.analyse); window 260 nt upstream of the last aligned codon to 160 nt downstream. Writes results/work/polyA/<acc>.pats"""
import sys,os
sys.path.insert(0,"src"); import uaa_caller as U
for acc in ("CM165726.1","CM163399.1"):
    fa=f"data/sweep/fa/{str(sum(map(ord,acc))%200).zfill(3)}/{acc}.fna"; n,L0,kept,*_=U.analyse(fa,True); s=U.load(fa)[1]; s2=s+s[:1000]; S=set()
    for d in kept:
        a,b=(d["ne"]-260,d["ne"]+160) if d["strand"]=="+" else (d["ns"]-160,d["ns"]+260)
        w=s2[max(0,a):b]
        for x in (w,U.rc(w)):
            for i in range(len(x)-24):
                k=x[i:i+25]
                if set(k)<=set("ACGT") and max(k.count(c) for c in "ACGT")<18: S.add(k)
    open(f"results/work/polyA/{acc}.pats","w").write("\n".join(sorted(S))+"\n"); print(acc,len(kept),"gene domains",len(S),"patterns")
