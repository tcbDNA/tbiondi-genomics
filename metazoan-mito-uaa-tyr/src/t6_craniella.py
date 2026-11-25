"""Species label of the independent Craniella RNA-seq run SRR33975184: exact species-diagnostic 31-mers (present once in one Craniella mitogenome, absent from the other mitogenomes), counts in the run."""
import os,subprocess,collections,sys
sys.path.insert(0,"src"); import kmer_support as KS
def load(p): return "".join(l.strip() for l in open(p) if l[0]!=">").upper()
def path(a):
    for p in("data/sweep/fa/"+str(sum(map(ord,a))%200).zfill(3)+f"/{a}.fna",f"data/sweep/posthoc/{a}.fna"):
        if os.path.exists(p): return p
G={"C.infrequens OZ318138":"OZ318138.1","C.zetlandica OZ407450":"OZ407450.1","C.villosa PX232378":"PX232378.1","C.arb QB209682":"QB209682.1","Ci.apion OZ465198":"OZ465198.1"}
S={k:load(path(a)) for k,a in G.items()}; K=31; rc=KS.rc
sets={k:{s[i:i+K] for i in range(len(s)-K+1)} for k,s in S.items()}
lab={}; N={}
for k in sets:
    oth=set().union(*[sets[j] for j in sets if j!=k]); oth|={rc(x) for x in oth}
    d=[x for x in sorted(sets[k]) if x not in oth and set(x)<=set("ACGT")]; step=max(1,len(d)//400); pick=d[::step][:400]; N[k]=len(pick); lab[k]=pick
c,e=KS.support("SRR33975184",lab,maxgb=25)
print("SRR33975184 (file results/work/SRR33975184.sra) diagnostic 31-mers used:",N,"| total exact counts:",{k:c[k] for k in lab} if c is not None else e)
