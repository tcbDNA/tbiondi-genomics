"""Candidate read checks (post-hoc, after Gate 3). For one mitogenome record: counts reads carrying the record's N33 31-mer, the same 31-mer with tRNA U33, and the 31-mers
with the other two bases, plus 3 control 31-mers from the same mitogenome (to show the run contains that mtDNA). Exact matches, both strands.
usage: gate3b_reads.py ACC RUN [RUN...]   (extra alt k-mers may be given as LABEL=KMER)"""
import sys,os,csv,json; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import kmer_support as KS
acc=sys.argv[1]; runs=[a for a in sys.argv[2:] if "=" not in a]; extra=dict(a.split("=") for a in sys.argv[2:] if "=" in a)
def load(a):
    for b in ("data/sweep/fa","data/sweep_unverified/fa","data/sweep/posthoc"):
        p=f"{b}/{str(sum(map(ord,a))%200).zfill(3)}/{a}.fna" if "posthoc" not in b else f"{b}/{a}.fna"
        if os.path.exists(p): return "".join(l.strip() for l in open(p) if l[0]!=">").upper()
s=load(acc); d=None
for fn in ("results/gate3/run/out.jsonl","results/gate3/posthoc_tetillidae/out.jsonl","results/gate3_unverified/run/out.jsonl"):
    if os.path.exists(fn):
        for l in open(fn):
            x=json.loads(l)
            if x["acc"]==acc: d=x["trny"]
g=int(d["g34_pos"]); m=d["strand"]=="-"; site=g if m else g-2; w=lambda p:s[p-15:p+16]; lab={}
for b in "ACGU":
    plus=(b.replace("U","T")).translate(str.maketrans("ACGT","TGCA")) if m else b.replace("U","T"); lab[f"N33={b}"]=[s[site-15:site]+plus+s[site+1:site+16]]
L=len(s)
for i,f in enumerate((0.15,0.45,0.75)):
    p=int(L*f); k=w(p)
    while s.count(k)!=1: p+=40; k=w(p)
    lab[f"control{i+1}@{p}"]=[k]
for k,v in extra.items(): lab[k]=[v]
print(acc,"record N33",d["N33"],"strand",d["strand"],"loop",d["loop"],"| 31-mers",{k:v[0] for k,v in lab.items() if k.startswith("N33")},flush=True)
for r in runs:
    c,e=KS.support(r,lab,maxgb=25); print("  run",r,{k:c[k] for k in lab} if c is not None else e,flush=True)
