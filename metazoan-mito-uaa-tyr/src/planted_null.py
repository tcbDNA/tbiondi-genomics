"""Planted-null sensitivity: replace TAT/TAC by TAA at conserved-Tyr columns (p(Y)>=0.5, unmasked) of standard-code controls at graded fractions."""
import sys,csv,random,json,os
sys.path.insert(0,"src"); import uaa_caller as U
FR=[0.0,0.05,0.10,0.25,0.50,1.0]; REPS=5; SEED=20261123
def plan():
    P=[r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")]
    matched=[r for r in P if r["role"]=="neg_matched"]; broad=[r for r in P if r["role"]=="neg_broad"]
    rng=random.Random(SEED); by={}
    for r in broad: by.setdefault(r["lineage"],[]).append(r)
    pick=[]
    for l in sorted(by): rng.shuffle(by[l])
    k=0
    while len(pick)<40 and k<40:
        for l in sorted(by):
            if k<len(by[l]) and len(pick)<40: pick.append(by[l][k])
        k+=1
    return [r["panel_id"] for r in matched]+[r["panel_id"] for r in pick],P
def one(a):
    pid,topo,frac,rep=a
    fa=f"data/panel/{pid}.fna"; circ=topo!="linear"
    name,L0,kept,recs,me,idx,H=U.analyse(fa,circ)
    sites=[r for r in recs if not r[5] and r[3]>=0.5 and r[1] in ("TAT","TAC")]
    rng=random.Random(f"{SEED}|{pid}|{frac}|{rep}"); n=len(sites); k=0 if frac==0 else max(1,round(frac*n)); sel=rng.sample(sites,min(k,n))
    seq=U.load(fa)[1]; s=list(seq)
    for r in sel:
        np_,st=r[6],r[7]
        if np_+3>len(s): continue
        cod="".join(s[np_:np_+3]) if st=="+" else U.rc("".join(s[np_:np_+3]))
        s[np_:np_+3]=list("TAA" if st=="+" else U.rc("TAA"))
    tf=f"results/work/planted/{pid}__{frac}__{rep}.fna"; open(tf,"w").write(f">{pid}\n{''.join(s)}\n")
    # summarise from a temp file whose id is panel-registered
    d=U.summarise(tf,circ); d.update(panel_id=f"PLANTED_{pid}",base=pid,frac=frac,rep=rep,sites=n,planted=len(sel)); return d
if __name__=="__main__":
    from concurrent.futures import ProcessPoolExecutor
    os.makedirs("results/work/planted",exist_ok=True); ids,P=plan(); topo={r["panel_id"]:r["topology"] for r in P}
    jobs=[(i,topo[i],f,r) for i in ids for f in FR for r in range(1 if f in (0.0,1.0) else REPS)]
    print("genomes",len(ids),"jobs",len(jobs),flush=True)
    with ProcessPoolExecutor(8) as ex, open("results/work/planted_panel.jsonl","w") as o:
        for d in ex.map(one,jobs,chunksize=4): o.write(json.dumps(d)+"\n")
    print("done")
