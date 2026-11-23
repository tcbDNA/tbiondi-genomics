"""Auxiliary mitochondrial checks (panel-guarded): terminal stops per PCG and codon O/E.
terminal(): walk in frame from the last profile-aligned codon of each gene (3'-most domain) to the next downstream gene start
 (same-strand PCG domain or tRNA; ARAGORN and tRNAscan-SE). Allowed overlap OV nt. Result per gene: TAA, TAG, T--(incomplete), none, unknown(no boundary).
 CODE-AGNOSTIC modes, applied identically to every genome (the only ones used for tests):
  'first' = first in-frame TAA/TAG after the last aligned codon (sense TAA in an unaligned tail counts as TAA: biased AGAINST a TAG shift in Tyr genomes);
  'tagfirst' = first in-frame TAG before boundary+OV, else a TAA abutting the boundary (biased TOWARD TAG).
 'Stop'/'Tyr' modes are code-aware and circular for the TAG-shift test; kept only for comparison with annotation."""
import sys,os,json,collections,math
sys.path.insert(0,"src"); import panel_guard, uaa_caller as U
OV=12; CAP=80
def feats(kept,fr,tr,L0,strand,H):
    st=[]
    for d in kept:
        s=d["ns"]-3*(d["hf"]-1) if d["strand"]=="+" else (d["ne"]+3*(d["hf"]-1))
        st.append((d["strand"],s if d["strand"]=="+" else L0-s,d["gene"]))
    for a,b,s,*_ in tr: st.append((s,a-1 if s=="+" else L0-b,"tRNA"))
    return [(p,g) for s,p,g in st if s==strand]
def terminal(fa,circular,mode,trna):
    name,L0,kept,recs,me,idx,H=U.analyse(fa,circular); seq=U.load(fa)[1]; fr=U.frames(seq,2000 if circular else 0); M={str(h.name):h.M for h in H}
    tr=[(a,b,s) for a,b,s,*_ in trna.get("aragorn",[]) if a>0]+[(a,b,s) for a,b,s,*_ in trna.get("trnascan",[])]
    out={}; genes=collections.defaultdict(list)
    for d in kept: genes[d["gene"]].append(d)
    for g,ds in genes.items():
        d=max(ds,key=lambda d:d["ht"]); st,f,cod,x,aa=fr[d["frame"]]; i=d["tt"]; rem=M[g]-d["ht"]
        if rem>30: out[g]="unaligned_3prime"; continue
        p0=f+3*i                                        # nt index in strand sequence x right after the last aligned codon
        off=0 if st=="+" else len(x)-L0                 # x for '-' = rc(seq+wrap); convert feature coords
        F=sorted(set(p+off if st=="-" else p for p,gg in feats(kept,fr,tr,L0,st,H) if gg!=g))
        F+= [p+L0 for p in F]; nxt=[p for p in F if p>=p0-15]; bound=nxt[0] if nxt else None
        if bound is None or bound-p0>3*(CAP+rem): b=None
        else: b=bound
        lim=(b+OV) if b is not None else p0+3*(rem+40); t=None; tag=None; taa_ab=None; k=i
        while k<len(cod) and f+3*k<lim:
            c=cod[k]; e=f+3*k+3
            if c in("TAA","TAG") and t is None: t=c
            if c=="TAG" and tag is None: tag=c
            if c=="TAA" and b is not None and b-4<=e<=b+OV: taa_ab=c
            k+=1
        if mode in("Stop","first"): r=t
        else: r=tag or taa_ab
        if r is None:
            if b is None: r="unknown"
            else:
                n=(b-p0)%3; tail=x[b-n:b] if n else ""; r={"T":"T--","TA":"TA-"}.get(tail,"none")
        out[g]=r
    return out
FOUR=[a+b for a,b in ("GC","GG","CC","CG","AC","GT","CT","TC")]
def codon_oe(fa,circular):
    name,L0,kept,recs,me,idx,H=U.analyse(fa,circular); c=collections.Counter(r[1] for r in recs if not r[5] and set(r[1])<=set("ACGT")); n=sum(c.values())
    F=[collections.Counter() for _ in range(3)]
    for k,v in c.items():
        for i in range(3): F[i][k[i]]+=v
    E=lambda k: n*math.prod(F[i][k[i]]/n for i in range(3)); oe={k:(round(c[k]/E(k),3) if E(k)>0 else None) for k in("TAA","TAT","TAC","TAG")}
    f4=collections.Counter()
    for k,v in c.items():
        if k[:2] in FOUR: f4[k[2]]+=v
    s=sum(f4[b] for b in "ATC"); exp=f4["A"]/s if s else None; y=c["TAA"]+c["TAT"]+c["TAC"]
    return dict(n_codons=n,TAA=c["TAA"],TAT=c["TAT"],TAC=c["TAC"],oe=oe,gc=None,a3_4fold=round(f4["A"]/max(1,sum(f4.values())),3),exp_taa_frac=round(exp,3) if exp is not None else None,obs_taa_frac=round(c["TAA"]/y,3) if y else None)
def one(a):
    f,circ,mode,trna=a; pid=os.path.basename(f).replace(".fna",""); d=dict(panel_id=pid)
    try: d["terminal"]=terminal(f,circ,mode,trna); d["first"]=terminal(f,circ,"first",trna); d["tagfirst"]=terminal(f,circ,"tagfirst",trna); d["oe"]=codon_oe(f,circ)
    except Exception as e: d["err"]=f"{type(e).__name__}: {e}"[:200]
    return d
if __name__=="__main__":
    import csv; from concurrent.futures import ProcessPoolExecutor
    P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; panel_guard.check(P)
    TR={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/trna_all.jsonl")}
    jobs=[(f"data/panel/{k}.fna",r["topology"]!="linear","Tyr" if r["role"]=="pos" else "Stop",TR.get(k,{})) for k,r in P.items()]
    with ProcessPoolExecutor(int(sys.argv[2])) as ex,open(sys.argv[1],"w") as o:
        for d in ex.map(one,jobs,chunksize=4): o.write(json.dumps(d)+"\n")
    print("finished",len(jobs))
