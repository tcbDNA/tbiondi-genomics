"""Pre-registered Gate 3 follow-ups whose code was written AFTER the PREREG lock (rules locked, code not) and BEFORE any sweep result was opened:
 weak-tier audit (2.1), A1 and B1 (section 5), T6 checks (NUMT, second trnY). Origin clusters for A1/B1 = families (genus if no family) holding >=1 primary Tyr species.
usage: gate3_post.py  (reads results/mitogenome_master.tsv and results/gate3/run/out.jsonl)"""
import sys,csv,json,collections,math
sys.path.insert(0,"src"); import uaa_caller as U
RANKS=("genus","family","order","class","phylum")
def path(r): return r.get("path") or "data/sweep/fa/"+str(sum(map(ord,r["acc"]))%200).zfill(3)+f"/{r['acc']}.fna"
_c={}
def cols(r):
    """(gene,col)->(codon,frame) over unmasked aligned columns; plus TAG columns and frames per gene"""
    if r["acc"] in _c: return _c[r["acc"]]
    n,L0,kept,recs,*_=U.analyse(path(r),r["topology"]!="linear"); fr=collections.defaultdict(set)
    for k in kept: fr[k["gene"]].add((k["strand"],k["frame"]))
    m={(x[0],x[2]):x[1] for x in recs if not x[5]}; _c[r["acc"]]=(m,[(x[0],x[2]) for x in recs if x[1]=="TAG" and not x[5]],fr); return _c[r["acc"]]
def closest(r,pool,n=10):
    out=[]; seen={r["acc"]}
    for k in RANKS:
        if not r.get(k): continue
        for x in sorted((x for x in pool if x.get(k)==r[k] and x["acc"] not in seen),key=lambda x:x["acc"]):
            out.append((k,x)); seen.add(x["acc"])
            if len(out)==n: return out
    return out
def weak_audit(r,stop_pool,others):
    rel=closest(r,stop_pool); m,tag,fr=cols(r); taa=[k for k,c in m.items() if c=="TAA"]
    if len(rel)<5: return dict(acc=r["acc"],relatives=len(rel),verdict="manual_alignment_audit (fewer than 5 Stop relatives)",sites=len(taa))
    RM=[cols(x)[0] for _,x in rel]; ok=[]
    for g,c in taa:
        v=[q[(g,c)] for q in RM if (g,c) in q]
        if len(v)<5 or sum(x in("TAT","TAC") for x in v)/len(v)<0.8: continue
        if len(fr[g])>1 or any(gg==g and abs(cc-c)<=30 for gg,cc in tag): continue
        ok.append((g,c))
    i=len({g for g,_ in ok})>=2; rep=""
    if i:
        for x in others:
            if x["acc"]==r["acc"] or not r.get("family") or x.get("family")!=r["family"] or x["species"]==r["species"]: continue
            if x["uaa"]=="Tyr": rep=x["acc"]+" primary"; break
            mm=cols(x)[0]
            if x["uaa"]=="weak" and any(mm.get(k)=="TAA" for k in ok): rep=x["acc"]+" weak, shared column"; break
    return dict(acc=r["acc"],species=r["species"],relatives=len(rel),closest_rank=rel[0][0],farthest_rank=rel[-1][0],sites=len(taa),supported_sites=len(ok),supported_genes=len({g for g,_ in ok}),replication=rep,verdict="supported_weak" if i and rep else "alignment_ok_no_replication" if i else "rejected")
LEVEL="family"   # "family" = most origins, "order" = fewest; both are reported (PREREG: origins as a range)
def clusters(C):
    d=collections.defaultdict(list)
    for r in C:
        if r["uaa"]=="Tyr": d[(r["phylum"],r["class"],r["order"])+((r["family"] or r["genus"],) if LEVEL=="family" or not r["order"] else ())].append(r)
    return d
def sisters(members,pool):
    r=members[0]; mem={x["acc"] for x in members}
    for k in RANKS:
        if len({x.get(k) for x in members})>1: continue      # only ranks shared by every member
        s=[x for x in pool if r.get(k) and x.get(k)==r[k] and x["acc"] not in mem]
        if s: return k,s
    return None,[]
def binom2(k,n):
    if n==0: return None
    p=[math.comb(n,i)/2**n for i in range(n+1)]; return min(1.0,sum(x for x in p if x<=p[k]+1e-12))
def A1(C,W):
    stop=[r for r in C if r["uaa"]=="Stop" and r["n33"]=="resolved"]; u=n=0
    for c,m in clusters(C).items():
        k,s=sisters(m,stop); b=collections.Counter(x["base"] for x in s); W("  A1",c,"sister rank",k,dict(b)); 
        if s: n+=1; u+=b["U"]>sum(b.values())/2
    W(f"  A1 origins with U33-majority sisters {u}/{n}; exact binomial p={binom2(u,n)}")
def B1(C,R,W,rule="first"):
    def sh(x):
        v=[t for t in (R[x["acc"]].get(rule) or {}).values() if t in("TAA","TAG")]; return (sum(t=="TAG" for t in v)/len(v),len(v)) if len(v)>=3 else None
    stop=[r for r in C if r["uaa"]=="Stop" and r["acc"] in R and sh(r)]; up=n=0
    for c,m in clusters(C).items():
        t=[sh(x) for x in m if x["acc"] in R and sh(x)]; k,s=sisters(m,stop)
        if not t or not s: W(f"  B1[{rule}]",c,"not testable: Tyr species with >=3 called genes",len(t),"sisters",len(s)); continue
        a=sum(x[0] for x in t)/len(t); q=[sh(x)[0] for x in s]; b=sum(q)/len(q); n+=a!=b; up+=a>b; W(f"  B1[{rule}]",c,f"Tyr TAG share {a:.2f} (species {len(t)}) vs sisters ({k}, {len(q)} species) {b:.2f}")
    W(f"  B1[{rule}] origins with higher TAG share {up}/{n}; sign test p={binom2(up,n)}")
def T6(rows,R,W):
    for r in rows:
        if r["order"]=="Albuliformes" and r["status"]=="included" and r["acc"] in R:
            d=R[r["acc"]]; g={x[3] for x in d["pcg"] or []}; W("  T6",r["acc"],r["species"],"topology",r["topology"],"PCG families",len(g),"internal TAA",r["taa_int"],"internal TAG",r["tag_int"],"GUA tRNA loci ARAGORN",len(d["second"]["aragorn"]),"tRNAscan",len(d["second"]["trnascan"]),"verified GUA arms (CM)",d["trny"].get("n_ok"),"N33",r["n33"],r["base"])
if __name__=="__main__":
    rows=list(csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")); R={json.loads(l)["acc"]:json.loads(l) for l in open("results/gate3/run/out.jsonl")}
    inc=[r for r in rows if r["status"]=="included"]; C=[r for r in inc if r["role"]=="confirmatory" and r["provenance"]!="named_extra"]; stop=[r for r in inc if r["uaa"]=="Stop" and r["provenance"]!="named_extra"]
    with open("results/gate3/post.txt","w") as o:
        W=lambda *x: print(*x,file=o)
        W("== weak-tier audit"); wk=[r for r in C if r["uaa"]=="weak"]; oth=[r for r in inc if r["uaa"] in("Tyr","weak")]
        for r in wk: W("  ",weak_audit(r,stop,oth))
        for lv in("family","order"):
            LEVEL=lv
            W(f"== A1 (clusters by {lv})"); A1(C,W); W(f"== B1 (clusters by {lv})"); B1(C,R,W,"first"); B1(C,R,W,"tagfirst")
        W("== T6"); T6(inc,R,W)
    print(open("results/gate3/post.txt").read()[:4000])
