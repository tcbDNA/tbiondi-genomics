"""Independent-origin counter on the NCBI taxonomy tree (all ranks from nodes.dmp).
 changes = minimum state changes (Hartigan parsimony, polytomies allowed, gains and losses equal).
 gains_min / gains_max = origins under no-reversal: gains_min treats every polytomy as resolvable in favour of one origin, gains_max treats it as hard.
Input: dict taxid -> 0/1 ('Stop'/'Tyr'). Unresolved genomes must not be passed in."""
import os,subprocess,collections
def load_nodes(p="data/ncbi/taxonomy/nodes.dmp"):
    if not os.path.exists(p): subprocess.run(["tar","xzf","new_taxdump.tar.gz","nodes.dmp","names.dmp","merged.dmp"],cwd="data/ncbi/taxonomy",check=True)
    return {a:b for a,b in ((f[0],f[1]) for f in (l.split("\t|\t",2) for l in open(p)))}
def build(tips,parent):
    ch=collections.defaultdict(set); 
    for t in tips:
        x=str(t); seen=0
        while x in parent and parent[x]!=x and seen<200: ch[parent[x]].add(x); x=parent[x]; seen+=1
    root="1"; return ch,root
def count(states,parent):
    S={str(k):v for k,v in states.items()}; ch,root=build(S,parent)
    import sys; sys.setrecursionlimit(10000)
    info={}
    def rec(n):
        # returns (Hartigan set, changes, gmin, gmax, all_tyr)
        kids=sorted(ch.get(n,())); own=S.get(n)
        items=[rec(k) for k in kids]+([({own},0,0,0,own==1)] if own is not None else [])   # a species node that also has sampled subtaxa
        if not items: return ({own},0,0,0,own==1)
        if len(items)==1 and not kids: return items[0]
        c=collections.Counter(s for it in items for s in it[0]); m=max(c.values()); A={s for s,v in c.items() if v==m}
        chg=sum(it[1] for it in items)+len(items)-m
        allt=all(it[4] for it in items)
        if allt: return (A,chg,0,0,True)                              # whole subtree Tyr: origin counted above
        full=sum(it[4] for it in items); gmin=sum(it[2] for it in items)+(1 if full else 0); gmax=sum(it[3] for it in items)+full
        return (A,chg,gmin,gmax,False)
    r=rec(root); 
    if r[4]: return dict(changes=r[1],gains_min=1,gains_max=1)
    return dict(changes=r[1],gains_min=r[2],gains_max=r[3])
def selftest():
    P={"1":"1","A":"1","B":"1","a1":"A","a2":"A","a3":"A","b1":"B","b2":"B","c":"1"}
    T=[({"a1":1,"a2":0,"b1":0,"b2":0},dict(changes=1,gains_min=1,gains_max=1)),
       ({"a1":0,"a2":0,"b1":0},dict(changes=0,gains_min=0,gains_max=0)),
       ({"a1":1,"a2":1,"b1":0,"b2":0,"c":0},dict(changes=1,gains_min=1,gains_max=1)),
       ({"a1":1,"a2":0,"b1":1,"b2":0},dict(changes=2,gains_min=2,gains_max=2)),
       ({"a1":1,"a2":1,"a3":0,"b1":0,"b2":0},dict(changes=2,gains_min=1,gains_max=2)),
       ({"a1":1,"b1":1},dict(changes=0,gains_min=1,gains_max=1))]
    bad=0
    for s,e in T:
        g=count(s,P); ok=all(g[k]==e[k] for k in e); bad+=not ok; print("ok" if ok else "MISMATCH",s,g,"expected",e)
    return bad
def brute(n=400,seed=1):
    """Exhaustive check on random trees (hard polytomies): min changes, and min 0->1 edges with no reversal (=gains_max)."""
    import random,itertools
    rng=random.Random(seed); bad=0
    for _ in range(n):
        P={"1":"1"}; inner=["1"]; tips=[]
        for i in range(rng.randint(1,4)): x=f"i{i}"; P[x]=rng.choice(inner); inner.append(x)
        for i in range(rng.randint(3,8)): x=f"t{i}"; P[x]=rng.choice(inner); tips.append(x)
        st={t:rng.randint(0,1) for t in tips}
        used=set()
        for t in tips:
            x=t
            while x!="1": used.add(x); x=P[x]
        nodes=[x for x in inner if x in used or x=="1"]; edges=[(P[x],x) for x in used]
        best=99; g=99
        for lab in itertools.product((0,1),repeat=len(nodes)):
            L=dict(zip(nodes,lab)); L.update(st); best=min(best,sum(L[a]!=L[b] for a,b in edges))
            if all(not(L[a]==1 and L[b]==0) for a,b in edges): g=min(g,sum(L[a]==0 and L[b]==1 for a,b in edges)+(1 if L["1"]==1 else 0))
        r=count(st,P)
        if r["changes"]!=best or r["gains_max"]!=g: bad+=1; print("BRUTE MISMATCH",st,P,r,best,g) if bad<4 else None
    print("brute force trees",n,"mismatches",bad); return bad
if __name__=="__main__": raise SystemExit(selftest()+brute())
