"""Cold review 2, check 8 (POST-HOC): reproducible versions of numbers that were produced by unshipped commands. Locked outputs are only read.
 b1: B1 recalculation (REPORT 4.4)   cols: shared UAA columns (REPORT 7.2.3)   stem: bird/turtle re-inclusion with and without AP002973.1   elig: the 'eligible' field of results/insdc_only_set.tsv"""
import sys,os,json,csv,collections,math
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; what=sys.argv[1]; out=open(f"results/posthoc2/repro_{what}.txt","w"); W=lambda *x:(print(*x),print(*x,file=out))
W("POST-HOC (cold review 2, check 8).")
G=list(csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")); inc=[r for r in G if r["status"]=="included"]
if what=="b1":
    import gate3_post as P
    R={}
    for l in open("results/gate3/run/out.jsonl"):
        x=json.loads(l); R[x["acc"]]=x
    C=[r for r in inc if r["role"]=="confirmatory" and r["provenance"]!="named_extra"]
    def sh(x,mn):
        v=[t for t in (R[x["acc"]].get("first") or {}).values() if t in("TAA","TAG")]; return (sum(t=="TAG" for t in v)/len(v),len(v)) if len(v)>=mn else None
    stop=[r for r in C if r["uaa"]=="Stop" and r["acc"] in R and sh(r,3)]; t=[r for r in C if r["uaa"]=="Tyr"]
    k,s=P.sisters(t,stop); q=[sh(x,3)[0] for x in s]; b=sum(q)/len(q); me=R[t[0]["acc"]]["first"]; W("Tyr genome",t[0]["acc"],"first-stop calls",{g:c for g,c in me.items()}); v=[c for c in me.values() if c in("TAA","TAG")]
    W("called genes",len(v),"TAG",v.count("TAG"),"| sister rank",k,"sisters with >=3 called genes",len(s),"distinct taxids",len({x["taxid"] for x in s}),"mean TAG share",round(b,4)); W("P(0 TAG of %d at that share)"%len(v),round((1-b)**len(v),3))
    W("as shipped (>=3 applied to the Tyr genome too): testable",sh(t[0],3) is not None)
    W("'first' takes the first in-frame TAA/TAG past the last aligned codon; it is a terminal-codon proxy. Under UAA=Tyr a TAA there need not be a stop (REPORT 7.2.5: the COX3 transcript passes its first UAA).")
    a=[r for r in C if r["uaa"]=="Stop" and r["n33"]=="resolved"]; k2,s2=P.sisters(t,a); W("A1 (as shipped): sister rank",k2,"records",len(s2),"taxids",len({x["taxid"] for x in s2}),dict(collections.Counter(x["base"] for x in s2)))
if what=="cols":
    import gate3_post as P, uaa_caller as U, numpy as np
    me=[r for r in G if r["acc"]=="CM003310.1"][0]; U3={r["acc"]:r for r in csv.DictReader(open("results/gate3_unverified/master.tsv"),delimiter="\t")}; pg=dict(U3["ON525593.1"])
    a=P.cols(me)[0]; b=P.cols(pg)[0]; H={str(h.name):h for h in U.hmms()}; Yi=U.ABC.symbols.index("Y"); pY=lambda k: float(np.array(H[k[0]].match_emissions)[k[1]][Yi])
    com=set(a)&set(b); cy={k for k in com if pY(k)>=0.5}; ta={k for k in com if a[k]=="TAA"}; tb={k for k in com if b[k]=="TAA"}; sh=ta&tb
    W("columns aligned in both",len(com),"| R. reniformis TAA",len(ta),"P. goodeyi TAA",len(tb),"shared",len(sh),"expected (independent placement)",round(len(ta)*len(tb)/len(com),2)); W("shared columns",sorted(sh),"of which conserved-Tyr",len(sh&cy))
    A=ta&cy; B=tb&cy; W("conserved-Tyr columns in both",len(cy),"| TAA",len(A),len(B),"shared",len(A&B),"expected",round(len(A)*len(B)/len(cy),2))
    hyp=lambda N,K,n,k: sum(math.comb(K,i)*math.comb(N-K,n-i) for i in range(k,min(K,n)+1))/math.comb(N,n); W("hypergeometric P(>=shared) within conserved-Tyr columns",round(hyp(len(cy),len(A),len(B),len(A&B)),4),"(descriptive; no matched evolutionary model)")
    W("non-Tyr shared columns",sorted(sh-cy),"expected outside Tyr columns",round((len(ta)-len(A))*(len(tb)-len(B))/(len(com)-len(cy)),2))
if what=="stem":
    import trny_decide as TD, gate3_posthoc_stem as S
    R={json.loads(l)["acc"]:json.loads(l) for l in open("results/gate3/run/out.jsonl")}; FE={json.loads(l)["acc"]:json.loads(l)["f"] for l in open("data/sweep/features_fixed.jsonl")}; cls=lambda r: r["class"] or r["order"]; d=collections.defaultdict(list)
    for r in inc:
        if r["n33"]=="stem_audit_review":
            t=R[r["acc"]]["trny"]; f=FE.get(r["acc"]); yy=[y for y in (f or {}).get("trna",[]) if y[4]]; s,e=sorted((t["start"],t["end"])); st="+" if t["strand"]!="-" else "-"
            d[cls(r)].append(dict(arm=S.arm("data/sweep/fa/"+str(sum(map(ord,r["acc"]))%200).zfill(3)+f"/{r['acc']}.fna",t),callers=sum(TD.agree(t,R[r["acc"]]["second"])),annotated=(any(min(e,y[1])-max(s,y[0])>0.5*min(e-s,y[1]-y[0]) and y[2]==st for y in yy) if yy else None)))
    keep={c for c,x in d.items() if S.evaluate(c,x)[0]}; W("re-included classes",sorted(keep))
    for nm,ex in (("as shipped (named extra AP002973.1 counted)",False),("named extra excluded (REPORT 7.1, GATE_3)",True)):
        st=[r for r in inc if r["role"]=="confirmatory" and r["uaa"]=="Stop" and (r["n33"]=="resolved" or (r["n33"]=="stem_audit_review" and cls(r) in keep)) and not(ex and r["provenance"]=="named_extra")]; nu=[r for r in st if r["base"]!="U"]; W(f"{nm}: non-U33 {len(nu)}/{len(st)}; resolved only {sum(r['n33']=='resolved' and r['base']!='U' for r in st)}/{sum(r['n33']=='resolved' for r in st)}")
if what=="elig":
    R=list(csv.DictReader(open("results/insdc_only_set.tsv"),delimiter="\t")); bad=[]
    for r in R:
        q=int(r["quality"]); t=r["title"].lower(); e="yes" if (q>=1 and "unverified" not in t and "partial" not in t and "nearly complete" not in t) else "no"
        if e!=r["eligible"]: bad.append(r["accession"])
    W("rows",len(R),"eligible",sum(r["eligible"]=="yes" for r in R),"| rule: quality>=1 and title without unverified / partial / nearly complete (one-off command of 2026-11-23 20:04 UTC, not in src/). Rows that disagree with the rule:",len(bad)); W("not eligible:",dict(collections.Counter("unverified" if "unverified" in r["title"].lower() else "partial" if "partial" in r["title"].lower() else "no_complete_or_circular_flag" for r in R if r["eligible"]=="no")))
    W("records",len(R),"distinct taxids",len({r["taxid"] for r in R}),"| eligible taxids",len({r["taxid"] for r in R if r["eligible"]=="yes"}))
    u=[r for r in R if "unverified" in r["title"].lower()]; W("UNVERIFIED rows",len(u),"distinct taxids",len({r["taxid"] for r in u}))
    M={r["acc"]:r for r in csv.DictReader(open("results/gate3_unverified/master.tsv"),delimiter="\t")}; ru=[r for r in M.values() if r.get("status")=="included"]; W("UNVERIFIED run rows",len(ru),"distinct taxids",len({r["taxid"] for r in ru}))
    W("note: gate1b_sweep_set.py as committed would write the table without 'eligible'; a rerun would then make gate3_build.py fail with a KeyError (it reads r['eligible']), not silently change the set.")
