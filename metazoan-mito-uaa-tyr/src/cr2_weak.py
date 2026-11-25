"""Cold review 2, item 10 (POST-HOC): weak-tier audit under stricter bookkeeping: relatives restricted to distinct taxids (and not the query's taxid), >=8 of the 10 relatives must align at a site, >=80% of those Tyr. Locked weak_audit is imported, not changed. writes results/posthoc2/weak_audit_strict.txt"""
import sys,os,csv,json; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import gate3_post as P
G=list(csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")); inc=[r for r in G if r["status"]=="included"]; S=[r for r in inc if r["provenance"]!="named_extra" and r["uaa"]=="Stop"]; oth=[r for r in inc if r["uaa"] in("Tyr","weak")]
PA={r["acc"]:r["path"] for r in csv.DictReader(open("results/gate3_unverified/manifest.tsv"),delimiter="\t")}; U={r["acc"]:r for r in csv.DictReader(open("results/gate3_unverified/master.tsv"),delimiter="\t")}
o=open("results/posthoc2/weak_audit_strict.txt","w"); W=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True)); W("POST-HOC (cold review 2, item 10)")
def distinct(r,n=10):
    out=[];seen={r["acc"]};tx={r["taxid"]}
    for k in P.RANKS:
        if not r.get(k): continue
        for x in sorted((x for x in S if x.get(k)==r[k] and x["acc"] not in seen),key=lambda x:x["acc"]):
            if x["taxid"] in tx: continue
            out.append((k,x)); seen.add(x["acc"]); tx.add(x["taxid"])
            if len(out)==n: return out
    return out
def sites(r,rel,need):
    m,tag,fr=P.cols(r); RM=[P.cols(x)[0] for _,x in rel]; ok=[]
    for g,c in [k for k,v in m.items() if v=="TAA"]:
        v=[q[(g,c)] for q in RM if (g,c) in q]
        if len(v)<need or sum(x in("TAT","TAC") for x in v)/len(v)<0.8: continue
        if len(fr[g])>1 or any(gg==g and abs(cc-c)<=30 for gg,cc in tag): continue
        ok.append((g,c))
    return ok
L=[r for r in inc if r["uaa"]=="weak" and r["role"]=="confirmatory"]+[dict(U["ON525593.1"],path=PA["ON525593.1"],uaa="weak")]
for r in L:
    a=P.closest(r,S); b=distinct(r); s1=sites(r,a,5); s2=sites(r,b,8); s3=sites(r,[x for x in b if x[1]["role"]=="confirmatory"],8) if sum(x[1]["role"]=="confirmatory" for x in b)>=8 else None
    W(r["acc"],r["species"],"| locked relatives",len(a),"distinct taxids",len({x['taxid'] for _,x in a}),"-> supported sites",len(s1),"genes",len({g for g,_ in s1}),"| strict relatives",len(b),"-> sites",len(s2),"genes",len({g for g,_ in s2}),sorted(s2)); W("    strict relatives:",[(k,x["species"]) for k,x in b])
    W("    same-family replication (locked rule):",P.weak_audit(r,S,oth+[r])["replication"] or "none","| verdict",P.weak_audit(r,S,oth+[r])["verdict"])
