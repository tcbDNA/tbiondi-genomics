"""EXPLORATORY / POST-HOC classification of the UNVERIFIED-flagged species. Uses the LOCKED functions (gate3_classify.uaa_class, n33; gate3_post.weak_audit) unchanged; writes only to results/gate3_unverified/.
Never merged with Gate 3 totals. Relatives for the weak audit come from the Gate 3 master table (Stop, included). usage: gate3u_classify.py [--count-only]"""
import sys,csv,json,collections,math,os
sys.path.insert(0,"src"); import gate3_classify as K, gate3_post as P
D="results/gate3_unverified/"; FE={}
for l in open("data/sweep_unverified/features.jsonl"): x=json.loads(l); FE[x["acc"]]=x["f"]
R={}
for l in open(D+"run/out.jsonl"): x=json.loads(l); R[x["acc"]]=x
PATH={r["acc"]:r["path"] for r in csv.DictReader(open(D+"manifest.tsv"),delimiter="\t")}; rows=[]
for m in csv.DictReader(open(D+"sweep_audit.tsv"),delimiter="\t"):
    a=m["acc"]
    if m["status"]!="included": rows.append(dict(m,uaa="excluded",n33="unresolved",base="",n33_reason=m["status"])); continue
    d=R.get(a)
    if d is None: rows.append(dict(m,uaa="not_evaluated",uaa_note="missing_from_run",n33="unresolved",base="",n33_reason="missing_from_run")); continue
    c,note=K.uaa_class(d); z=K.n33(d,FE.get(a)); u=d.get("uaa") or {}
    rows.append(dict(m,path=PATH[a],uaa=c,uaa_note=note,n33=z["status"],base=z["N33"],n33_reason=z["reason"],pass0=z["pass0"],pass5=z["pass5"],pass10=z["pass10"],score=z["score"],stem=z["stem"],loop=z["loop"],n_agree=z["n_agree"],
      taa_int=u.get("n_taa_internal",""),taa_genes=u.get("n_taa_genes",""),taa_tyr=u.get("n_taa_tyrcol",""),llr=u.get("llr_Y",""),tag_int=u.get("n_tag_internal",""),aligned=u.get("aligned_codons",""),genes=u.get("n_genes_hit",""),codetta=(d.get("codetta") or {}).get("TAA",""),module_errors=";".join(sorted(d["err"])) if d.get("err") else ""))
miss=sum(r["uaa"]=="not_evaluated" and r.get("uaa_note")=="missing_from_run" for r in rows); print("rows",len(rows),"evaluated",sum(r["status"]=="included" for r in rows)-miss,"missing from run",miss)
if "--count-only" in sys.argv: sys.exit()
if miss: print("NOTE: incomplete run; results below cover evaluated genomes only")
KEYS=list(dict.fromkeys(k for r in rows for k in r)); w=csv.DictWriter(open(D+"master.tsv","w"),fieldnames=KEYS,delimiter="\t",restval=""); w.writeheader(); w.writerows(rows)
G3={r["species"] for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included"}; g3=[r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included"]
E=[r for r in rows if r["status"]=="included" and r["uaa"]!="not_evaluated"]
with open(D+"tests.txt","w") as o:
    W=lambda *x: print(*x,file=o)
    W("EXPLORATORY / POST-HOC. UNVERIFIED-flagged INSDC species, locked pipeline. Never added to Gate 3 totals.")
    W("evaluated",len(E),"| of these, species also represented in Gate 3:",sum(r["species"] in G3 for r in E))
    W("UAA classes",dict(collections.Counter(r["uaa"] for r in E))); W("N33",dict(collections.Counter(r["n33"] for r in E)))
    st=[r for r in E if r["uaa"]=="Stop" and r["n33"]=="resolved"]; nu=[r for r in st if r["base"]!="U"]; lo,hi=K.wil(len(nu),len(st)); W(f"non-U33 among resolved Stop: {len(nu)}/{len(st)} = {100*len(nu)/max(1,len(st)):.3f}% (Wilson {100*(lo or 0):.3f}-{100*(hi or 0):.3f}%)",dict(collections.Counter(r["base"] for r in nu)))
    W("phylum (n, resolved)",{p:(sum(r['phylum']==p for r in E),sum(r['phylum']==p and r['n33']=='resolved' for r in E)) for p,_ in collections.Counter(r['phylum'] for r in E).most_common(12)})
    W("\n== Tyr / weak calls"); S=[r for r in g3 if r["uaa"]=="Stop" and r["provenance"]!="named_extra"]; oth=[r for r in E if r["uaa"] in("Tyr","weak")]
    for r in oth: W("  ",r["uaa"],r["acc"],r["species"],r["phylum"],r["order"],r["family"],"TAA",r["taa_int"],"genes",r["taa_genes"],"tyr",r["taa_tyr"],"llr",r["llr"],"TAG",r["tag_int"],"codetta",r["codetta"],"N33",r["n33"],r["base"],"also in Gate3:",r["species"] in G3,"note",r["uaa_note"])
    W("\n== weak audit (relatives = Gate 3 Stop)")
    G3o=[r for r in g3 if r["uaa"] in("Tyr","weak")]
    for r in oth:
        if r["uaa"]=="weak":
            try: W("  ",P.weak_audit(r,S,G3o+oth))
            except Exception as e: W("   audit error",r["acc"],str(e)[:120])
    W("\n== non-U33 among Stop, resolved")
    for r in sorted(nu,key=lambda r:(r["phylum"],r["class"],r["order"],r["species"])): W("  ",r["acc"],r["species"],r["phylum"],r["class"],r["order"],r["family"],r["base"],r["loop"],"score",r["score"],"stem",r["stem"],"agree",r["n_agree"],r["n33_reason"],"also in Gate3:",r["species"] in G3)
    W("\n== stem_audit_review non-U33 (Stop)"); [W("  ",r["acc"],r["species"],r["phylum"],r["class"],r["order"],r["base"],r["loop"],r["score"],r["stem"]) for r in E if r["n33"]=="stem_audit_review" and r["base"]!="U"]
    W("module errors",sum(bool(r.get("module_errors")) for r in E))
print(open(D+"tests.txt").read()[:2500])
