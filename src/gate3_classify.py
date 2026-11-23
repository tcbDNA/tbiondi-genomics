"""Gate 3 classification, exactly as PREREG: UAA class, weak tier, trnY N33 status with locus audit, named tests. Reads out.jsonl + sweep_audit.tsv + features.jsonl.
usage: gate3_classify.py outdir [--panel]   (--panel: test mode on the panel run)"""
import sys,json,csv,math,collections,os
sys.path.insert(0,"src"); import trny_decide as TD, origins
def wil(k,n):
    if n==0: return (None,None)
    z=1.96; p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return (max(0.0,c-h),min(1.0,c+h))
def uaa_class(d):
    u,c=d.get("uaa"),d.get("codetta")
    if not u: return "not_evaluated","uaa_module_error"
    if not c or c.get("status")!="ok": cod=None
    else: cod=c.get("TAA")
    if u["aligned_codons"]<800 or u["n_genes_hit"]<5 or u["tyrcols"]<15: return "insufficient",""
    if cod is None: return "not_evaluated","codetta_error"
    T=u["n_taa_internal"]>=5 and u["n_taa_genes"]>=3 and u["n_taa_tyrcol"]>=3 and u["llr_Y"]>=3
    if T and cod=="Y": return "Tyr",""
    if u["n_taa_internal"]<=1 and u["n_taa_tyrcol"]==0 and cod!="Y": return "Stop",""
    if (u["n_taa_internal"]>=2 and u["n_taa_tyrcol"]>=1) or cod=="Y": return "weak",("counter_Tyr_codetta_not_Y" if T else "")
    return "ambiguous",""
ov=lambda a,b: min(a[1],b[1])-max(a[0],b[0])+1
def locus(d,feats):
    t=d["trny"]; a,b=sorted((int(t["start"]),int(t["end"]))); st="+" if t["strand"]!="-" else "-"; worst=0; why=""
    if feats:
        L=[(f[:2],"CDS_"+str(f[3])) for f in feats["cds"]]+[(f[:2],"tRNA_"+str(f[3])) for f in feats["trna"] if not f[4]]
    else: L=[((x[0],x[1]),"tRNAcall_"+x[4]) for x in (d.get("trna") or {}).get("aragorn",[]) if x[4]!="GTA" and x[2]==st]+[((x[0],x[1]),"PCG_"+x[3]) for x in (d.get("pcg") or [])]
    for iv,nm in L:
        o=ov((a,b),iv)
        if o>worst: worst,why=o,nm
    yy=[f for f in (feats or {}).get("trna",[]) if f[4]]; dis=bool(yy) and not any(ov((a,b),f[:2])>0.5*min(b-a,f[1]-f[0]) and f[2]==st for f in yy)
    return worst,why,dis,bool(yy)
def n33(d,feats):
    """returns dict(status, N33, reason, pass0, pass5, pass10). status: resolved / stem_audit_review / unresolved"""
    t=d.get("trny"); s=d.get("second") or {}; o=dict(status="unresolved",N33="",reason="",pass0="",pass5="",pass10="",gate=False,score="",stem="",loop="",n_agree="")
    if not t: o["reason"]="trny_module_error"; return o
    if t.get("status")!="called": o["reason"]="not_called:"+str(t.get("reason") or t.get("status")); return o
    sc=float(t["score"]); ag=sum(TD.agree(t,s)); stem=int(t["stem_pairs"]); o.update(N33=t["N33"],score=sc,stem=stem,loop=t["loop"],n_agree=ag)
    nonU=t["N33"]!="U"; base=(sc>=10 and ag>=1) if nonU else (sc>=20 or (sc>=10 and ag>=1))
    if not base: o["reason"]="gate_fail"; return o
    w,why,dis,has=locus(d,feats); o["gate"]=(stem>=4) if nonU else True
    for tol in(0,5,10): o[f"pass{tol}"]=bool(w<=tol and not dis and stem>=4)
    lo=w<=10 and not dis
    if dis: o["reason"]="wrong_locus:disagrees_with_annotated_trnY"
    elif w>10: o["reason"]=f"wrong_locus:overlaps_{why}_{w}nt"
    elif stem<4: o["status"]="stem_audit_review"; o["reason"]="stem_lt4_only"
    if lo and stem>=4: o["status"]="resolved"; o["reason"]=f"overlap_{w}nt"+("" if has else ";no_annotated_trnY")
    return o
def main(od):
    FE={}
    for l in open("data/sweep/features.jsonl"): x=json.loads(l); FE[x["acc"]]=x["f"]
    R={json.loads(l)["acc"]:json.loads(l) for l in open(od+"/out.jsonl")}; M={r["acc"]:r for r in csv.DictReader(open("results/gate3/sweep_audit.tsv"),delimiter="\t")}
    rows=[]
    for a,m in M.items():
        if m["status"]!="included": rows.append(dict(m,uaa="excluded",uaa_note=m["status"],n33="unresolved",base="",n33_reason=m["status"])); continue
        d=R.get(a)
        if d is None: rows.append(dict(m,uaa="not_evaluated",uaa_note="missing_from_run",n33="unresolved",base="",n33_reason="missing_from_run")); continue
        c,note=uaa_class(d); z=n33(d,FE.get(a)); u=d.get("uaa") or {}
        rows.append(dict(m,uaa=c,uaa_note=note,n33=z["status"],base=z["N33"],n33_reason=z["reason"],pass0=z["pass0"],pass5=z["pass5"],pass10=z["pass10"],score=z["score"],stem=z["stem"],loop=z["loop"],n_agree=z["n_agree"],taa_int=u.get("n_taa_internal",""),taa_genes=u.get("n_taa_genes",""),taa_tyr=u.get("n_taa_tyrcol",""),llr=u.get("llr_Y",""),tag_int=u.get("n_tag_internal",""),aligned=u.get("aligned_codons",""),genes=u.get("n_genes_hit",""),codetta=(d.get("codetta") or {}).get("TAA",""),module_errors=";".join(d["err"])))
    K=list(dict.fromkeys(k for r in rows for k in r)); w=csv.DictWriter(open("results/mitogenome_master.tsv","w"),fieldnames=K,delimiter="\t",restval=""); w.writeheader(); w.writerows(rows); return rows
def species(rows,conf=True): return [r for r in rows if r["status"]=="included" and r["provenance"]!="named_extra" and (r["role"]=="confirmatory" or not conf)]
def fet(a,b,c,d):
    from scipy.stats import fisher_exact
    return fisher_exact([[a,b],[c,d]])[1]
def tests(rows,out):
    W=lambda *x: print(*x,file=out); C=species(rows); par=origins.load_nodes()
    W("included genomes",sum(r["status"]=="included" for r in rows),"| confirmatory",len(C)); W("UAA classes (confirmatory)",dict(collections.Counter(r["uaa"] for r in C))); W("UAA classes (all roles)",{ro:dict(collections.Counter(r["uaa"] for r in rows if r["status"]=="included" and r["role"]==ro)) for ro in sorted({r["role"] for r in rows if r["status"]=="included"})})
    W("N33 status (confirmatory)",dict(collections.Counter(r["n33"] for r in C)))
    for tol in(0,5,10): W(f"locus audit pass at {tol} nt (confirmatory):",sum(str(r.get(f"pass{tol}"))=="True" for r in C))
    W("\n== phylum-stratified trnY coverage (confirmatory)"); ph=collections.defaultdict(collections.Counter)
    for r in C: ph[r["phylum"] or "(none)"][r["n33"]]+=1; ph[r["phylum"] or "(none)"]["wrong_locus"]+=r["n33_reason"].startswith("wrong_locus")
    for p,c in sorted(ph.items(),key=lambda x:-sum(v for k,v in x[1].items() if k!="wrong_locus")): n=c["resolved"]+c["unresolved"]+c["stem_audit_review"]; W(f"{p}\tn={n}\tresolved {c['resolved']} ({100*c['resolved']/n:.1f}%)\tstem_audit_review {c['stem_audit_review']}\tunresolved {c['unresolved']} (wrong-locus {c['wrong_locus']})")
    res=[r for r in C if r["n33"]=="resolved"]; st=[r for r in res if r["uaa"]=="Stop"]; ty=[r for r in res if r["uaa"]=="Tyr"]; nu=[r for r in st if r["base"]!="U"]
    states={r["taxid"]:(1 if r["uaa"]=="Tyr" else 0) for r in C if r["uaa"] in("Tyr","Stop")}; og=origins.count(states,par) if states else {}
    W("\n== new Tyr (confirmatory, primary):",sum(r["uaa"]=="Tyr" for r in C),"| origins on confirmatory genomes",og); fam=collections.defaultdict(set)
    for r in C:
        if r["uaa"]=="Tyr": fam[(r["phylum"],r["order"],r["family"] or r["genus"])].add(r["species"]); W("  Tyr",r["acc"],r["species"],r["phylum"],r["order"],r["family"],r["provenance"],"N33",r["n33"],r["base"],"TAA",r["taa_int"],"genes",r["taa_genes"],"tyr",r["taa_tyr"],"llr",r["llr"],"internal TAG",r["tag_int"])
    for k,v in fam.items(): W("  clade",k,"species",len(v),"->","NEW CLADE (two species)" if len(v)>=2 else "single-genome: needs independent evidence (RNA-seq/TSA, proteomics, or separate BioProject)")
    W("\n== T0",{(r["uaa"],r["n33"],r["base"]):0 for r in []}); W("  Tyr by N33 (confirmatory):",dict(collections.Counter((r["n33"],r["base"]) for r in C if r["uaa"]=="Tyr")))
    lo,hi=wil(len(nu),len(st)); W(f"\n== T1 Stop species resolved {len(st)}; non-U33 {len(nu)} = {100*len(nu)/max(1,len(st)):.3f}% (95% CI {100*(lo or 0):.3f}-{100*(hi or 0):.3f})",dict(collections.Counter(r["base"] for r in nu)))
    for p,c in collections.Counter(r["phylum"] for r in nu).most_common(): n=sum(r["phylum"]==p for r in st); W(f"   {p}: {c}/{n}")
    W("   by provenance:",{pv:f"{sum(r['provenance']==pv for r in nu)}/{sum(r['provenance']==pv for r in st)}" for pv in sorted({r['provenance'] for r in st})})
    for tol in(0,5,10): x=[r for r in C if r["uaa"]=="Stop" and str(r.get(f"pass{tol}"))=="True"]; W(f"   at {tol}-nt tolerance: non-U33 {sum(r['base']!='U' for r in x)}/{len(x)}")
    g=og.get("gains_min",0); W("   Fisher (Tyr x non-U33):",f"p={fet(sum(r['base']!='U' for r in ty),sum(r['base']=='U' for r in ty),len(nu),len(st)-len(nu)):.3g}" if g>=3 else f"not run (new independent origins {g} < 3)")
    W("   independent origins of non-U33 among Stop (resolved):",origins.count({r["taxid"]:int(r["base"]!="U") for r in st},par) if st else {})
    W("\n== sufficiency counterexample candidates (confirmatory, Stop, resolved, non-U33)"); [W("  ",r["acc"],r["species"],r["phylum"],r["class"],r["order"],r["family"],r["provenance"],r["base"],r["loop"],"score",r["score"],"stem",r["stem"],"agree",r["n_agree"],r["n33_reason"]) for r in sorted(nu,key=lambda r:(r["phylum"],r["order"],r["species"]))]
    tp=sum(r["base"]!="U" for r in ty); pp=tp+len(nu); W("\n== T2 precision",f"{tp}/{pp}","recall",f"{tp}/{len(ty)}","| inferential" if (len(ty)>=5 and g>=3) else f"| descriptive only (new Tyr species {len(ty)}, origins {g}; needs >=5 and >=3)")
    inc=[r for r in rows if r["status"]=="included"]
    E=[r for r in inc if r["phylum"]=="Echinodermata"]; W("\n== T4 Echinodermata by class: (N33 status, base, UAA)")
    for c in sorted({r["class"] for r in E}): W("  ",c or "(none)",dict(collections.Counter((r["n33"],r["base"] if r["n33"]=="resolved" else "",r["uaa"]) for r in E if r["class"]==c)))
    AS=[r for r in E if r["class"]=="Asteroidea" and r["n33"]=="resolved" and r["provenance"]!="named_extra"]; c33=[r for r in AS if r["base"]=="C" and r["uaa"]=="Stop"]
    rd="inconclusive" if len(AS)<2 else "Li-consistent" if (len(c33)>=2 and len({r["family"] for r in c33})>=2) else "Zhao-consistent" if (len(AS)>=5 and sum(r["base"]=="U" for r in AS)/len(AS)>=.8) else "mixed"
    both=(len(c33)>=2 and len({r["family"] for r in c33})>=2) and (len(AS)>=5 and sum(r["base"]=="U" for r in AS)/len(AS)>=.8)
    W("  asteroid resolved",len(AS),dict(collections.Counter(r["base"] for r in AS)),"READING:",rd,"(BOTH criteria met)" if both else ""); [W("   Asterias rubens:",r["acc"],r["n33"],r["base"],r["loop"],r["uaa"],r["n33_reason"]) for r in inc if r["species"]=="Asterias rubens"]
    EL=[r for r in inc if r["order"] in("Albuliformes","Anguilliformes","Elopiformes","Notacanthiformes")]; W("\n== T6 Elopomorpha")
    for o in sorted({r["order"] for r in EL}): W("  ",o,dict(collections.Counter((r["n33"],r["base"] if r["n33"]=="resolved" else "",r["uaa"]) for r in EL if r["order"]==o)))
    [W("   Albuliformes:",r["acc"],r["species"],r["role"],r["n33"],r["base"],r["loop"],r["score"],r["stem"],r["n_agree"],r["uaa"],r["n33_reason"]) for r in EL if r["order"]=="Albuliformes"]
    oth=[r for r in EL if r["order"]=="Albuliformes" and r["n33"]=="resolved" and r["acc"]!="NC_082992.1"]; W("   READING:","lineage-level" if any(r["base"]=="C" for r in oth) else "singleton (NUMT/read checks pending)" if oth and all(r["base"]=="U" for r in oth) else "inconclusive")
    W("\n== weak tier (confirmatory):",sum(r["uaa"]=="weak" for r in C)); [W("  ",r["acc"],r["species"],r["phylum"],r["order"],"TAA",r["taa_int"],"genes",r["taa_genes"],"tyr",r["taa_tyr"],"codetta",r["codetta"],r["uaa_note"],"N33",r["n33"],r["base"]) for r in C if r["uaa"]=="weak"]
    W("\n== stem_audit_review (all included):",sum(r["n33"]=="stem_audit_review" for r in inc),"non-U33 among them:",[(r["acc"],r["species"],r["base"],r["uaa"]) for r in inc if r["n33"]=="stem_audit_review" and r["base"]!="U"])
    W("unresolved reasons (confirmatory):",collections.Counter(r["n33_reason"].split(":")[0]+(":"+r["n33_reason"].split(":")[1].split("_")[0] if ":" in r["n33_reason"] else "") for r in C if r["n33"]=="unresolved").most_common(12)); W("module errors:",sum(bool(r.get("module_errors")) for r in inc))
if __name__=="__main__":
    rows=main(sys.argv[1]); 
    with open("results/gate3/tests.txt","w") as o: tests(rows,o)
    print(open("results/gate3/tests.txt").read()[:3000])
