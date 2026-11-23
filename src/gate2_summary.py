"""Reproduce every number quoted in GATE_2_CALIBRATION.md from the work files (rules as proposed for PREREG)."""
import json,csv,collections,sys
sys.path.insert(0,"src"); import trny_decide as TD
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
D={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/work/uaa_panel.jsonl")}
CO={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/work/codetta_mitodb_lin.jsonl")}
CD={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/work/codetta_mitodb_default.jsonl")}
PL=[json.loads(l) for l in open("results/work/planted_panel.jsonl")]
C=dict(min_codons=800,min_genes_hit=5,min_tyrcols=15,ti=5,tg=3,ty=3,llr=3,stop_max=1,stop_tyr=0)
def counter(d,c=C):
    if d["aligned_codons"]<c["min_codons"] or d["n_genes_hit"]<c["min_genes_hit"] or d["tyrcols"]<c["min_tyrcols"]: return "insufficient"
    if d["n_taa_internal"]>=c["ti"] and d["n_taa_genes"]>=c["tg"] and d["n_taa_tyrcol"]>=c["ty"] and d["llr_Y"]>=c["llr"]: return "Tyr"
    if d["n_taa_internal"]<=c["stop_max"] and d["n_taa_tyrcol"]<=c["stop_tyr"]: return "Stop"
    return "ambiguous"
def final(d,co):
    x=counter(d)
    if x=="Tyr": return "Tyr" if co=="Y" else "ambiguous"
    if x=="Stop": return "Stop" if co!="Y" else "ambiguous"
    return x
tab=collections.defaultdict(collections.Counter)
for k in D: tab[P[k]["role"]][final(D[k],CO[k].get("TAA"))]+=1
print("UAA final",{r:dict(v) for r,v in tab.items()}); print("codetta Lin TAA calls",{r:dict(collections.Counter(CO[k].get("TAA") for k in D if P[k]["role"]==r)) for r in tab})
print("codetta default TAA calls",{r:dict(collections.Counter(CD[k].get("TAA") for k in D if P[k]["role"]==r)) for r in tab})
pc=[json.loads(l) for l in open("results/work/planted_codetta.jsonl")]; pc={(d["base"],d["frac"]):d.get("TAA") for d in pc}
res=collections.defaultdict(collections.Counter); cnt=collections.Counter()
for d in PL:
    res[d["frac"]][counter(d)]+=1
    if d["rep"]==0: cnt[(d["frac"],final(d,pc.get((d["base"],d["frac"]),"?") if d["frac"]>0 else "?"))]+=1
print("planted counter-only (all reps)",{f:dict(v) for f,v in sorted(res.items())}); print("planted final rep0 (counter & codetta)",dict(cnt))
T=TD.load("results/work/trny_panel.jsonl"); S=TD.load("results/work/trny_second.jsonl")
def tr(k,lo=10,hi=20):
    d,s=T[k],S[k]
    if d["status"]!="called": return ("no_call",d.get("reason"),"")
    sc=float(d["score"]); ag=sum(TD.agree(d,s))
    if d["N33"]!="U": return ("called","corrob",d["N33"]) if (sc>=lo and ag>=1 and d["stem_pairs"]>=4) else ("no_call","nonU33_uncorroborated",d["N33"])
    return ("called","std",d["N33"]) if (sc>=hi or (sc>=lo and ag>=1)) else ("no_call","score_low_uncorroborated",d["N33"])
for role in("pos","neg_matched","neg_broad"):
    x=[tr(k) for k in P if P[k]["role"]==role]; print("trnY",role,len(x),collections.Counter(a for a,_,_ in x),"N33",collections.Counter(c for a,_,c in x if a=="called"),"reasons",collections.Counter(b for a,b,_ in x if a!="called"))
print("pos detail:",[(k,tr(k)) for k in P if P[k]["role"]=="pos" and tr(k)[0]!="called"])
ann=json.load(open("results/work/annotated_trnY.json"))
neg=[k for k in P if P[k]["role"]!="pos" and P[k]["lineage"]!="Chordata"]; known=[k for k in neg if ann.get(k,0)>0 or P[k]["role"]=="neg_matched"]
for nm,s in (("all negatives non-Chordata",neg),("with annotated/known trnY",known),("with annotated trnY excl Cnidaria",[k for k in known if P[k]["lineage"]!="Cnidaria"])):
    f=[k for k in s if tr(k)[0]!="called"]; print(nm,len(s),"fail",len(f),f"{100*len(f)/len(s):.1f}%",collections.Counter(P[k]["lineage"] for k in f).most_common(6))
json.dump({k:tr(k) for k in P},open("results/work/trny_final_panel.json","w"))
