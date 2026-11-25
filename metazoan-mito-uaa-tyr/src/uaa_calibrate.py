"""Evaluate candidate UAA rules on the panel and planted series; print counts. Rules are compared, the lock happens in PREREG."""
import json,csv,collections
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
D=[json.loads(l) for l in open("results/work/uaa_panel.jsonl")]; PL=[json.loads(l) for l in open("results/work/planted_panel.jsonl")]
def call(d,c):
    if d["aligned_codons"]<c["min_codons"] or d["n_genes_hit"]<c["min_genes_hit"] or d["tyrcols"]<c["min_tyrcols"]: return "insufficient"
    if d["n_taa_internal"]>=c["ti"] and d["n_taa_genes"]>=c["tg"] and d["n_taa_tyrcol"]>=c["ty"] and d["llr_Y"]>=c["llr"]: return "Tyr"
    if d["n_taa_internal"]<=c["stop_max"] and d["n_taa_tyrcol"]<=c["stop_tyr"]: return "Stop"
    return "ambiguous"
base=dict(min_codons=800,min_genes_hit=5,min_tyrcols=15,ti=8,tg=3,ty=3,llr=3,stop_max=1,stop_tyr=0)
for nm,c in (("base",base),("ti5",{**base,"ti":5}),("ti10_tg4",{**base,"ti":10,"tg":4}),("ti8_llr10",{**base,"llr":10})):
    x=collections.defaultdict(collections.Counter)
    for d in D: x[P[d["panel_id"]]["role"]][call(d,c)]+=1
    print(nm,{k:dict(v) for k,v in x.items()})
    for f in sorted({d["frac"] for d in PL}):
        y=[call(d,c) for d in PL if d["frac"]==f]; print("   planted",f,len(y),collections.Counter(y))
print("insufficient in broad:",[(d["panel_id"],P[d["panel_id"]]["lineage"],d["aligned_codons"],d["n_genes_hit"],d["tyrcols"]) for d in D if call(d,base)=="insufficient" and P[d["panel_id"]]["role"]=="neg_broad"])
print("min of positives:",min(d["n_taa_internal"] for d in D if P[d["panel_id"]]["role"]=="pos"),min(d["n_taa_genes"] for d in D if P[d["panel_id"]]["role"]=="pos"),min(d["n_taa_tyrcol"] for d in D if P[d["panel_id"]]["role"]=="pos"),min(d["llr_Y"] for d in D if P[d["panel_id"]]["role"]=="pos"),"max of negatives (int,genes,tyr,llr):",max(d["n_taa_internal"] for d in D if P[d["panel_id"]]["role"]!="pos"),max(d["n_taa_genes"] for d in D if P[d["panel_id"]]["role"]!="pos"),max(d["n_taa_tyrcol"] for d in D if P[d["panel_id"]]["role"]!="pos"),max(d["llr_Y"] for d in D if P[d["panel_id"]]["role"]!="pos"))
