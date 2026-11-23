"""Leave-one-lineage-out and leave-one-origin-out recall of the UAA=Tyr call. Codetta setting is also re-chosen in each fold.
Procedures (all reported): N = negatives only (max negative + 1; LLR>0; never sees a positive); M = midpoint between max negative and min training positive; R = 0.8 x min training positive.
Fixed = the rule in PREREG, which was hand-set after seeing all positives (in-sample)."""
import json,csv,math,collections
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
J=lambda f:{json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/"+f)}
D=J("uaa_panel.jsonl"); CL=J("codetta_mitodb_lin.jsonl"); CD=J("codetta_mitodb_default.jsonl"); PL=[json.loads(l) for l in open("results/gate2/planted_panel.jsonl")]
F=("n_taa_internal","n_taa_genes","n_taa_tyrcol","llr_Y")
pos=[k for k in P if P[k]["role"]=="pos"]; neg=[k for k in P if P[k]["role"]!="pos"]
DUP={"ZHAO_SRR35647597__Anomia_chinensis":"LIN_Anomia_chinensis","LIN_Anomia_trigonopsis_NZ":"LIN_Anomia_trigonopsis_Mauritania","LIN_Enigmonia_aenigmatica_S_yellow":"LIN_Enigmonia_aenigmatica_S_purple","ZHAO_SRR35651820_Placuna_placenta":"LIN_Placuna_vitream","ZHAO_ERR15194044_Pseudotanais_sp_a":"OZ281199.1"}
ORI={"Anomia":"Anomioidea","Enigmonia":"Anomioidea","Placunidae":"Anomioidea","Plicatulidae":"Plicatulidae","Cephalodiscidae":"Cephalodiscidae","Radopholus":"Radopholus","Pseudotanais":"Pseudotanais"}
mxn=[max(D[k][f] for k in neg) for f in F]
def thr(proc,tr):
    mn=[min(D[k][f] for k in tr) for f in F]
    if proc=="N": return [mxn[0]+1,mxn[1]+1,mxn[2]+1,1e-9]
    if proc=="M": return [(a+b)/2 for a,b in zip(mxn,mn)]
    if proc=="R": return [0.8*b for b in mn]
    return [5,3,3,3]
cod=lambda C,k: C[k].get("TAA")=="Y"
def codset(tr):
    best=None
    for nm,C in(("default",CD),("lin",CL)):
        if any(cod(C,k) for k in neg): continue
        r=sum(cod(C,k) for k in tr)/len(tr)
        if best is None or r>best[0]: best=(r,nm,C)
    return best
call=lambda k,t,C: all(D[k][f]>=x for f,x in zip(F,t)) and cod(C,k)
out=[]
for unit,key in(("lineage",lambda k:P[k]["lineage"]),("origin",lambda k:ORI[P[k]["lineage"]])):
    for proc in("Fixed","N","M","R"):
        for u in sorted({key(k) for k in pos}):
            te=[k for k in pos if key(k)==u]; tr=[k for k in pos if key(k)!=u]; t=thr(proc,tr)
            r=codset(tr) if proc!="Fixed" else (None,"lin",CL); C=r[2]
            td=[k for k in te if k not in DUP]; fp=sum(call(k,t,C) for k in neg)
            out.append(dict(unit=unit,proc=proc,held_out=u,genomes=len(te),recalled=sum(call(k,t,C) for k in te),counter_only=sum(all(D[k][f]>=x for f,x in zip(F,t)) for k in te),distinct=len(td),distinct_recalled=sum(call(k,t,C) for k in td),codetta=r[1],thresholds=[round(x,2) for x in t],false_pos_340=fp))
w=csv.DictWriter(open("results/gate2/lolo.tsv","w"),fieldnames=list(out[0]),delimiter="\t"); w.writeheader(); w.writerows(out)
for o in out: print(o["unit"],o["proc"],o["held_out"],f"{o['recalled']}/{o['genomes']}",f"counter-only {o['counter_only']}",f"distinct {o['distinct_recalled']}/{o['distinct']}",o["codetta"],o["thresholds"],"FP",o["false_pos_340"])
print("max over negatives",mxn)
def weak(d,codY): return (d["n_taa_internal"]>=2 and d["n_taa_tyrcol"]>=1) or codY
print("weak tier on negatives:",[(k,P[k]["species"]) for k in neg if weak(D[k],cod(CL,k)) and not call(k,[5,3,3,3],CL)])
pc={(json.loads(l)["base"],json.loads(l)["frac"]):json.loads(l).get("TAA") for l in open("results/gate2/planted_codetta.jsonl")}
for f in sorted({d["frac"] for d in PL}):
    x=[d for d in PL if d["frac"]==f and d["rep"]==0]; pr=[all(d[a]>=b for a,b in zip(F,[5,3,3,3])) and (f>0 and pc.get((d["base"],f))=="Y") for d in x]; wk=[(not p) and weak(d,f>0 and pc.get((d["base"],f))=="Y") for d,p in zip(x,pr)]
    print("planted",f,"n",len(x),"primary",sum(pr),"weak",sum(wk),"neither",len(x)-sum(pr)-sum(wk))
