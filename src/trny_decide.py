"""Combine CM call + ARAGORN + tRNAscan-SE -O (GUA anticodon) into a final trnY status. Final rule (locked in PREREG):
 called_high : CM score>=20 AND N33 known  OR  (score>=10 AND >=1 independent caller agrees at the same locus)
 non-U33 calls additionally require (a) >=1 independent caller agreeing AND (b) >=4/5 stem pairs.
 Everything else -> no_confident_call (reason)."""
import json,csv,collections,sys
def load(p): return {json.loads(l)["panel_id"]:json.loads(l) for l in open(p)}
def agree(d,s):
    if d["status"]!="called" or not d.get("g34_pos"): return (False,False)
    a,b=sorted((int(d["start"]),int(d["end"])))
    ar=any(abs(x[0]-d["g34_pos"])<=2 and x[1]==("+" if d["strand"]!="-" else "-") for x in s.get("aragorn",[]))
    ts=any(min(b,max(p,q))-max(a,min(p,q))>0.5*(b-a) for p,q,sc in s.get("trnascan",[]))
    return (ar,ts)
def decide(d,s):
    ag=agree(d,s); n=sum(ag); sc=float(d.get("score") or 0)
    if d["status"]!="called": return dict(final="no_confident_call",reason=d.get("reason","") or d["status"],n_agree=0,N33="")
    good=(sc>=20) or (sc>=10 and n>=1)
    if not good: return dict(final="no_confident_call",reason="score_below_gate_no_corroboration",n_agree=n,N33=d["N33"],score=sc)
    if d["N33"]!="U" and not (n>=1 and d["stem_pairs"]>=4): return dict(final="no_confident_call",reason="nonU33_uncorroborated",n_agree=n,N33=d["N33"],score=sc)
    return dict(final="called",reason=d["reason"],n_agree=n,N33=d["N33"],score=sc,loop=d["loop"],stem=d["stem_pairs"])
if __name__=="__main__":
    T=load(sys.argv[1]); S=load(sys.argv[2]); P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
    out=csv.writer(open(sys.argv[3],"w"),delimiter="\t"); out.writerow("panel_id species lineage role final reason n_agree N33 score".split()); R={}
    for k,d in T.items():
        s=S.get(k,{}); r=decide(d,s); R[k]=r; out.writerow([k,P[k]["species"],P[k]["lineage"],P[k]["role"],r["final"],r["reason"],r["n_agree"],r["N33"],r.get("score","")])
    for role in ("pos","neg_matched","neg_broad"):
        x={k:r for k,r in R.items() if P[k]["role"]==role}
        print(role,len(x),collections.Counter(r["final"] for r in x.values()),"N33",collections.Counter(r["N33"] for r in x.values() if r["final"]=="called"))
    nonV=[k for k in R if P[k]["role"]!="pos" and P[k]["lineage"]!="Chordata"]; nc=[k for k in nonV if R[k]["final"]!="called"]; print("non-vertebrate(non-Chordata) negatives",len(nonV),"no_confident_call",len(nc),round(100*len(nc)/len(nonV),1),"%")
    print("  by lineage",collections.Counter(P[k]["lineage"] for k in nc).most_common()); print("  reasons",collections.Counter(R[k]["reason"] for k in nc))
    print("non-U33 called:",[(k,P[k]["species"],P[k]["role"],R[k]["N33"],R[k]["score"]) for k in R if R[k]["final"]=="called" and R[k]["N33"]!="U" and P[k]["role"]!="pos"])
    print("positives not called:",[(k,R[k]["reason"],R[k]["n_agree"]) for k in R if P[k]["role"]=="pos" and R[k]["final"]!="called"])
