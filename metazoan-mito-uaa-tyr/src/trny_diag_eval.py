"""Evaluate extractor variants on the panel with the locked gate (U33: score>=20 or >=10+corroboration; non-U33: >=10 + corroboration + stem>=4)."""
import json,csv,sys,collections,os
sys.path.insert(0,"src"); import trny_decide as TD
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
ann=json.load(open("results/gate2/annotated_trnY.json")); AN=json.load(open("results/gate2/trny_diag/annotated_trny_panel.json"))
S=TD.load("results/gate2/trny_second.jsonl")
neg=[k for k in P if P[k]["role"]!="pos"]; nv=[k for k in neg if P[k]["lineage"]!="Chordata"]; known=[k for k in nv if ann.get(k,0)>0 or P[k]["role"]=="neg_matched"]
def gate(d,s):
    if d["status"]!="called": return False
    sc=float(d["score"]); ag=sum(TD.agree(d,s))
    if d["N33"]!="U": return sc>=10 and ag>=1 and d["stem_pairs"]>=4
    return sc>=20 or (sc>=10 and ag>=1)
def disc(k,d):
    A=AN.get(k) or []
    if not A: return None
    a,b=sorted((int(d["start"]),int(d["end"]))); return not any(min(b,x["end"])-max(a,x["start"])>0.5*min(b-a,x["end"]-x["start"]) for x in A)
def report(nm,T):
    ok={k:gate(T[k],S[k]) for k in P if k in T}
    f=[k for k in known if not ok.get(k)]; fa=[k for k in nv if not ok.get(k)]
    raw=[k for k in neg if T[k]["status"]=="called" and T[k]["N33"]!="U"]; acc=[k for k in neg if ok.get(k) and T[k]["N33"]!="U"]
    dd=[disc(k,T[k]) for k in nv if ok.get(k)]; dd=[x for x in dd if x is not None]
    nem=[k for k in known if P[k]["lineage"]=="Nematoda"]; pos=[k for k in P if P[k]["role"]=="pos"]
    print(f"{nm}: fail known-gene {len(f)}/{len(known)}={100*len(f)/len(known):.1f}% | all {len(fa)}/{len(nv)} | Nematoda {sum(not ok[k] for k in nem)}/{len(nem)} | accepted non-U33 on negatives {[(k,P[k]['species'],T[k]['N33']) for k in acc]} | raw non-U33 {len(raw)} | discordant {sum(dd)}/{len(dd)} | no-gene accepted {sum(ok[k] for k in nv if k not in known)}/{len(nv)-len(known)} | positives called {sum(ok[k] for k in pos)}/20 C33 {sum(ok[k] and T[k]['N33']=='C' for k in pos)} exc {sum(T[k]['status']=='exception' for k in T)}")
    return ok
if __name__=="__main__":
    R={"baseline":report("baseline",TD.load("results/gate2/trny_panel.jsonl"))}
    for v in sys.argv[1:]: R[v]=report(v,TD.load(f"results/gate2/trny_diag/{v}.jsonl"))
    ch=collections.Counter()
    for v in sys.argv[1:]:
        b=R["baseline"]; print(v,"newly recovered:",[(k,P[k]["lineage"]) for k in known if R[v].get(k) and not b[k]],"| newly lost:",[(k,P[k]["lineage"]) for k in known if b[k] and not R[v].get(k)])
