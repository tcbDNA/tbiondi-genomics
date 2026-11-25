"""Evaluate the armless nematode models. Nematode recovery uses only the model that EXCLUDED the genome's order; false non-U33 is counted for every model on all negatives."""
import json,csv,sys,os,collections,re
sys.path.insert(0,"src"); import trny_decide as TD, trny_diag_eval as E
P,S,AN,known,neg=E.P,E.S,E.AN,E.known,E.neg
F=json.load(open("results/gate2/trny_diag/nem_folds.json")); order=F["order"]
A=TD.load("results/gate2/trny_diag/A.jsonl"); M={}
for o in F["folds"]:
    p=f"results/gate2/trny_diag/C_nem_minus_{re.sub(chr(92)+'W','',o)}.jsonl"
    if os.path.exists(p):
        x=TD.load(p)
        if len(x)==360: M[o]=x
p="results/gate2/trny_diag/C_nem_all.jsonl"; ALL=TD.load(p) if os.path.exists(p) else {}
print("models complete:",list(M),"nem_all",len(ALL)); nem=[k for k in known if P[k]["lineage"]=="Nematoda"]
def tophit(T,k): return T[k]["status"]=="called"
for thr in (20,25,30,35):
    rec=0;wrong=0;n=0
    for k in nem:
        o=order.get(k)
        if o not in M: continue
        n+=1; d=M[o][k]
        if tophit(d and M[o],k) and float(d["score"])>=thr:
            if E.disc(k,d): wrong+=1
            else: rec+=1
    print(f"held-out nematodes (n={n}) score>={thr}: right locus {rec}, wrong locus {wrong}, not called {n-rec-wrong}")
base_fail=[k for k in nem if not E.gate(A[k],S[k])]
for thr in (20,25,30,35):
    g=[(k,M[order[k]][k]) for k in base_fail if order.get(k) in M and M[order[k]][k]["status"]=="called" and float(M[order[k]][k]["score"])>=thr]
    print(f"  of {len(base_fail)} nematodes failing variant A: recovered at thr {thr}: right {sum(not E.disc(k,d) for k,d in g)} wrong {sum(bool(E.disc(k,d)) for k,d in g)}")
print("\nfalse non-U33 on negatives (model hit alone, no corroboration):")
for nm,T in list(M.items())+([("nem_all",ALL)] if len(ALL)==360 else []):
    print("  ",nm,{thr:sum(1 for k in neg if T[k]["status"]=="called" and T[k]["N33"]!="U" and float(T[k]["score"])>=thr) for thr in(15,20,25,30)},"| nematodes only",{thr:sum(1 for k in neg if P[k]["lineage"]=="Nematoda" and T[k]["status"]=="called" and T[k]["N33"]!="U" and float(T[k]["score"])>=thr) for thr in(15,20,25,30)},"| wrong-locus accepted>=20 all neg",sum(1 for k in neg if T[k]["status"]=="called" and float(T[k]["score"])>=20 and E.disc(k,T[k])))
