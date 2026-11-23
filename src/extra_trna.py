"""Extra-tRNA check: tRNAs whose anticodon could read UAA/UAG/UAU (TTA, CTA, ATA) and duplicate GUA copies, per genome. Panel only.
corroborated = ARAGORN and tRNAscan-SE agree on locus and anticodon; single = one caller only."""
import json,csv,collections,sys
sys.path.insert(0,"src"); import panel_guard
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; panel_guard.check(P)
T={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/trna_all.jsonl")}
def ov(a,b): return min(a[1],b[1])-max(a[0],b[0])>0.5*min(a[1]-a[0],b[1]-b[0]) and a[2]==b[2]
def scan(t):
    A=[x for x in t["aragorn"]]; S=[x for x in t["trnascan"]]; o={}
    for ac in("TTA","CTA","ATA","GTA"):
        a=[x for x in A if x[4]==ac]; s=[x for x in S if x[4]==ac]
        cor=sum(any(ov(x,y) for y in s) for x in a); o[ac]=dict(corroborated=cor,single=len(a)+len(s)-2*cor,hits=[("A",x[0],x[1],x[2],x[3]) for x in a]+[("S",y[0],y[1],y[2],y[3],y[5]) for y in s])
    return o
R={k:scan(t) for k,t in T.items()}
w=csv.writer(open("results/gate2/extra_trna.tsv","w"),delimiter="\t"); w.writerow("panel_id species lineage role TTA_cor TTA_single CTA_cor CTA_single ATA_cor ATA_single GTA_cor GTA_single".split())
for k,o in R.items(): w.writerow([k,P[k]["species"],P[k]["lineage"],P[k]["role"]]+[o[a][b] for a in("TTA","CTA","ATA","GTA") for b in("corroborated","single")])
for ac in("TTA","CTA","ATA"):
    for role in("pos","neg_matched","neg_broad"):
        ks=[k for k in R if P[k]["role"]==role]; print(ac,role,"n",len(ks),"corroborated",sum(R[k][ac]["corroborated"]>0 for k in ks),"single-caller only",sum(R[k][ac]["corroborated"]==0 and R[k][ac]["single"]>0 for k in ks))
print("GUA copies >1 corroborated:",collections.Counter((P[k]["role"],R[k]["GTA"]["corroborated"]>1) for k in R))
print("\npositives with TTA/CTA/ATA hits:"); 
for k in R:
    if P[k]["role"]=="pos":
        h={a:R[k][a]["hits"] for a in("TTA","CTA","ATA") if R[k][a]["hits"]}
        if h: print(k,h)
print("errors",[k for k in T if T[k]["err"]])
