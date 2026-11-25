import json,csv,random,sys,collections; sys.path.insert(0,"src"); import trny_decide as TD
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
ann=json.load(open("results/work/annotated_trnY.json")); T=TD.load("results/work/trny_panel.jsonl"); S=TD.load("results/work/trny_second.jsonl")
con={r["panel_id"]:r for r in csv.DictReader(open("results/work/trny_concordance.tsv"),delimiter="\t")}
def rule(k,lo,hi):
    d,s=T[k],S[k]
    if d["status"]!="called": return False
    sc=float(d["score"]); ag=sum(TD.agree(d,s))
    if d["N33"]!="U": return sc>=lo and ag>=1 and d["stem_pairs"]>=4
    return sc>=hi or (sc>=lo and ag>=1)
A=[k for k in P if P[k]["role"]=="neg_broad" and P[k]["lineage"]!="Chordata" and ann.get(k,0)>0]
random.Random(7).shuffle(A); h1,h2=A[::2],A[1::2]
def ev(ks,lo,hi):
    f=sum(not rule(k,lo,hi) for k in ks)/len(ks); acc=[k for k in ks if rule(k,lo,hi) and con[k]["overlap"] in("yes","no")]; d=sum(con[k]["overlap"]=="no" for k in acc)/max(1,len(acc)); return f,d
G=[(lo,hi) for lo in (8,10,12) for hi in (15,18,20,22,25)]
for name,tr,te in (("train half A, test B",h1,h2),("train B, test A",h2,h1)):
    best=min((g for g in G if ev(tr,*g)[1]<=0.03),key=lambda g:ev(tr,*g)[0],default=None)
    print(name,"chosen (discordance<=3% on train):",best,"train fail/disc",tuple(round(x,3) for x in ev(tr,*best)) if best else None,"TEST fail/disc",tuple(round(x,3) for x in ev(te,*best)) if best else None)
for g in [(10,15),(10,18),(10,20),(10,25)]: print("full annotated set",g,tuple(round(x,3) for x in ev(A,*g)),"n",len(A))
import math
for k,n in ((22,273),(26,273),(28,273)):
    p=k/n; print(k,n,"95% Wilson",[round((p+1.96**2/(2*n)+s*1.96*math.sqrt(p*(1-p)/n+1.96**2/(4*n*n)))/(1+1.96**2/n),3) for s in(-1,1)])
