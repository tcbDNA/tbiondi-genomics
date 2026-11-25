import json,csv,collections
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
A=json.load(open("results/gate2/annot_terminals.json")); X={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/aux_panel.jsonl")}
print("Error of each symmetric rule on genes with an annotated COMPLETE stop (negatives only; annotation independent of the rule):")
for m in("first","tagfirst"):
    c=collections.Counter()
    for k,d in X.items():
        if P[k]["role"]=="pos": continue
        for g,t in A.get(k,{}).items():
            v=d[m].get(g)
            if t in("TAA","TAG") and v in("TAA","TAG"): c[(t,v)]+=1
    print("  ",m,dict(c),"TAA called TAG:",f"{100*c[('TAA','TAG')]/max(1,c[('TAA','TAG')]+c[('TAA','TAA')]):.1f}%")
ORI={"Anomia":"Anomioidea(Tyr)","Enigmonia":"Anomioidea(Tyr)","Placunidae":"Anomioidea(Tyr)"}
print("\nTAG fraction among TAA/TAG calls, per lineage (genes; near-duplicate genomes included):")
rows=collections.defaultdict(lambda: collections.Counter())
for k,d in X.items():
    if P[k]["role"]=="neg_broad": ln="BROAD non-Chordata" if P[k]["lineage"]!="Chordata" else "BROAD Chordata"
    else: ln=ORI.get(P[k]["lineage"],P[k]["lineage"]+("(Tyr)" if P[k]["role"]=="pos" else "(Stop)"))
    for m in("first","tagfirst"):
        for v in d[m].values():
            if v in("TAA","TAG"): rows[ln][(m,v)]+=1
for ln in sorted(rows):
    c=rows[ln]; f=lambda m: f"{c[(m,'TAG')]}/{c[(m,'TAG')]+c[(m,'TAA')]}={100*c[(m,'TAG')]/max(1,c[(m,'TAG')]+c[(m,'TAA')]):.0f}%"; print("  ",ln,"first",f("first"),"| tagfirst",f("tagfirst"))
