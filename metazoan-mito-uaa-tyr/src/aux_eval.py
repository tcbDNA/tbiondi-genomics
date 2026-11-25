import json,csv,collections
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
A=json.load(open("results/gate2/annot_terminals.json")); X={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/aux_panel.jsonl")}
print("errors",[(k,d["err"]) for k,d in X.items() if "err" in d][:3])
G={"unknown","unaligned_3prime","none"}
inc=lambda t:"inc" if t in("T--","TA-") else t
for role in("neg_matched","neg_broad","pos"):
    cm=collections.Counter(); cov=collections.Counter()
    for k,d in X.items():
        if P[k]["role"]!=role or k not in A or "terminal" not in d: continue
        for g,t in A[k].items():
            if t=="none": continue
            m=d["terminal"].get(g,"gene_not_found"); cov["total"]+=1
            if m in G or m=="gene_not_found": cov[m]+=1
            else: cm[(inc(t),inc(m))]+=1
    n=sum(cm.values()); c=sum(v for (a,b),v in cm.items() if a==b)
    print(f"{role}: annotated genes {cov['total']}; called {n} ({100*n/max(1,cov['total']):.0f}%); concordant {c}/{n}={100*c/max(1,n):.1f}% | uncalled {dict((a,b) for a,b in cov.items() if a!='total')}"); print("   (annotated,module):",{k:v for k,v in cm.most_common()})
print("\nper-lineage, positives and matched: module terminal counts vs annotation (called genes)")
for ln in sorted({P[k]["lineage"] for k in X if P[k]["role"]!="neg_broad"}):
    m=collections.Counter(); a=collections.Counter()
    for k,d in X.items():
        if P[k]["lineage"]==ln and P[k]["role"]!="neg_broad":
            for g,t in d.get("terminal",{}).items():
                if t not in G: m[inc(t)]+=1
            for g,t in A.get(k,{}).items(): a[inc(t)]+=1
    print("  ",ln,"module",dict(m),"| annotation",dict(a))
print("\ncodon O/E: A. chinensis (Zhao re-assembly) TAA/TAT/TAC",[X["ZHAO_SRR35647597__Anomia_chinensis"]["oe"][x] for x in("TAA","TAT","TAC")],"Zhao Table S4 mean counts 5.4/3.5/7 -> fractions 0.34/0.22/0.44; mine",end=" ")
o=X["ZHAO_SRR35647597__Anomia_chinensis"]["oe"]; s=o["TAA"]+o["TAT"]+o["TAC"]; print([round(o[x]/s,2) for x in("TAA","TAT","TAC")])
for k in X:
    if P[k]["role"]=="pos": o=X[k]["oe"]; print("  ",k,"obs TAA frac",o["obs_taa_frac"],"exp(4-fold)",o["exp_taa_frac"],"O/E",o["oe"])
