"""Compare the Gate 3 runner on the panel with the Gate 2 outputs, module by module."""
import json,sys,csv,collections
J=lambda p,k="panel_id":{json.loads(l)[k]:json.loads(l) for l in open(p)}
G=J("results/work/gate3_paneltest/out.jsonl","acc"); T=J("results/gate2/trny_panel_final.jsonl"); S=J("results/gate2/trny_second.jsonl"); U=J("results/gate2/uaa_panel.jsonl"); C=J("results/gate2/codetta_mitodb_lin.jsonl"); X=J("results/gate2/aux_panel.jsonl"); A=J("results/gate2/trna_all.jsonl")
n=collections.Counter(); ex=collections.defaultdict(list)
for k,g in G.items():
    chk={"trny":all(g["trny"].get(f)==T[k].get(f) for f in("status","start","end","strand","score","N33","loop","stem_pairs","g34_pos")),
     "second":sorted(map(str,g["second"]["aragorn"]))==sorted(map(str,S[k]["aragorn"])) and sorted(map(str,g["second"]["trnascan"]))==sorted(map(str,S[k]["trnascan"])),
     "trna":g["trna"]["aragorn"]==A[k]["aragorn"] and g["trna"]["trnascan"]==A[k]["trnascan"],
     "uaa":all(g["uaa"][f]==U[k][f] for f in g["uaa"] if f!="panel_id"),
     "codetta":g["codetta"].get("code")==C[k].get("code"),"first":g["first"]==X[k]["first"],"tagfirst":g["tagfirst"]==X[k]["tagfirst"],"oe":g["oe"]==X[k]["oe"]}
    for m,ok in chk.items():
        n[(m,ok)]+=1
        if not ok: ex[m].append(k)
for m in("trny","second","trna","uaa","codetta","first","tagfirst","oe"): print(m,"identical",n[(m,True)],"different",n[(m,False)],ex[m][:4])
print("genomes",len(G),"module errors",sum(bool(g["err"]) for g in G.values()))
s=sorted(g["sec"] for g in G.values()); print("sec per genome median",s[len(s)//2],"max",s[-1])
