"""[P] Post-hoc (2026-11-25), for the preprint Discussion: C33 among resolved, confirmatory UAA=Stop vertebrate genomes, Elopomorpha against the rest. Counts only; no test, nothing locked changes.
Reads results/mitogenome_master.tsv (locked sweep) and data/ncbi/taxonomy (nodes, merged, names); Elopomorpha = descendants of NCBI taxon 32521. Gains: the T3 counter (src/origins.py, minimum gains on the NCBI hierarchy, no reversals; gains_min and gains_max), run as in src/tiers.py on one state per taxid: 1 if any genome of the taxid is C33, 0 if U33; G33 genomes are left out. Writes results/posthoc2/elopomorpha_c33.tsv"""
import csv,collections,sys
sys.path.insert(0,"src"); import origins
R=[r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included" and r["role"]=="confirmatory" and r["acc"]!="AP002973.1" and r["uaa"]=="Stop" and r["n33"]=="resolved" and r["phylum"]=="Chordata"]
need={r["taxid"] for r in R}; LIN={}
PA={}
for l in open("data/ncbi/taxonomy/nodes.dmp"):
    f=l.split("\t|\t",2); PA[f[0].strip()]=f[1].strip()
MG={}
for l in open("data/ncbi/taxonomy/merged.dmp"):
    f=[x.strip() for x in l.split("|")]; MG[f[0]]=f[1]
for t in need:
    x=MG.get(t,t); a=set()
    while x in PA and x not in a: a.add(x); x=PA[x]
    if a: LIN[t]=a
ELO,VERT,OPH="32521","7742",None
miss=[r["acc"] for r in R if r["taxid"] not in LIN]
V=[r for r in R if VERT in LIN.get(r["taxid"],())]; E=[r for r in V if ELO in LIN[r["taxid"]]]; O=[r for r in V if ELO not in LIN[r["taxid"]]]
def row(name,x):
    c=[r for r in x if r["base"]=="C"]; g=[r for r in x if r["base"]=="G"]
    return dict(group=name,resolved_stop_genomes=len(x),taxa=len({r["taxid"] for r in x}),C33=len(c),G33=len(g),other_nonU=sum(r["base"] not in "UCG" for r in x),C33_pct=round(100*len(c)/len(x),3) if x else "",C33_species="; ".join(sorted(r["organism"] for r in c)))
PAR=origins.load_nodes()
def gains(x):
    st={r["taxid"]:0 for r in x if r["base"]=="U"}; st.update({r["taxid"]:1 for r in x if r["base"]=="C"})
    mixed=sorted({r["taxid"] for r in x if r["base"]=="C"}&{r["taxid"] for r in x if r["base"]=="U"})
    return (origins.count(st,PAR),mixed) if any(st.values()) else (dict(gains_min=0,gains_max=0),mixed)
out=[row("Vertebrata",V),row("Elopomorpha",E),row("other vertebrates",O)]
for od in sorted({r["order"] for r in E}): out.append(row("Elopomorpha: "+(od or "no order"),[r for r in E if r["order"]==od]))
for fa in sorted({r["family"] for r in E if any(q["family"]==r["family"] and q["base"]=="C" for q in E)}): out.append(row("Elopomorpha family: "+fa,[r for r in E if r["family"]==fa]))
SUB={}
for l in open("data/ncbi/taxonomy/nodes.dmp"):
    f=[x.strip() for x in l.split("|")]
    if f[2]=="subfamily": SUB[f[0]]=1
NM={}
for l in open("data/ncbi/taxonomy/names.dmp"):
    f=[x.strip() for x in l.split("|")]
    if f[3]=="scientific name" and f[0] in SUB: NM[f[0]]=f[1]
for nm in ("Ophichthinae",):
    ids=[k for k,v in NM.items() if v==nm]; x=[r for r in E if ids and ids[0] in LIN[r["taxid"]]]; out.append(row("Elopomorpha subfamily: "+nm,x))
for g in ("Megalops","Elops","Anguilla","Ophichthus","Conger"): out.append(row("genus "+g,[r for r in E if r["genus"]==g]))
G={"Vertebrata":V,"Elopomorpha":E,"other vertebrates":O}; G.update({"Elopomorpha: "+od:[r for r in E if r["order"]==od] for od in {r["order"] for r in E}})
for o in out:
    if o["group"] in G: g,mx=gains(G[o["group"]]); o["gains_min"],o["gains_max"],o["taxids_with_both_C_and_U"]=g["gains_min"],g["gains_max"],";".join(mx)
    else: o["gains_min"]=o["gains_max"]=o["taxids_with_both_C_and_U"]=""
# the same without the one unconfirmed C33 call (Xyrias revulsus)
g2,_=gains([r for r in E if r["organism"]!="Xyrias revulsus"]); print("Elopomorpha without Xyrias:",g2)
g3,_=gains([r for r in E if r["organism"] not in("Xyrias revulsus",)]+[r for r in O]); print("vertebrates without Xyrias:",g3)
with open("results/posthoc2/elopomorpha_c33.tsv","w",newline="") as f: w=csv.DictWriter(f,list(out[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(out)
print("no lineage:",len(miss),miss[:5])
for o in out: print(o["group"],o["resolved_stop_genomes"],"taxa",o["taxa"],"C33",o["C33"],"G33",o["G33"],"pct",o["C33_pct"],"gains",o["gains_min"],o["gains_max"],o["taxids_with_both_C_and_U"],"|",o["C33_species"])
