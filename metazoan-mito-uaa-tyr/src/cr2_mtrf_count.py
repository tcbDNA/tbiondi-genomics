"""Pre-push fix 3 (POST-HOC, 2026-11-25): the number of UAA=Stop species with a UAA call in this study and a readable mtRF1a (MTRF1L) position 216. Computed from results/gate2/mtrf_loop.tsv, results/gate4/loop.tsv, the Gate 2 fetch audit and the UAA classes (mitogenome_master.tsv; Gate 2 panel counter + Codetta for panel genomes that are not in the sweep; UNVERIFIED master for M. hapla).
writes results/posthoc2/mtrf1a_stop_denominator.tsv"""
import csv,json,collections,sys
sys.path.insert(0,"src"); import gate3_classify as K
M={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")}; U={r["acc"]:r for r in csv.DictReader(open("results/gate3_unverified/master.tsv"),delimiter="\t")}
up={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/uaa_panel.jsonl")}; cp={}
for l in open("results/gate2/codetta_mitodb_lin.jsonl"): x=json.loads(l); cp[x["panel_id"]]=x
def panel_class(pid):
    if pid in M and M[pid]["status"]=="included": return M[pid]["uaa"],"sweep table (panel role)"
    return K.uaa_class(dict(uaa=up.get(pid),codetta=dict(cp.get(pid,{}),status=cp.get(pid,{}).get("status","ok"))))[0],"Gate 2 panel output"
F={r["accession"]:r for r in csv.DictReader(open("results/audit/mtrf_bg_fetch.tsv"),delimiter="\t")}
rows=[]; d=collections.defaultdict(list)
for r in csv.DictReader(open("results/gate2/mtrf_loop.tsv"),delimiter="\t"):
    if r["covers"]=="True" and r["source"].startswith("BG_"): d[r["source"][3:]].append(r)
for a,v in sorted(d.items()):
    L=[r for r in v if r["cls"]=="HUMAN_MTRF1L"]; c,src=panel_class(F[a]["panel_id"]); rows.append(["Gate 2 background",F[a]["species"],F[a]["panel_id"],a,c,src,"mtRF1a (MTRF1L)" if L else "only MTRF1 (paralogue)","/".join(sorted({r["p216"] for r in (L or v)})),"yes" if L and c=="Stop" else "no"])
for sp,pids,srcn in (("Heteranomia squamula",["LIN_Heteranomia_squamula_S_USA","LIN_Heteranomia_squamula_S_Russia","LIN_Heteranomia_squamula_UK"],"FIG_Hsq"),("Pododesmus macrochisma",["LIN_Pododesmus_macrochisma"],"GENOME_PodMac")):
    c={panel_class(p)[0] for p in pids}; p216={r["p216"] for r in csv.DictReader(open("results/gate2/mtrf_loop.tsv"),delimiter="\t") if r["source"]==srcn and r["covers"]=="True"}; rows.append(["Gate 2 matched Stop sister",sp,";".join(pids),srcn,"/".join(sorted(c)),"Gate 2 panel output","mtRF1a (MTRF1L)","/".join(sorted(p216)),"yes" if c=={"Stop"} else "no"])
G4={"PROT_Albula_glossodonta_BYU":("Albula glossodonta","CM165726.1"),"GEN_Lycenchelys":("Lycenchelys verrillii","CM163399.1"),"GEN_Cr_infrequens":("Craniella infrequens","OZ318138.1"),"GEN_Cr_zetlandica":("Craniella zetlandica","OZ407450.1"),"PROT_Megalops":("Megalops cyprinoides",None),"PROT_Anguilla":("Anguilla anguilla",None),"PROT_Lycodopsis":("Lycodopsis pacificus",None),"PROT_Zoarces":("Zoarces viviparus",None),"GEN_C_apion":("Cinachyrella apion",None),"GEN_P_vulnus":("Pratylenchus vulnus","GQ332425.1"),"PROT_M_hapla":("Meloidogyne hapla",None),"GEN_R_robustus":("Rotylenchus robustus",None),"GEN_H_columbus":("Hoplolaimus columbus",None),"GEN_R_reniformis":("Rotylenchulus reniformis","CM003310.1")}
g=collections.defaultdict(list)
for r in csv.DictReader(open("results/gate4/loop.tsv"),delimiter="\t"):
    if r["covered"]=="True" and r["label"]=="MTRF1L": g[r["source"]].append(r)
for s,(sp,acc) in G4.items():
    if not g[s]: rows.append(["Gate 4",sp,acc or "",s,"","","no covered mtRF1a","","no"]); continue
    H=[m for m in M.values() if m["status"]=="included" and m["species"]==sp] or [dict(m,uaa=m["uaa"]) for m in U.values() if m.get("status")=="included" and m["species"]==sp]; c="/".join(sorted({m["uaa"] for m in H})) or "no call"; src="UNVERIFIED sweep (exploratory)" if H and not any(m["acc"] in M for m in H) else "sweep table" if H else ""
    note="yes" if c=="Stop" and sp!="Hoplolaimus columbus" else "no"; rows.append(["Gate 4",sp,";".join(m["acc"] for m in H)[:60],s,c,src,"mtRF1a (MTRF1L)" if sp!="Hoplolaimus columbus" else "bacterial-like contaminant","/".join(sorted({r["p216"] for r in g[s]})),note])
rows.append(["Gate 2 reads","Rhabdopleura annulata","","reads","no call (Stop in the literature)","","mtRF1a (reads)","T","no"])
w=csv.writer(open("results/posthoc2/mtrf1a_stop_denominator.tsv","w"),delimiter="\t"); w.writerow("set species genome source uaa_class uaa_source paralogue p216 counted".split()); w.writerows(rows)
for k in ("Gate 2 background","Gate 2 matched Stop sister","Gate 4"): x=[r for r in rows if r[0]==k]; print(k,"rows",len(x),"| counted",sum(r[8]=="yes" for r in x),"| not counted:",[(r[1],r[4],r[6]) for r in x if r[8]!="yes"])
y=[r for r in rows if r[8]=="yes"]; print("TOTAL counted",len(y),"distinct species",len({r[1] for r in y}),"residues at 216",collections.Counter(r[7] for r in y))
