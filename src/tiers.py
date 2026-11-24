"""Post-hoc (2026-11-24), cold review items 4 and 5: ONE mechanical evidence rule for the 17 non-U33 Stop hits. Evidence fields are typed below with their source file; the tier is computed, not typed. Writes results/cold_review/tiers.tsv and origins per tier."""
import csv,sys,collections
sys.path.insert(0,"src"); import origins
# for each hit: independent=(specimen/record evidence of the same species from a different individual with the same base; text), ind_reads=(n same base, n U33) from an independent individual, own=(n same base, n U33) own specimen or record run, lineage=(taxon, n other species with a confident call, all same base?, best tier of others)
E={
"NC_037360.1":dict(indep="4 other isolates PV357791-PV357794 (one study) and AB937094; run SRR39401262 572/1",ind_rec=True,ind_reads=(572,1),own=(10,0)),
"NC_015546.1":dict(indep="none (RefSeq copy of AP011402; other kubotai only)",ind_rec=False,ind_reads=None,own=None),
"CM165726.1":dict(indep="HI-2016 fish, RNA-seq SRR13036352 1797/1; own fish PacBio SRR38145622 32/0",ind_rec=False,ind_reads=(1797,1),own=(32,0)),
"NC_005796.1":dict(indep="none; no SRA data",ind_rec=False,ind_reads=None,own=None),
"NC_013627.1":dict(indep="PX939991, PX754814, PZ284732; run SRR36227421 272/0 (other fish, own sequence)",ind_rec=True,ind_reads=(272,0),own=None),
"NC_013436.1":dict(indep="PZ743035, PZ285086, PZ285087, OR482581; run SRR22396609 363/0 (other fish)",ind_rec=True,ind_reads=(363,0),own=None),
"PZ495457.1":dict(indep="none",ind_rec=False,ind_reads=None,own=(419,0)),
"NC_013628.2":dict(indep="PZ285083, PZ285084, PZ285088, PZ285090; run SRR38094758 85/0 (1 A)",ind_rec=True,ind_reads=(85,0),own=None),
"CM163399.1":dict(indep="QB104173, PX940805, PV339947; run SRR37037826 622/1",ind_rec=True,ind_reads=(622,1),own=None),
"PV357903.1":dict(indep="none",ind_rec=False,ind_reads=None,own=(88,0)),
"PX804232.1":dict(indep="none",ind_rec=False,ind_reads=None,own=(200,0)),
"KY021065.1":dict(indep="none; no SRA data",ind_rec=False,ind_reads=None,own=None),
"NC_082198.1":dict(indep="none (RefSeq copy of OR241493, same isolate); no SRA data",ind_rec=False,ind_reads=None,own=None),
"KU726971.1":dict(indep="two populations (reads 176/0, 271/0) but second U33 locus on chromosome II",ind_rec=False,ind_reads=(271,0),own=None,drop="another tRNA-Tyr locus in the genome with a different base (chr II, U33), placed identically by Infernal and ARAGORN"),
"OZ318138.1":dict(indep="independent RNA run SRR33975184 (label confirmed by 1,086,464 diagnostic 31-mers): 49-55 G33/0 U33 aligned; own DNA 40/0",ind_rec=False,ind_reads=(49,0),own=(40,0)),
"OZ407450.1":dict(indep="none",ind_rec=False,ind_reads=None,own=(55,0)),
"NC_047225.1":dict(indep="none (RefSeq copy of MK855100); RAD runs hold no mitochondrial reads",ind_rec=False,ind_reads=None,own=None),
}
# lineage-backed facts are COMPUTED (cold review follow-up, 2026-11-24; the earlier version typed them by hand and had Toxocara, Bufotes and Ophichthidae counts wrong): the smallest NCBI taxon that
# contains another species with a confident (resolved) tRNA-Tyr call in the sweep table; the named extra AP002973 is excluded, panel genomes are included as lineage evidence only. Tier 3 needs every other species in that taxon to have the hit's base, and one of them in Tier 1 or Tier 2.
def lineage():
    par={};rank={}
    for l in open("data/ncbi/taxonomy/nodes.dmp"):
        f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; par[f[0]]=f[1]; rank[f[0]]=f[2]
    nm={}
    for l in open("data/ncbi/taxonomy/names.dmp"):
        f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]
        if f[3].rstrip("\t|").strip()=="scientific name": nm[f[0]]=f[1]
    def sp(t):
        x=t
        while x!="1" and rank.get(x)!="species": x=par.get(x,"1")
        return x if x!="1" else t
    def anc(t):
        x=par.get(t,"1"); out=[]
        while x!="1": out.append(x); x=par.get(x,"1")
        return out
    by=collections.defaultdict(lambda: collections.defaultdict(set)); sid={}
    for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t"):
        if r["status"]!="included" or r["n33"]!="resolved" or r["provenance"]=="named_extra": continue
        s_=sp(r["taxid"]); sid[r["acc"]]=s_
        for x in anc(s_): by[x][s_].add((r["base"],r["role"]))
    out={}
    for a in E:
        s_=sid[a]
        for x in anc(s_):
            oth={k:v for k,v in by[x].items() if k!=s_}
            if oth:
                bases=collections.Counter("/".join(sorted({p[0] for p in v})) for v in oth.values()); panel_only=sum(all(p[1].startswith("panel") for p in v) for v in oth.values())
                out[a]=dict(taxon=f"{rank[x]} {nm[x]}",n=len(oth),bases=dict(bases),panel_only=panel_only,species={nm[k]:sorted({p[0] for p in v}) for k,v in oth.items()},counted=len(oth)-panel_only,uniform=set(bases)=={next(iter(bases))} and len(bases)==1); break
    return out
LIN=lineage()
BASE={a:[r["base"] for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["acc"]==a][0] for a in E}
MIN=10
def t12(e):
    if e.get("drop"): return "dropped"
    ok=lambda r: r and r[0]>=MIN and r[1]<=0.01*(r[0]+r[1])
    if e["ind_rec"] or ok(e["ind_reads"]): return "Tier 1"
    if ok(e["own"]): return "Tier 2"
    return ""
tier={a:t12(e) for a,e in E.items()}
def tier3(a,nmin=1,panel=True):
    L=LIN[a]; n=L["n"] if panel else L["counted"]
    if not (n>=nmin and L["uniform"] and next(iter(L["bases"]))==BASE[a]): return False
    oth={R["acc"]:R for R in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if R["acc"] in E and R["acc"]!=a and R["species"] in L["species"]}
    return any(tier[x] in("Tier 1","Tier 2") for x in oth)
for a in E:
    if tier[a]=="" and tier3(a): tier[a]="Tier 3"
for a in E:
    if tier[a]=="": tier[a]="unconfirmed"
R={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")}
N={"Tier 1":"Tier 1 same-species independent","Tier 2":"Tier 2 same-specimen reads only","Tier 3":"Tier 3 lineage-backed only","unconfirmed":"unconfirmed","dropped":"dropped"}
w=csv.writer(open("results/cold_review/tiers.tsv","w"),delimiter="\t"); w.writerow("acc species family base ncbi_mt_code tier evidence lineage_taxon lineage_other_species lineage_bases lineage_panel_only".split())
code={}
for l in open("data/ncbi/taxonomy/nodes.dmp"):
    f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; code[f[0]]=f[8]
for a in sorted(E,key=lambda a:(R[a]["base"],list(N).index(tier[a]),R[a]["species"])):
    w.writerow([a,R[a]["species"],R[a]["family"],R[a]["base"]+"33",code.get(R[a]["taxid"]),N[tier[a]],E[a]["indep"],LIN[a]["taxon"] if a in LIN else "",LIN[a]["n"] if a in LIN else "",";".join(f"{k}:{v}" for k,v in LIN[a]["bases"].items()) if a in LIN else "",LIN[a]["panel_only"] if a in LIN else ""])
print(collections.Counter((R[a]["base"],tier[a]) for a in E))
par=origins.load_nodes(); C=[r for r in R.values() if r["status"]=="included" and r["role"]=="confirmatory" and r["provenance"]!="named_extra" and r["uaa"]=="Stop" and r["n33"]=="resolved"]; U={r["taxid"] for r in C if r["base"]=="U"}
def og(sel):
    st={t:0 for t in U}; st.update({R[a]["taxid"]:1 for a in E if sel(a)}); return origins.count(st,par)["gains_min"]
for nm,f in (("C33 Tier 1",lambda a:R[a]["base"]=="C" and tier[a]=="Tier 1"),("C33 Tier 1+2",lambda a:R[a]["base"]=="C" and tier[a] in("Tier 1","Tier 2")),("C33 Tier 1+2+3",lambda a:R[a]["base"]=="C" and tier[a] in("Tier 1","Tier 2","Tier 3")),("G33 Tier 1",lambda a:R[a]["base"]=="G" and tier[a]=="Tier 1"),("G33 Tier 1+2",lambda a:R[a]["base"]=="G" and tier[a] in("Tier 1","Tier 2"))): print(nm,"genomes",sum(f(a) for a in E),"origins",og(f),"families",sorted({R[a]["family"] for a in E if f(a)}))

print("Sensitivity (Tier 3 cutoff >=2 other species, panel counted): ",{BASE[a]+" "+R[a]["species"]:("Tier 3" if tier3(a,2) else "unconfirmed") for a in E if tier[a]=="Tier 3" or tier[a]=="unconfirmed"})
print("Sensitivity (>=2, panel genomes not counted): ",{BASE[a]+" "+R[a]["species"]:("Tier 3" if tier3(a,2,False) else "unconfirmed") for a in E if tier[a]=="Tier 3"})
for a in ("NC_015546.1","NC_005796.1","KY021065.1","NC_082198.1","NC_047225.1"): print(R[a]["species"],BASE[a]+"33",LIN[a]["taxon"],"others",LIN[a]["n"],LIN[a]["bases"],"panel-only",LIN[a]["panel_only"],"C/G others:",{k:v for k,v in LIN[a]["species"].items() if v!=["U"]})
