"""Gate 1: dedupe mitogenomes and intersect with nuclear assemblies / OrthoDB / RNA-seq. Metadata only."""
import csv, gzip, collections, sys, re
csv.field_size_limit(10**9)
TAX={}
for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
    f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; TAX[f[0]]=f  # 0 taxid 1 name 2 species 3 genus 4 family 5 order 6 class 7 phylum 8 kingdom
name2tid={}
sp_name=lambda f: f[2] or f[1]
sp2tid=collections.defaultdict(set)
for t,f in TAX.items(): name2tid[f[1]]=t
def info(t):
    f=TAX.get(str(t)); 
    if not f: return None
    s=sp_name(f); return dict(species=s,genus=f[3] or s.split()[0],phylum=f[7],cls=f[6],order=f[5])
def ph(i): return i["phylum"] or "(no phylum rank)"
R=list(csv.DictReader(open("results/refseq_mito_inventory.tsv"),delimiter="\t"))
O=list(csv.DictReader(open("results/organelle_mito_inventory.tsv"),delimiter="\t"))
aud=csv.writer(open("results/audit/gate1_dedup.tsv","w"),delimiter="\t"); aud.writerow(["accession","organism","taxid","set","decision","reason"])
derived={r["derived_from"] for r in R if r["derived_from"]}
recs=[]
for r in R: recs.append(dict(acc=r["accession"],org=r["organism"],tax=r["taxid"],set="RefSeq",topo=r["topology"],len=int(r["length"]),comp=r["completeness"],ncds=int(r["n_cds"]),tt=r["transl_tables"]))
drop=0
for r in O:
    if r["accession"].split(".")[0] in derived: aud.writerow([r["accession"],r["organism"],r["taxid"],"organelle","drop","duplicate_of_RefSeq_accession"]); drop+=1; continue
    recs.append(dict(acc=r["accession"],org=r["organism"],tax=r["taxid"],set=r["dtol"],topo=r["topology"],len=int(r["length"]),comp="circular" if r["topology"]=="circular" else "linear",ncds=int(r["n_cds"]),tt=""))
for r in recs:
    i=info(r["tax"]); r["i"]=i
    if not i: aud.writerow([r["acc"],r["org"],r["tax"],r["set"],"keep_no_taxonomy","taxid_not_in_dump"])
groups=collections.defaultdict(list)
OPEN=re.compile(r"\b(sp|cf|aff|nr|gen)\.|\bsp\b|uncultured|environmental|hybrid| x ")
for r in recs:
    k=(r["i"] or {}).get("species",r["org"])
    if OPEN.search(r["org"]) or OPEN.search(k): k=k+"|"+r["acc"]
    r["open"]="|" in k; groups[k].append(r)
keep=[]
for sp,g in groups.items():
    g.sort(key=lambda r:(r["set"]!="RefSeq",r["topo"]!="circular",-r["len"])); keep.append(g[0])
    for r in g[1:]: aud.writerow([r["acc"],r["org"],r["tax"],r["set"],"redundant","species_already_represented"])
    if not (g[0]["i"]): pass
print("records",len(R),len(O),"dup accession dropped",drop,"species-level representatives",len(keep))
miss=[r for r in recs if not r["i"]]; print("no taxonomy:",len(miss))
# nuclear resources by species and genus
A=collections.defaultdict(set); Ag=collections.defaultdict(set); Aann=set(); Achr=set()
for a in csv.DictReader(open("results/assemblies_metazoa.tsv"),delimiter="\t"):
    i=info(a["taxid"]); 
    if not i: continue
    A[i["species"]].add(a["accession"]); Ag[i["genus"]].add(a["accession"])
    if a["annotated"]=="yes": Aann.add(i["species"])
    if a["level"] in ("Chromosome","Complete Genome"): Achr.add(i["species"])
Odb=set(); Odbg=set()
for l in gzip.open("data/orthodb/odb13v0_species.tab.gz","rt"):
    i=info(l.split("\t")[0]); 
    if i: Odb.add(i["species"]); Odbg.add(i["genus"])
RNA=collections.Counter()
for l in open("data/ena/invert_rnaseq.tsv"):
    f=l.rstrip("\n").split("\t")
    if f[0].startswith("run"): continue
    i=info(f[1]); 
    if i: RNA[i["species"]]+=1
RNAg=collections.Counter()
for s,c in RNA.items(): pass
for l in ():
    pass
g_rna=collections.defaultdict(int)
for s,c in RNA.items(): g_rna[s.split()[0]]+=c
w=csv.writer(open("results/gate1_species.tsv","w"),delimiter="\t")
w.writerow("accession organism species genus phylum class order set topology length completeness n_cds asm_species asm_genus annotated chromosome_level orthodb rnaseq_runs".split())
T=collections.OrderedDict()
def bump(p,k): T.setdefault(p,collections.Counter())[k]+=1
for r in keep:
    i=r["i"] or dict(species=r["org"],genus=r["org"].split()[0],phylum="",cls="",order="")
    p=ph(i) if r["i"] else "(unmapped)"
    sp,gn=i["species"],i["genus"]
    f=dict(asm_sp=len(A.get(sp,())),asm_gn=len(Ag.get(gn,())),ann=sp in Aann,chr=sp in Achr,odb=sp in Odb,rna=RNA.get(sp,0))
    w.writerow([r["acc"],r["org"],sp,gn,p,i["cls"],i["order"],r["set"],r["topo"],r["len"],r["comp"],r["ncds"],f["asm_sp"],f["asm_gn"],f["ann"],f["chr"],f["odb"],f["rna"]])
    bump(p,"mito"); bump(p,"circular") if r["topo"]=="circular" else None
    for k,c in (("asm_species",f["asm_sp"]>0),("asm_genus",f["asm_gn"]>0),("annotated",f["ann"]),("chrom",f["chr"]),("orthodb",f["odb"]),("rna",f["rna"]>0)):
        if c: bump(p,k)
    if f["asm_sp"] or f["odb"] or f["rna"]: bump(p,"any_nuc_source")
    elif f["asm_gn"]: bump(p,"genus_only")
H=["mito","circular","asm_species","asm_genus","annotated","chrom","orthodb","rna","any_nuc_source","genus_only"]
o=csv.writer(open("results/gate1_phylum_summary.tsv","w"),delimiter="\t"); o.writerow(["phylum"]+H)
for p,c in sorted(T.items(),key=lambda x:-x[1]["mito"]): o.writerow([p]+[c[k] for k in H]); print(p,[c[k] for k in H])
