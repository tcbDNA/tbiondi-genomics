"""Build the Gate 3 sweep: one FASTA per representative genome + manifest + annotation features. No caller is run here.
Every representative gets a row in results/gate3/sweep_audit.tsv with included / excluded:reason."""
import csv,gzip,glob,json,os,re,sys,warnings,collections; warnings.filterwarnings("ignore")
from Bio import SeqIO
TAX={}
for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
    f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; TAX[f[0]]=f
S=list(csv.DictReader(open("results/gate1_species.tsv"),delimiter="\t")); tid={r["accession"]:r["taxid"] for fn in("results/refseq_mito_inventory.tsv","results/organelle_mito_inventory.tsv") for r in csv.DictReader(open(fn),delimiter="\t")}
REP={}
for r in S: REP[r["accession"]]=dict(acc=r["accession"],org=r["organism"],species=r["species"],taxid=tid[r["accession"]],prov={"RefSeq":"RefSeq","DToL":"DToL"}.get(r["set"],"other_assembly"))
for r in csv.DictReader(open("results/insdc_only_set.tsv"),delimiter="\t"):
    if r["eligible"]=="yes": REP[r["accession"]]=dict(acc=r["accession"],org=r["species"],species=r["species"],taxid=r["taxid"],prov="INSDC_only")
NAME={f[1]:t for t,f in TAX.items()}      # taxids are looked up by name, never typed in
for a,sp in(("AP002973.1","Albula glossodonta"),("LR745847.1","Asterias rubens")):
    if a not in REP: REP[a]=dict(acc=a,org=sp,species=sp,taxid=NAME[sp],prov="named_extra")
P=list(csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")); PID={r["panel_id"] for r in P}; bi=lambda s:" ".join(re.sub(r"[^A-Za-z. ]"," ",s).split()[:2]).lower()
PSP={bi(r["species"]) for r in P}|{"plicatula australis","anomia chinensis","placuna placenta"}; KNOWN={"anomia","enigmonia","placuna","plicatula","cephalodiscus","radopholus","pseudotanais"}
BAN={"PP708078","PP708091","OZ576062"}
seqs=0; os.makedirs("data/sweep/fa",exist_ok=True); os.makedirs("results/gate3",exist_ok=True); man=[]; aud=[]; F=open("data/sweep/features.jsonl","w")
def feats(r):
    d={"cds":[],"trna":[]}
    for f in r.features:
        if f.type=="CDS": d["cds"].append([int(f.location.start)+1,int(f.location.end),"+" if f.location.strand==1 else "-",(f.qualifiers.get("gene") or f.qualifiers.get("product") or ["?"])[0]])
        elif f.type=="tRNA":
            t=str(f.qualifiers.get("product",""))+str(f.qualifiers.get("gene","")); d["trna"].append([int(f.location.start)+1,int(f.location.end),"+" if f.location.strand==1 else "-",(f.qualifiers.get("product") or f.qualifiers.get("gene") or ["?"])[0],bool(re.search(r"Tyr|trnY",t,re.I))])
    return d if d["cds"] or d["trna"] else None
done=set()
def put(acc,seq,topo,ft):
    if acc not in REP or acc in done: return
    done.add(acc); R=REP[acc]; seq=re.sub(r"[^ACGTN]","N",seq.upper()); L=len(seq); nN=seq.count("N"); ex=""
    if acc.split(".")[0] in BAN: ex="excluded_by_prereg"
    elif not 8000<=L<=200000: ex="length_outside_8_200kb"
    elif nN>0.01*L: ex="ambiguous_gt1pct"
    f=TAX.get(R["taxid"]) or [""]*9; sp=f[2] or f[1] or R["species"]
    panel=acc in PID; psp=bi(sp) in PSP or bi(R["org"]) in PSP; kn=sp.split(" ")[0].lower() in KNOWN
    sh="data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)
    if not ex:
        os.makedirs(sh,exist_ok=True); p=f"{sh}/{acc}.fna"; open(p,"w").write(f">{acc}\n{seq}\n"); man.append([acc,p,topo]); F.write(json.dumps(dict(acc=acc,f=ft))+"\n")
    aud.append([acc,R["org"],sp,R["taxid"],f[7],f[6],f[5],f[4],f[3],R["prov"],topo,L,nN,"panel" if panel else "panel_species" if psp else "known_reassigned_genus" if kn else "confirmatory","included" if not ex else ex])
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id in REP:
            try: s=str(r.seq)
            except Exception: continue
            put(r.id,s,r.annotations.get("topology",""),feats(r))
for fn in sorted(glob.glob("data/sweep/gb/*.gb")):
    for r in SeqIO.parse(fn,"genbank"):
        if r.id in REP:
            try: s=str(r.seq)
            except Exception: continue
            put(r.id,s,r.annotations.get("topology",""),feats(r))
topo={r["accession"]:r["topology"] for fn in("results/organelle_mito_inventory.tsv","results/refseq_mito_inventory.tsv") for r in csv.DictReader(open(fn),delimiter="\t")}
for fn in sorted(glob.glob("data/ncbi/organelle_fasta/*.fna")):
    for r in SeqIO.parse(fn,"fasta"):
        a=r.id.split("|")[1] if "|" in r.id else r.id
        if a in REP: put(a,str(r.seq),topo.get(a,""),None)
for a,R in REP.items():
    if a not in done: f=TAX.get(R["taxid"]) or [""]*9; aud.append([a,R["org"],R["species"],R["taxid"],f[7],f[6],f[5],f[4],f[3],R["prov"],"","","","","excluded:no_sequence_obtained"])
F.close(); w=csv.writer(open("results/gate3/manifest.tsv","w"),delimiter="\t"); w.writerow(["acc","path","topology"]); w.writerows(man)
w=csv.writer(open("results/gate3/sweep_audit.tsv","w"),delimiter="\t"); w.writerow("acc organism species taxid phylum class order family genus provenance topology length nN role status".split()); w.writerows(aud)
print("representatives",len(REP),"included",len(man)); print(collections.Counter(a[-1] for a in aud)); print("roles",collections.Counter(a[-2] for a in aud if a[-1]=="included")); print("provenance",collections.Counter(a[9] for a in aud if a[-1]=="included"))
