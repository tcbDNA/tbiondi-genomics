"""Gate 2: freeze the control panel. Allowlist for every Gate 2 caller. Selection uses documented code and taxonomy only."""
import csv,glob,gzip,os,random,re,sys,collections,hashlib
sys.path.insert(0,"src"); import ncbi
from Bio import SeqIO
import xml.etree.ElementTree as ET, warnings; warnings.filterwarnings("ignore")
SEED=20261123; random.seed(SEED); os.makedirs("data/panel",exist_ok=True); rows=[]
def add(pid,sp,lin,role,code,src,seq,topo,note=""):
    seq=re.sub(r"[^ACGTN]","N",seq.upper()); open(f"data/panel/{pid}.fna","w").write(f">{pid}\n{seq}\n")
    rows.append([pid,sp,lin,role,code,src,len(seq),topo,hashlib.sha256(seq.encode()).hexdigest()[:16],note])
def gbseq(f):
    t=open(f,errors="replace").read(); return re.sub(r"[^A-Za-z]","",t.split("ORIGIN",1)[1].split("//")[0]), ("circular" if "circular" in t.split("\n",1)[0] else "linear")
# Lin Figshare
LIN=[("Anomia chinensis","Anomia","pos"),("Anomia trigonopsis Mauritania","Anomia","pos"),("Anomia trigonopsis NZ","Anomia","pos"),("Anomia_ephippium_S-France","Anomia","pos"),
("Enigmonia aenigmatica S-purple","Enigmonia","pos"),("Enigmonia_aenigmatica_S-yellow","Enigmonia","pos"),("Placuna ephippium","Placunidae","pos"),("Placuna quadrangula","Placunidae","pos"),("Placuna vitream","Placunidae","pos"),
("Plicatuala_australis","Plicatulidae","pos"),("Plicatula muricata","Plicatulidae","pos"),("Pododesmus patelliformis","Pododesmus","neg"),("Pododesmus_glaucus","Pododesmus","neg"),("Pododesmus_macrochisma","Pododesmus","neg"),
("Heteranomia aculeata F-type","Heteranomia","neg"),("Heteranomia aculeata M-type","Heteranomia","neg")]
for n,l,r in LIN:
    s,t=gbseq(f"data/figshare/mito/{n}.gb"); add("LIN_"+re.sub(r"\W+","_",n),n.replace("_"," "),l,r+"_matched" if r=="neg" else "pos","UAA=Tyr" if r=="pos" else "UAA=Stop","Lin2026 Figshare 33250323",s,t,"+1 frameshift at TAG per Lin" if l=="Heteranomia" else "")
for n in ("Heteranomia_squamula S-USA.fa","Heteranomia_squamula_S-Russia.fa","Heteranomia_squamula-UK.fasta"):
    r=next(SeqIO.parse("data/figshare/mito/"+n,"fasta")); add("LIN_"+re.sub(r"\W+","_",n.rsplit(".",1)[0]),"Heteranomia squamula","Heteranomia","neg_matched","UAA=Stop","Lin2026 Figshare 33250323",str(r.seq),"unstated","+1 frameshift at TAG per Lin; UK = OZ181722.1" if "UK" in n else "+1 frameshift at TAG per Lin")
# Zhao MOESM4 re-assemblies
t=open("data/sources/fulltext/zhao2026_MOESM4.txt").read()
for blk in [b for b in t.split(">") if b.strip()]:
    head=blk.split("\n")[0]
    if "Pseudotanais" in head: nm="ERR15194044_Pseudotanais_sp_a"; seq=re.sub(r"\s","",blk)[len("ERR15194044_Pseudotanais_sp._a"):]
    else: nm=re.sub(r"\W+","_",head.strip()); seq=re.sub(r"\s","","".join(blk.split("\n")[1:]))
    l="Pseudotanais" if "Pseudot" in nm else "Anomia" if "Anomia" in nm else "Placunidae"
    add("ZHAO_"+nm.replace("_complete_mitochondrial_genome",""),nm,l,"pos","UAA=Tyr","Zhao2026 MOESM4",seq,"unstated","re-assembly from raw reads")
# GenBank controls
GB={"MF374802.1":("Cephalodiscidae","pos","UAA=Tyr"),"FN313571.1":("Radopholus","pos","UAA=Tyr"),"OZ281197.1":("Pseudotanais","pos","UAA=Tyr"),"OZ281198.1":("Pseudotanais","pos","UAA=Tyr"),"OZ281199.1":("Pseudotanais","pos","UAA=Tyr"),"OZ281200.1":("Pseudotanais","pos","UAA=Tyr"),
"FN908482.1":("Rhabdopleura","neg_matched","UAA=Stop"),"ON600795.1":("Rhabdopleura","neg_matched","UAA=Stop"),"ON600796.1":("Rhabdopleura","neg_matched","UAA=Stop"),
"LC597489.1":("Tanaidacea","neg_matched","UAA=Stop"),"LC923628.1":("Tanaidacea","neg_matched","UAA=Stop"),"LC923627.1":("Tanaidacea","neg_matched","UAA=Stop"),"LC923626.1":("Tanaidacea","neg_matched","UAA=Stop"),"LC923625.1":("Tanaidacea","neg_matched","UAA=Stop"),"LC923624.1":("Tanaidacea","neg_matched","UAA=Stop"),
"NC_023208.1":("Nematoda","neg_matched","UAA=Stop"),"NC_021120.1":("Nematoda","neg_matched","UAA=Stop"),"NC_056772.1":("Nematoda","neg_matched","UAA=Stop"),"NC_026556.1":("Nematoda","neg_matched","UAA=Stop"),"NC_024097.1":("Nematoda","neg_matched","UAA=Stop"),"NC_025291.1":("Nematoda","neg_matched","UAA=Stop"),"NC_061210.1":("Nematoda","neg_matched","UAA=Stop"),"NC_001328.1":("Nematoda","neg_matched","UAA=Stop")}
root=ET.fromstring(ncbi.call("efetch",db="nuccore",id=",".join(GB),rettype="gb",retmode="xml"))
for s in root.iter("GBSeq"):
    a=s.findtext("GBSeq_accession-version"); l,r,c=GB[a]; add(a,s.findtext("GBSeq_organism"),l,r,c,"GenBank",s.findtext("GBSeq_sequence"),s.findtext("GBSeq_topology"))
# Broad canonical panel
S={r["accession"]:r for r in csv.DictReader(open("results/gate1_species.tsv"),delimiter="\t")}
R=[r for r in csv.DictReader(open("results/refseq_mito_inventory.tsv"),delimiter="\t")]
used={r[0] for r in rows}
def ok(r):
    lin=r["lineage"]; tt=r["transl_tables"]
    return (r["accession"] in S and r["completeness"]=="complete" and r["topology"]=="circular" and int(r["n_cds"])>=12 and r["seq_in_gbff"]=="yes" and int(r["n_ambiguous"])==0
        and tt in ("2","4","5","9","13") and not any(x in lin for x in ("Echinodermata","Hemichordata","Anomiidae","Placunidae","Plicatulidae","Tanaidacea")) and r["accession"] not in used)
QUOTA=dict(Arthropoda=70,Chordata=40,Mollusca=45,Nematoda=30,Platyhelminthes=25,Annelida=25,Cnidaria=22,Porifera=14,Nemertea=8,Bryozoa=6,Acanthocephala=6,Xenacoelomorpha=5,Ctenophora=4,Brachiopoda=4,Chaetognatha=3,Nematomorpha=3,Onychophora=2,Kinorhyncha=2,Tardigrada=2,Rotifera=1,Priapulida=1,Placozoa=1,Entoprocta=1,Gnathostomulida=1,Phoronida=1,Gastrotricha=1,Orthonectida=1)
FORCE=["Pecten maximus","Argopecten purpuratus","Magallana gigas","Crassostrea gigas","Mytilus edulis","Drosophila melanogaster","Homo sapiens","Mus musculus","Danio rerio","Ciona intestinalis","Nematostella vectensis","Hydra vulgaris","Lineus longissimus","Penaeus vannamei","Aplysia californica","Daphnia pulex"]
cand=collections.defaultdict(list)
for r in R:
    if ok(r): cand[S[r["accession"]]["phylum"]].append(r)
pick=[]
for ph,q in QUOTA.items():
    f=[r for r in cand[ph] if r["organism"] in FORCE]; by=collections.defaultdict(list)
    for r in cand[ph]:
        if r not in f: by[(S[r["accession"]]["class"],S[r["accession"]]["order"])].append(r)
    for v in by.values(): random.shuffle(v)
    ks=sorted(by); random.shuffle(ks); sel=list(f)
    while len(sel)<q and any(by[k] for k in ks):
        for k in ks:
            if by[k] and len(sel)<q: sel.append(by[k].pop())
    pick+=sel[:max(q,len(f))]
want={r["accession"] for r in pick}; lin={r["accession"]:r for r in pick}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id in want:
            add(r.id,r.annotations["organism"],S[r.id]["phylum"],"neg_broad","UAA=Stop (annotated table %s)"%lin[r.id]["transl_tables"],"RefSeq r238",str(r.seq),"circular",S[r.id]["class"]+"/"+S[r.id]["order"])
w=csv.writer(open("results/gate2_panel.tsv","w"),delimiter="\t"); w.writerow("panel_id species lineage role documented_code source length topology sha16 note".split()); w.writerows(rows)
print(len(rows),collections.Counter(r[3] for r in rows)); print("broad by phylum",collections.Counter(r[2] for r in rows if r[3]=="neg_broad"))
