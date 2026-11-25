"""Pre-push fix 4 (POST-HOC, 2026-11-25): are the five Tier 1 fish whose independent support is another GenBank record really pairs of distinct specimens of one species? GenBank metadata (results/work/r2/spec/records.gb), whole-mitogenome and COX1 divergence, N33, and independent COI/CYTB records by name (barcodes.gb). Also the panel A. vulpes species check.
writes results/posthoc2/conspecific_specimen_audit.tsv and conspecific_specimen_notes.txt"""
import sys,os,subprocess,re,csv,collections,warnings,hashlib; warnings.filterwarnings("ignore")
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import trny_extract as TE
from Bio import SeqIO
B="env/prefix/bin/"; W="results/work/r2/spec/"; o=open("results/posthoc2/conspecific_specimen_notes.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True)); P("POST-HOC (pre-push fix 4)")
R={r.id:r for r in SeqIO.parse(W+"records.gb","genbank")}
from Bio.Seq import Seq
for i,r in R.items():
    if not r.seq.defined:      # contig-type GenBank records carry no sequence; the sweep FASTA has it
        f=[x for x in ("data/sweep/fa/"+str(sum(map(ord,i))%200).zfill(3)+f"/{i}.fna",f"results/work/cons/{i}.fna") if os.path.exists(x)]; r.seq=Seq("".join(l.strip() for l in open(f[0]) if l[0]!=">"))
Bc=[r for r in SeqIO.parse(W+"barcodes.gb","genbank") if r.seq.defined]
def meta(r):
    q={}
    for f in r.features:
        if f.type=="source": q={k:v[0] for k,v in f.qualifiers.items()}; break
    refs=r.annotations.get("references",[]); sub=[x for x in refs if x.title.lower().startswith("direct")]; j=sub[-1].journal if sub else ""; m=re.match(r"Submitted \(([^)]*)\)\s*(.*)",j,re.S); inst=(m[2] if m else "").replace("\n"," "); inst=re.sub(r"Contact:[A-Z][a-z]+ [A-Z][a-z]+ ","",inst)
    return dict(organism=q.get("organism",""),voucher=(q.get("specimen_voucher") or "").replace(" ",""),isolate=q.get("isolate",""),geo=(q.get("geo_loc_name") or q.get("country") or ""),lat_lon=q.get("lat_lon",""),date=q.get("collection_date",""),biosample=";".join(x.split(":")[1] for x in r.dbxrefs if x.startswith("BioSample")),bioproject=";".join(x.split(":")[1] for x in r.dbxrefs if x.startswith("BioProject")),submitted=(m[1] if m else ""),institution=inst[:70],tech=";".join(str(dict(v).get("Sequencing Technology","")) for v in r.annotations.get("structured_comment",{}).values()))
def blast(q,s,task="megablast",ev="1e-20"):
    open(W+"q.fa","w").write(f">q\n{q}\n"); open(W+"s.fa","w").write(f">s\n{s}\n"); return [l.split("\t") for l in subprocess.run([B+"blastn","-task",task,"-query",W+"q.fa","-subject",W+"s.fa","-evalue",ev,"-dust","no","-outfmt","6 pident length mismatch gaps qstart qend sstart send bitscore"],capture_output=True,text=True).stdout.splitlines()]
def mito(a,b):
    """identity over non-overlapping query intervals (best HSPs first); returns aligned bp, mismatches, gap bases"""
    used=[];L=mm=gp=0
    for h in sorted(blast(a,b),key=lambda h:-float(h[8])):
        u,v=int(h[4]),int(h[5])
        if any(min(v,y)-max(u,x)>20 for x,y in used): continue
        used.append((u,v)); L+=int(h[1]); mm+=int(h[2]); gp+=int(h[3])
    return L,mm,gp
def cox1(r):
    for f in r.features:
        if f.type=="CDS" and re.search(r"cox1|coi\b|co1\b|cytochrome c? ?oxidase subunit (1|I)$",((f.qualifiers.get("gene") or [""])[0]+" "+(f.qualifiers.get("product") or [""])[0]).lower().strip(),re.I): return str(f.extract(r.seq)).upper()
COX={i:cox1(r) for i,r in R.items()}
def cox_pair(a,b):
    qa=COX.get(a) or ""; sb=str(R[b].seq).upper()
    if not qa: qa,sb=COX.get(b) or "",str(R[a].seq).upper()
    h=sorted(blast(qa,sb,"blastn","1e-10"),key=lambda h:-float(h[8]))
    return (int(h[0][1]),int(h[0][2])+int(h[0][3])) if h else (0,0)
def n33(i):
    p=f"results/work/cons/{i}.fna"; q=W+"t.fna"; open(q,"w").write(f">{i}\n{str(R[i].seq).upper()}\n"); d=TE.call(q); return (d.get("N33"),d.get("score")) if d.get("status")=="called" else ("not called","")
SETS=[("Serrivomer sector","NC_013436.1",["PZ743035.1","PZ285086.1","PZ285087.1","OR482581.1"]),("Serrivomer beanii","NC_013627.1",["PX939991.1","PX754814.1","PZ284732.1"]),("Lycenchelys verrillii","CM163399.1",["QB104173.1","PX940805.1","PV339947.1"]),("Stemonidium hypomelas","NC_013628.2",["PZ285083.1","PZ285084.1","PZ285088.1","PZ285090.1"]),("Microdevario kubotai","NC_037360.1",["PV357791.1","PV357792.1","PV357793.1","PV357794.1","AB937094.1"])]
CONG={"Serrivomer sector":["NC_013627.1","PZ495457.1","NC_013628.2"],"Serrivomer beanii":["NC_013436.1","PZ495457.1","NC_013628.2"],"Stemonidium hypomelas":["NC_013436.1","NC_013627.1"],"Microdevario kubotai":["NC_015546.1"],"Lycenchelys verrillii":[]}
w=csv.writer(open("results/posthoc2/conspecific_specimen_audit.tsv","w"),delimiter="\t"); w.writerow("species role record length filed_as voucher isolate biosample bioproject geo lat_lon collection_date submitted institution sequencing N33 trnY_score mito_aligned_bp mito_mismatches mito_gap_bases mito_identity_pct cox1_bp cox1_differences cox1_identity_pct distinct_specimen same_species note".split())
def row(sp,role,i,key):
    m=meta(R[i]); n=n33(i); s=str(R[i].seq).upper()
    if key is None: L=mm=gp=cb=cd=""; pid=cp=""
    else:
        L,mm,gp=mito(str(R[key].seq).upper(),s); pid=round(100*(1-mm/L),2) if L else ""; cb,cd=cox_pair(key,i); cp=round(100*(1-cd/cb),2) if cb else ""
    return m,n,[sp,role,i,len(s),m["organism"],m["voucher"],m["isolate"],m["biosample"],m["bioproject"],m["geo"],m["lat_lon"],m["date"],m["submitted"],m["institution"],m["tech"],n[0],n[1],L,mm,gp,pid,cb,cd,cp]
for sp,key,ind in SETS:
    km,kn,kr=row(sp,"counted record (Tier 1)",key,None); w.writerow(kr+["","",""]); P(f"\n== {sp}: counted record {key}; voucher {km['voucher'] or '-'}, isolate {km['isolate'] or '-'}, {km['institution']}, {km['submitted']}; N33 {kn[0]}")
    for i in ind:
        m,n,r=row(sp,"other specimen",i,key); ids=[x for x in (m["voucher"],m["isolate"],m["biosample"]) if x]; kid=[x for x in (km["voucher"],km["isolate"],km["biosample"]) if x]
        identical=str(R[i].seq)==str(R[key].seq); diff=r[18]+r[19]; dist="yes" if not identical and not set(ids)&set(kid) and diff>0 else "CHECK"; same="yes" if r[23]!="" and r[23]>=98 and r[20]>=97.5 else "not confirmed (COX1 <98%)"
        note=[]
        if m["organism"]!=sp: note.append("filed as "+m["organism"])
        if not ids: note.append("no voucher, isolate or BioSample in the record")
        if not kid: note.append("the counted record has no voucher, isolate or BioSample; distinctness rests on submitter, date and sequence")
        w.writerow(r+[dist,same,"; ".join(note)]); P("  ",i,"| filed as",m["organism"],"| voucher",m["voucher"] or "-","| isolate",m["isolate"] or "-","|",m["biosample"] or "-","|",m["geo"],m["lat_lon"],m["date"],"|",m["institution"][:40],"| mito",f"{r[20]}% ({r[18]} mismatches, {r[19]} gap bases over {r[17]})","| COX1",f"{r[23]}% ({r[22]} of {r[21]})","| N33",n[0],"| distinct",dist,"same species",same,note)
    # among the other specimens: are any two the same sequence?
    H=collections.defaultdict(list); [H[hashlib.md5(str(R[i].seq).upper().encode()).hexdigest()].append(i) for i in [key]+ind]; P("   identical sequences within the set:",[v for v in H.values() if len(v)>1] or "none")
    for c in CONG[sp]:
        L,mm,gp=mito(str(R[key].seq).upper(),str(R[c].seq).upper()); cb,cd=cox_pair(key,c); P("   nearest relatives for scale:",c,R[c].annotations["organism"],f"mito {100*(1-mm/L):.1f}% over {L}" if L else "-",f"COX1 {100*(1-cd/cb):.1f}%" if cb else "-"); w.writerow([sp,"relative, for scale",c,len(R[c].seq),R[c].annotations["organism"]]+[""]*12+[L,mm,gp,round(100*(1-mm/L),2) if L else "",cb,cd,round(100*(1-cd/cb),2) if cb else "","","",""])
    # independent barcodes by name
    qa=COX.get(key) or (next((COX[i] for i in ind if COX.get(i)),"")); skip={x.split(".")[0] for x in [key]+ind}|{"AP007250","AP010857","AP010858","AP011403"}
    open(W+"bc.fa","w").write("".join(f">{r.id}\n{str(r.seq).upper()}\n" for r in Bc if r.id.split(".")[0] not in skip)); open(W+"q.fa","w").write(f">q\n{qa}\n"); org={r.id:r for r in Bc}; best=collections.defaultdict(list)
    for l in subprocess.run([B+"blastn","-task","blastn","-query",W+"q.fa","-subject",W+"bc.fa","-evalue","1e-20","-dust","no","-max_target_seqs","2000","-outfmt","6 sseqid pident length"],capture_output=True,text=True).stdout.splitlines():
        a,p,n=l.split("\t")
        if int(n)>=400: best[org[a].annotations["organism"]].append((float(p),int(n),a))
    B5=sorted(((max(v)[0],k,len(v),max(v)[2],round(min(x[0] for x in v),1)) for k,v in best.items()),reverse=True)[:5]; P("   independent COI records (>=400 bp), best identity by name:",[(k,f"best {p:.1f}% {a}",f"n={n}",f"worst {mn}%") for p,k,n,a,mn in B5])
P("\n== Panel A. vulpes NC_082992.1 (= OR546140.1)"); m=meta(R["OR546140.1"]); P("  ",m)
qs=str(R["NC_082992.1"].seq).upper(); P("   copy of OR546140.1 identical:",qs==str(R["OR546140.1"].seq).upper(),"| N33",n33("NC_082992.1"))
for c in ("CM165726.1","OZ556087.1","AP002973.1"): L,mm,gp=mito(qs,str(R[c].seq).upper()); P("   mitogenome against",c,R[c].annotations["organism"],f"{100*(1-mm/L):.1f}% over {L}")
open(W+"bc.fa","w").write("".join(f">{r.id}\n{str(r.seq).upper()}\n" for r in Bc if r.annotations["organism"].startswith("Albula"))); open(W+"q.fa","w").write(f">q\n{qs}\n"); org={r.id:r for r in Bc}; best=collections.defaultdict(list)
for l in subprocess.run([B+"blastn","-task","blastn","-query",W+"q.fa","-subject",W+"bc.fa","-evalue","1e-20","-dust","no","-max_target_seqs","3000","-outfmt","6 sseqid pident length qstart qend"],capture_output=True,text=True).stdout.splitlines():
    a,p,n,u,v=l.split("\t")
    if int(n)>=300: best[(org[a].annotations["organism"],re.sub(r".*(cytochrome b|cytb|COI|cytochrome (c )?oxidase subunit (1|I)|12S|16S|control region|D-loop|NADH).*",lambda z:z[1],org[a].description,flags=re.I)[:24].lower())].append((float(p),int(n),a))
T=sorted(((max(v)[0],k,len(v),max(v)[1],max(v)[2],round(min(x[0] for x in v),1)) for k,v in best.items()),reverse=True)[:16]
for p,k,n,ln,a,mn in T: P("   ",k,f"best {p:.1f}% over {ln} nt ({a}), n={n}, worst {mn}%")
geo=collections.Counter()
for (sp_,g_),v in best.items():
    if sp_=="Albula vulpes":
        for p,n,a in v:
            q=[f for f in org[a].features if f.type=="source"][0].qualifiers; geo[(g_,"same" if p>=97.5 else "other",(q.get("geo_loc_name") or q.get("country") or ["?"])[0].split(":")[0])]+=1
P("   records filed as A. vulpes by place and whether >=97.5% to the panel genome:",geo.most_common(14))
P("\nHow the links were made: the record lists came from scripts that searched GenBank by genus, called tRNA-Tyr on every record and printed the target species (src/gate3b_conspecific.py -> results/gate3/conspecific_b.txt; src/conspecific_filter.py for Microdevario). Whether a record is another specimen of the same species was then read by hand from titles and vouchers, and typed into src/tiers.py (ind_rec=True and the 'indep' text) and src/dossiers.py. It was not computed until this audit.")
# ---- extras: which M. kubotai lineage the read specimen (SRR39401262) belongs to; distances among the M. kubotai records; A. vulpes against A. goreensis
P("\n== Extras")
K=["NC_037360.1","PV357791.1","PV357792.1","PV357793.1","PV357794.1","AB937094.1"]
for i,a in enumerate(K):
    P("   M. kubotai mitogenome identity (%):",a,[(b,round(100*(1-mito(str(R[a].seq).upper(),str(R[b].seq).upper())[1]/mito(str(R[a].seq).upper(),str(R[b].seq).upper())[0]),1)) for b in K[i+1:]])
f="results/work/r2/hc/M_kubotai/all.fa"
if os.path.exists(f):
    s="";k=False
    for l in open(f):
        if l[0]==">": k=l.startswith(">NODE_206418_")
        elif k: s+=l.strip()
    P("   assembled mitochondrial contig of SRR39401262 (",len(s),"bp) against each record:",[(b,round(100*(1-mito(s,str(R[b].seq).upper())[1]/mito(s,str(R[b].seq).upper())[0]),2),mito(s,str(R[b].seq).upper())[1]) for b in K])
open(W+"q.fa","w").write(f">q\n{str(R['NC_082992.1'].seq).upper()}\n"); g=collections.defaultdict(list)
for l in subprocess.run([B+"blastn","-task","blastn","-query",W+"q.fa","-subject",W+"bc.fa","-evalue","1e-20","-dust","no","-max_target_seqs","3000","-outfmt","6 sseqid pident length"],capture_output=True,text=True).stdout.splitlines():
    a,p,n=l.split("\t"); ds=org[a].description.lower()
    if int(n)>=400: g[(org[a].annotations["organism"],"COI" if re.search(r"coi|co1|cox1|oxidase subunit (1|i)\b",ds) else "CYTB" if re.search(r"cytb|cytochrome b",ds) else "other")].append((float(p),a))
for m_ in ("COI","CYTB"): P("  ",m_,"of the panel A. vulpes, best identity by name:",sorted(((max(v)[0],k[0],len(v)) for k,v in g.items() if k[1]==m_),reverse=True)[:8])

P("   S. sector: all five records sit in one COI group (97.5-99.8% to 11 independent S. sector barcodes); that group also holds one record filed S. jesperseni (GU440522), one S. sp. (JN641138) and one misfiled Stemonidium (PZ575687, India); the next named group is S. beanii at about 92%.")
P("   Read specimens against the counted record (mitogenome): S. sector SRR22396609 = OR482581 98.7%; S. beanii SRR36227421 = PX754814 99.2%; St. hypomelas SRR38094758 = PZ285090 99.5%; L. verrillii SRR37037826 = PX940805 99.6%; M. kubotai SRR39401262 97.6% (99.8% to PV357792, kub2).")
