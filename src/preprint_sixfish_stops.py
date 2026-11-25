"""[P] Post-hoc (2026-11-25), for the preprint: terminal stops of the 13 mitochondrial protein genes in the six Tier 1 C33 fish. Nothing locked changes.
Rules: results/posthoc2/sixfish_polyA_criteria.md (written before any RNA read was looked at).
Inputs: results/gate3/run/out.jsonl (locked sweep outputs), data/sweep FASTA, GenBank annotations (results/work/r2/spec/records.gb), data/sweep/features_fixed.jsonl (reference fish),
 results/posthoc2/sixfish_polya_junctions.tsv (src/preprint_sixfish_polya.py), if present.
Outputs: results/posthoc2/sixfish_terminal_uaa.tsv (one row per species), sixfish_terminal_genes.tsv (one row per gene and annotation), sixfish_terminal_refs.tsv (reference fish)."""
import sys,os,csv,json,re,statistics,subprocess,collections
sys.path.insert(0,"src"); import uaa_caller as U
from Bio import SeqIO
from Bio.Seq import Seq
FISH=[("Albula glossodonta","CM165726.1",["OZ556087.1"]),("Lycenchelys verrillii","CM163399.1",["QB104173.1","PX940805.1","PV339947.1"]),("Microdevario kubotai","NC_037360.1",["NC_037360.1"]),
 ("Serrivomer beanii","NC_013627.1",["NC_013627.1"]),("Serrivomer sector","NC_013436.1",["NC_013436.1"]),("Stemonidium hypomelas","NC_013628.2",["NC_013628.2"])]
GENES="ND1 ND2 COX1 COX2 ATP8 ATP6 COX3 ND3 ND4L ND4 ND5 ND6 CYTB".split()
AL={"COI":"COX1","COII":"COX2","COIII":"COX3","CO1":"COX1","CO2":"COX2","CO3":"COX3","CYT B":"CYTB","COB":"CYTB"}
def norm(x): x=x.upper(); return AL.get(x,x)
def path(a): return f"data/sweep/fa/{str(sum(map(ord,a))%200).zfill(3)}/{a}.fna"
def rc(s): return U.rc(s)
SEQ={a:U.load(path(a))[1] for _,a,_ in FISH}
GB={r.id:r for r in SeqIO.parse("results/work/r2/spec/records.gb","gb")}
STOPS={"TAA":"UAA","TAG":"UAG","AGA":"AGA","AGG":"AGG"}
def cds(acc):
    """gene -> (strand, start0, end0, nt length in frame) of an annotated record; the record's own sequence"""
    r=GB[acc]; out={}
    for f in r.features:
        if f.type!="CDS": continue
        g=norm((f.qualifiers.get("gene") or f.qualifiers.get("product"))[0]); n=len(f.location)-(int(f.qualifiers.get("codon_start",["1"])[0])-1)
        out[g]=("+" if f.location.strand==1 else "-",int(f.location.start),int(f.location.end),n,str(f.extract(r.seq)).upper())
    return out
def strandseq(s,st): return s if st=="+" else rc(s)
def place(tail,x):
    """best ungapped placement of tail in x; returns (end0, identity)"""
    L=len(tail); best=(-1,0)
    seeds=[(i,tail[i:i+12]) for i in range(0,L-11,6)]; cand=set()
    for i,k in seeds:
        p=x.find(k)
        while p>=0: cand.add(p-i); p=x.find(k,p+1)
    for c in cand:
        if c<0 or c+L>len(x): continue
        m=sum(a==b for a,b in zip(tail,x[c:c+L]))
        if m>best[1]: best=(c+L,m)
    return best[0],best[1]/L
def classify(x,e,n):
    """x = strand sequence of the counted genome, e = end of CDS in x (0-based, exclusive), n = CDS length; class of the stop"""
    k=n%3
    if k==0:
        c=x[e-3:e]; return STOPS.get(c,"other:"+c),c
    c=x[e-k:e]; w=x[e-k:e-k+3]
    return ({"T":"T--","TA":"TA-"}.get(c,"other:"+c)),c+"|"+x[e:e+6]
# ---- sweep outputs
SW={}
for l in open("results/gate3/run/out.jsonl"):
    d=json.loads(l)
    if d["acc"] in SEQ: SW[d["acc"]]=d
# ---- reference fish for the C-terminal check
M=[r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")]
F={}
for l in open("data/sweep/features_fixed.jsonl"):
    d=json.loads(l)
    if d.get("f"): F[d["acc"]]=d["f"]
refs=[]
for od in ("Elopiformes","Anguilliformes","Cypriniformes","Perciformes"):
    gen=set(); n=0
    for r in sorted((r for r in M if r["order"]==od and r["role"]=="confirmatory" and r["status"]=="included" and r["provenance"]=="RefSeq" and r["uaa"]=="Stop" and r["n33"]=="resolved" and r["base"]=="U"),key=lambda r:(0 if r["family"]=="Zoarcidae" or r["genus"] in("Danio","Anguilla") else 1,r["acc"])):
        c=F.get(r["acc"],{}).get("cds") or []
        if r["genus"] in gen or sorted({x[3] for x in c})!=sorted(GENES) or len(c)!=13: continue
        gen.add(r["genus"]); refs.append((od,r["acc"],r["organism"],{x[3]:x for x in c})); n+=1
        if n==4: break
def prot(x,a,b):
    s=x[a:b]; n=len(s); k=n%3; s=s[:n-k] if k else (s[:-3] if s[-3:] in STOPS else s); return str(Seq(s).translate(table=2)).replace("*","X")
RP=collections.defaultdict(dict)
for od,acc,org,c in refs:
    s=U.load(path(acc))[1]
    for g,(a,b,st,_) in c.items():
        x=strandseq(s,st); L=len(s); aa,bb=(a-1,b) if st=="+" else (L-b,L-a+1); RP[g][acc]=prot(x,aa,bb)
def mafft(d):
    r=subprocess.run(["env/prefix/bin/mafft","--auto","--quiet","-"],input="".join(f">{k}\n{v}\n" for k,v in d.items()),capture_output=True,text=True,check=True); out={}
    for b in r.stdout.split(">")[1:]: h,s=b.split("\n",1); out[h.split()[0]]=s.replace("\n","")
    return out
def offset(aln,q,refids):
    ends=[max(i for i,c in enumerate(aln[r]) if c!="-") for r in refids]; med=int(statistics.median(ends)); qe=max(i for i,c in enumerate(aln[q]) if c!="-")
    if qe>=med: return sum(aln[q][i]!="-" for i in range(med+1,qe+1))
    return -sum(sum(aln[r][i]!="-" for r in refids)*2>len(refids) for i in range(qe+1,med+1))
# ---- polyA results
PA={}
if os.path.exists("results/posthoc2/sixfish_polya_junctions.tsv"):
    for r in csv.DictReader(open("results/posthoc2/sixfish_polya_junctions.tsv"),delimiter="\t"): PA[(r["acc"],r["gene"])]=r
rows=[]; summ=[]
for sp,acc,ann in FISH:
    s=SEQ[acc]; L=len(s); name,L0,kept,recs,*_=U.analyse(path(acc),True); fr=U.frames(s,2000); H={str(h.name):h.M for h in U.hmms()}
    dom={}
    for d in kept:
        if d["gene"] not in dom or d["ht"]>dom[d["gene"]]["ht"]: dom[d["gene"]]=d
    tr=SW[acc]["trna"]; T=[(a,b,st) for a,b,st,*_ in tr["aragorn"]+tr["trnascan"]]
    final={}
    for g in GENES:
        calls=[]
        for an in ann:
            c=cds(an)
            if g not in c: rows.append(dict(species=sp,acc=acc,gene=g,annotation=an,note="gene not annotated in this record")); continue
            st,a,b,n,cs=c[g]; x=strandseq(s,st)
            if an==acc: e=b if st=="+" else L-a; idn=1.0
            else:
                e,idn=place(cs[-90:],x+x[:300])
                if idn<0.95: rows.append(dict(species=sp,acc=acc,gene=g,annotation=an,lift_identity=round(idn,3),note="lift-over below 95%, not used")); continue
                e=e if e<=L else e-L
            cl,cod=classify(x,e,n); flag="TAA_overlaps_next_gene" if cl in("TA-",) and x[e:e+3]=="ATG" else ""
            if cl=="UAA" and g=="ATP6" and x[e-1:e+2]=="ATG": flag="TAA_overlaps_next_gene"
            # independent check 1: profile frame
            d=dom.get(g); st2,f,codl,xx,_=fr[d["frame"]]; off=0 if st2=="+" else len(xx)-L; p0=f+3*d["tt"]   # index in xx after last aligned codon
            e_xx=e+off; k=n%3; stop_start=e_xx-(k or 3); inframe=(st2==st and (stop_start-f)%3==0); first=None; j=p0
            while j<stop_start+3*60:
                if xx[j:j+3] in STOPS: first=(xx[j:j+3],j); break
                j+=3
            ahead=[j for j in range(p0,stop_start,3) if xx[j:j+3] in STOPS]
            if k==0: prof="agrees" if inframe and first and first[1]==stop_start else ("frame agrees; earlier stop" if inframe and ahead else "differs")
            else:
                ts=[t for t in T if t[2]==st and abs((t[0]-1 if st=="+" else L-t[1])-e)<=2]; nxt=(g=="ATP6" and x[e:e+3]=="ATG")
                prof="agrees" if inframe and not ahead and (ts or nxt) else ("frame agrees; no tRNA or gene start called at the end" if inframe and not ahead else ("frame agrees; earlier stop" if inframe else "differs"))
            # independent check 2: C-terminal alignment
            a0=e-n; p=prot(x if a0>=0 else x[a0:]+x,max(a0,0),e); p=str(Seq(x[a0:e-(k or 3)]).translate(table=2)).replace("*","X"); D=dict(RP[g]); D["Q"]=p; al=mafft(D); o=offset(al,"Q",list(RP[g])); ncon="conserved" if abs(o)<=5 else "shifted"
            rows.append(dict(species=sp,acc=acc,gene=g,annotation=an,lift_identity=round(idn,3),strand=st,end_1based=e if st=="+" else L-e+1,stop_class=cl,stop_sequence=cod,flag=flag,tail_dist_codons=(stop_start-p0)//3,profile_check=prof,cterm_offset_aa=o,cterm=ncon,note=""))
            calls.append((cl,flag,prof,ncon,an,rows[-1]))
        final[g]=calls
    SWU=SW[acc]["uaa"]; cnt=collections.Counter(); det=collections.defaultdict(list); dis=[]
    for g in GENES:
        cs={c[0] for c in final[g]}
        if len(cs)>1: dis.append(g+":"+"/".join(sorted(f"{c[0]}({c[4]})" for c in final[g])))
        # ATP6 TA|ATG counts as incomplete (the last base of UAA belongs to COX3 in the genome)
        best=[c for c in final[g] if not (c[1]=="TAA_overlaps_next_gene" and c[0]=="UAA")] or final[g]; cl="TA-" if any(c[1] for c in final[g]) else collections.Counter(c[0] for c in best).most_common(1)[0][0]
        ok=any(c[3]=="conserved" and c[2]=="agrees" and (c[0]==cl) for c in final[g]) or (cl=="TA-" and any(c[3]=="conserved" and c[2]=="agrees" for c in final[g]))
        k={"UAA":"uaa","T--":"inc","TA-":"inc","UAG":"uag","AGA":"agr","AGG":"agr"}.get(cl,"oth"); cnt[k]+=1; det[k].append(g+("" if k!="inc" else " "+cl)+("" if ok else "*")); cnt[k+"_ok"]+=ok
    pa=[(g,PA[(acc,g)]) for g in GENES if (acc,g) in PA]; hasrna=bool(pa)
    sup=[g for g,r in pa if r["verdict"].startswith("transcript-supported UAA")]; bic=[g for g,r in pa if r["verdict"].startswith("UAA in bicistronic")]; oth=[f"{g}: {r['verdict']}" for g,r in pa if g not in sup+bic]
    summ.append(dict(species=sp,accession=acc,annotation_source=",".join(ann)+("" if ann==[acc] else " (lifted over)"),internal_UAA=SWU["n_taa_internal"],UAA_at_conserved_Tyr=SWU["n_taa_tyrcol"],conserved_Tyr_sites=SWU["tyrcols"],UAU=SWU["tat"],UAC=SWU["tac"],other_codons=SWU["tyrcol_other"],aligned_codons=SWU["aligned_codons"],
      PCG=13,encoded_UAA=cnt["uaa"],encoded_UAA_boundary_confirmed=cnt["uaa_ok"],encoded_UAA_genes=" ".join(det["uaa"]),incomplete=cnt["inc"],incomplete_boundary_confirmed=cnt["inc_ok"],incomplete_genes=", ".join(det["inc"]),UAG=cnt["uag"],UAG_genes=" ".join(det["uag"]),AGR=cnt["agr"],AGR_genes=" ".join(det["agr"]),other=cnt["oth"],other_genes=" ".join(det["oth"]),
      annotation_disagreements="; ".join(dis) or "none",RNA_evidence=("RNA-seq "+pa[0][1]["run"]) if hasrna else "annotation only (no RNA-seq or TSA found 2026-11-25)",transcript_supported_UAA=len(sup) if hasrna else "NA",transcript_supported_genes=" ".join(sup),UAA_in_bicistronic_transcript=" ".join(bic),transcript_other="; ".join(oth)))
def w(p,R):
    ks=[]
    for r in R:
        for k in r:
            if k not in ks: ks.append(k)
    with open(p,"w",newline="") as f: o=csv.DictWriter(f,ks,delimiter="\t",lineterminator="\n"); o.writeheader(); o.writerows(R)
w("results/posthoc2/sixfish_terminal_genes.tsv",rows); w("results/posthoc2/sixfish_terminal_uaa.tsv",summ)
w("results/posthoc2/sixfish_terminal_refs.tsv",[dict(order=a,accession=b,organism=c) for a,b,c,_ in refs])
for r in summ: print(r["species"],"| internal UAA",r["internal_UAA"],"of",r["conserved_Tyr_sites"],"Tyr sites | UAA",r["encoded_UAA"],"(",r["encoded_UAA_genes"],") inc",r["incomplete"],"(",r["incomplete_genes"],") UAG",r["UAG"],r["UAG_genes"],"AGR",r["AGR"],r["AGR_genes"],"other",r["other"],r["other_genes"],"|",r["annotation_disagreements"],"|",r["RNA_evidence"],r["transcript_supported_UAA"])
