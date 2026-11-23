"""HMMER filters are OFF (F1=F2=F3=1, no bias filter): with default filters a 723-bit COX1 hit on a translated frame of an AT-rich genome was discarded.
Frame-aware UAA (TAA) meaning caller from profile-HMM alignments of the 13 OXPHOS PCGs (pyhmmer).
Records, per genome: internal TAA codons at aligned match columns, at conserved-Tyr columns, per gene; Tyr log-likelihood.
Masking: first/last 10 codons of each gene alignment (chained) and codons within 10 match columns of the profile ends."""
import sys,os,json,math
sys.path.insert(0,"src"); import panel_guard
import pyhmmer, numpy as np
from collections import defaultdict
HMMFILE="data/pfam/mito/metazoan_mito_proteins.hmm"
STD="KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
B="TCAG"; CODONS=[a+b+c for a in B for b in B for c in B]
# NCBI table 5 (invertebrate mito) in TCAG order
T5="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
TAB={c:T5[i] for i,c in enumerate(CODONS)}
for c in ("TAA","TAG"): TAB[c]="X"
AA="ACDEFGHIKLMNPQRSTVWY"; ABC=pyhmmer.easel.Alphabet.amino()
_HM=None
def hmms():
    global _HM
    if _HM is None:
        with pyhmmer.plan7.HMMFile(HMMFILE) as f: _HM=list(f)
    return _HM
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def load(fa):
    t=open(fa).read().split("\n",1); return t[0][1:].split()[0],"".join(t[1].split()).upper()
def frames(seq,wrap=0):
    s=seq+seq[:wrap] if wrap else seq; L=len(s); out=[]
    for strand,x in (("+",s),("-",rc(s))):
        for f in range(3):
            cod=[x[i:i+3] for i in range(f,L-2,3)]; out.append((strand,f,cod,x,"".join(TAB.get(c,"X") if set(c)<=set("ACGT") else "X" for c in cod)))
    return out
def analyse(fa,circular=True,E=1e-5,tyr_thr=0.5):
    name,seq=load(fa); L0=len(seq); wrap=2000 if circular else 0
    fr=frames(seq,wrap); seqs=[]
    for i,(st,f,c,x,aa) in enumerate(fr):
        seqs.append(pyhmmer.easel.TextSequence(name=f"f{i}",sequence=aa).digitize(ABC))
    H=hmms(); bg=None
    pl=pyhmmer.plan7.Pipeline(ABC,bg_filter=False,E=E,domE=E,incE=E,incdomE=E) if False else None
    doms=[]
    for hi,h in enumerate(H):
        res=next(pyhmmer.hmmsearch(h,seqs,cpus=1,E=E,domE=E,incE=E,incdomE=E,bias_filter=False,F1=1.0,F2=1.0,F3=1.0,Z=1,domZ=1))
        for hit in res:
            for d in hit.domains:
                if d.i_evalue>E: continue
                al=d.alignment; doms.append(dict(gene=str(h.name),frame=int(str(hit.name)[1:]),score=d.score,ie=d.i_evalue,tf=al.target_from,tt=al.target_to,hf=al.hmm_from,ht=al.hmm_to,hs=al.hmm_sequence,ts=al.target_sequence,h=hi))
    # nucleotide coordinates on the forward strand (0-based), then greedy best non-overlapping domains per locus
    for d in doms:
        st,f,cod,x,_=fr[d["frame"]]; d["strand"]=st; d["c0"],d["c1"]=d["tf"]-1,d["tt"]
        a=f+3*d["c0"]; b=f+3*d["c1"]
        d["ns"],d["ne"]=(a,b) if st=="+" else (len(x)-b,len(x)-a)
    doms.sort(key=lambda d:-d["score"]); kept=[]
    for d in doms:
        if any(k["strand"]==d["strand"] and min(d["ne"],k["ne"])-max(d["ns"],k["ns"])>30 for k in kept): continue
        if wrap and d["ns"]>=L0: continue
        kept.append(d)
    g=defaultdict(list)
    for d in kept: g[d["gene"]].append(d)
    me={hi:np.array(h.match_emissions) for hi,h in enumerate(H)}  # [M+1,K] row 0 = entry
    idx={a:ABC.symbols.index(a) for a in AA}; Yi=idx["Y"]
    recs=[]  # (gene,codon,col,pY,aaVec,masked)
    for gene,ds in g.items():
        M=H[ds[0]["h"]].M; 
        gs,ge=min(x["ns"] for x in ds),max(x["ne"] for x in ds)
        for d in ds:
            st,f,cod,x,_=fr[d["frame"]]; hp,tp=d["hf"]-1,d["tf"]-1
            for hc,tc in zip(d["hs"],d["ts"]):
                if hc!=".": hp+=1
                if tc!="-":
                    tp+=1
                if hc!="." and tc!="-" and hc.isalpha():
                    col=hp; codon=cod[tp-1]; ntpos=f+3*(tp-1); npos=ntpos if st=="+" else len(x)-ntpos-3
                    masked=(col<=10 or col>M-10 or npos<gs+30 or npos+3>ge-30)
                    recs.append((gene,codon,col,float(me[d["h"]][col,Yi]),d["h"],masked,npos,st))
    return name,L0,kept,recs,me,idx,H
def summarise(fa,circular=True,E=1e-5,thr=0.5):
    name,L0,kept,recs,me,idx,H=analyse(fa,circular,E)
    use=[r for r in recs if not r[5]]
    taa=[r for r in use if r[1]=="TAA"]; ntaa=len(taa)
    tyrcols=[r for r in use if r[3]>=thr]
    n_t=lambda c: sum(1 for r in tyrcols if r[1]==c)
    cnt={c:n_t(c) for c in ("TAA","TAT","TAC","TAG")}
    other=len(tyrcols)-sum(cnt.values())
    bgv=np.array(pyhmmer.plan7.Background(ABC).residue_frequencies); ll={}
    for a,i in idx.items():
        ll[a]=sum(math.log(max(me[r[4]][r[2],i],1e-6)/bgv[i]) for r in taa)
    llY=ll["Y"]; second=max(v for k,v in ll.items() if k!="Y") if taa else 0
    genes=sorted({r[0] for r in taa}); gY=sorted({r[0] for r in tyrcols if r[1]=="TAA"}); 
    tag=sum(1 for r in use if r[1]=="TAG")
    return dict(panel_id=name,length=L0,n_domains=len(kept),n_genes_hit=len({k["gene"] for k in kept}),aligned_codons=len(use),n_taa_internal=ntaa,n_taa_genes=len(genes),n_taa_tyrcol=cnt["TAA"],n_tyrcol_genes=len(gY),tyrcols=len(tyrcols),tat=cnt["TAT"],tac=cnt["TAC"],tag_tyrcol=cnt["TAG"],tyrcol_other=other,n_tag_internal=tag,llr_Y=round(llY,2),llr_Y_vs_2nd=round(llY-second,2) if taa else 0.0,frac_taa_tyr=round(cnt["TAA"]/max(1,cnt["TAA"]+cnt["TAT"]+cnt["TAC"]),4))
if __name__=="__main__":
    import csv
    fs=sys.argv[1:]; ids=[os.path.basename(f).replace(".fna","") for f in fs]; panel_guard.check(ids)
    topo={r["panel_id"]:r["topology"] for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
    w=None
    for i,f in zip(ids,fs):
        d=summarise(f,circular=topo.get(i,"circular")!="linear")
        if w is None: w=list(d); print("\t".join(w))
        print("\t".join(str(d[k]) for k in w),flush=True)
