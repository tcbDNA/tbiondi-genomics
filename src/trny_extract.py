"""Anticodon-anchored trnY extractor. Reads canonical N33 (base 5' of verified G34-U35-A36) from MITOS trnY.cm consensus columns.
Usage: trny_extract.py fasta [fasta...] > tsv  (control panel only unless GATE3_UNLOCK=1 and PREREG.md committed)"""
import sys,subprocess,re,os,tempfile
sys.path.insert(0,"src"); import panel_guard
CM="data/mitos_ref/refseq89m/ncRNA/trnY.cm"; CMS="env/prefix/bin/cmsearch"
# consensus columns (1-based, verified against model consensus): stem5' 25-29, loop 30-36 (N32..N38), stem3' 37-41
C33,C34,C35,C36,S5,S3,LOOP=31,32,33,34,range(25,30),range(37,42),range(30,37)
PAIR={("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def parse_sto(path):
    seqs={};rf="";cur=[]
    for l in open(path):
        if l.startswith("#=GC RF"): rf+=l.split(None,2)[2].strip()
        elif l.startswith("#=GC") or l.startswith("#=GR") or l.startswith("#=GS") or l.startswith("#=GF") or l.startswith("#") or l.strip() in ("","//"): continue
        else:
            a,b=l.split(None,1); seqs[a]=seqs.get(a,"")+b.strip()
    return seqs,rf
def col_map(rf):
    m={};c=0
    for i,ch in enumerate(rf):
        if ch not in ".~": c+=1; m[c]=i
    return m
def run(fa,maxsens=False):
    with tempfile.TemporaryDirectory() as d:
        sto=f"{d}/a.sto"; tb=f"{d}/t"
        cmd=[CMS,"--notextw","-E","10","--incE","10","--noali","-A",sto,"--tblout",tb,"--cpu","2"]+(["--max"] if maxsens else [])+[CM,fa]
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
        if r.returncode: return None,("cmsearch_error",r.stderr[:150])
        tbl={}
        for l in open(tb):
            if l.startswith("#"): continue
            f=l.split(); tbl[f"{f[0]}/{f[7]}-{f[8]}"]=(int(f[7]),int(f[8]),f[9],float(f[14]),float(f[15]))
        if not os.path.exists(sto): return [],None
        seqs,rf=parse_sto(sto); cm=col_map(rf); out=[]
        for nm,al in seqs.items():
            key=nm if nm in tbl else next((k for k in tbl if k.split("/")[-1]==nm.split("/")[-1] and k.split("/")[0]==nm.split("/")[0]),None)
            s,e,strand,score,ev=tbl[key] if key else (0,0,"?",0,0)
            g=lambda c: al[cm[c]].upper().replace("T","U") if c in cm else "-"
            idx=[cm[c] for c in LOOP]; ins=any(al[i].islower() or False for i in range(idx[0],idx[-1]+1))
            contig=idx[-1]-idx[0]==6 and all(al[i] not in "-." for i in idx)
            anti="".join(g(c) for c in (C34,C35,C36)); stem=sum((g(a),g(b)) in PAIR for a,b in zip(S5,reversed(list(S3))))
            pre=sum(1 for ch in al[:cm[C34]] if ch not in "-.~"); g34=(s+pre if strand!="-" else s-pre) if cm.get(C34) is not None and al[cm[C34]] not in "-." else 0
            out.append(dict(hit=nm,start=s,end=e,strand=strand,score=score,evalue=ev,g34_pos=g34,anticodon=anti,N33=g(C33),loop="".join(g(c) for c in LOOP),stem_pairs=stem,loop_ok=contig,g37=g(35)))
        return out,None
def call(fa):
    best=None; tried=[]
    for sens in (False,True):
        hits,err=run(fa,sens)
        if err: return dict(status=err[0],reason=err[1])
        ok=[h for h in hits if h["anticodon"]=="GUA" and h["loop_ok"] and h["stem_pairs"]>=3]
        if ok:
            ok.sort(key=lambda h:-h["score"]); b=ok[0]; return dict(status="called",reason="sens" if sens else "std",n_hits=len(hits),n_ok=len(ok),**b)
        tried.append(hits)
    h=sorted([x for t in tried for x in t],key=lambda h:-h["score"])
    if not h: return dict(status="no_hit",reason="no_cm_hit_E<=10")
    b=h[0]; r="anticodon_not_GUA" if b["anticodon"]!="GUA" else "loop_not_7nt" if not b["loop_ok"] else "stem_lt3"
    return dict(status="failed",reason=r,n_hits=len(h),n_ok=0,**b)
if __name__=="__main__":
    ids=[os.path.basename(f).replace(".fna","") for f in sys.argv[1:]]; panel_guard.check(ids)
    cols="status reason n_hits n_ok hit start end strand score evalue anticodon N33 loop stem_pairs loop_ok g37".split()
    print("panel_id\t"+"\t".join(cols))
    for i,f in zip(ids,sys.argv[1:]):
        r=call(f); print(i+"\t"+"\t".join(str(r.get(c,"")) for c in cols),flush=True)
