"""Diagnostic variants of the extractor, panel only. usage: trny_variants.py VARIANT CMFILE outfile workers [idfile]
A = loop-check bug fix (only this sequence's own inserted residues count); B = A plus glocal (-g); C = B with another CM.
All hits with a verified GUA arm are kept so that gates can be compared afterwards."""
import sys,os,json,subprocess,tempfile,csv
from concurrent.futures import ProcessPoolExecutor,as_completed
sys.path.insert(0,"src"); import panel_guard, trny_extract as T
def hits(fa,cm,glocal,maxs,cols):
    with tempfile.TemporaryDirectory() as d:
        sto,tb=f"{d}/a",f"{d}/t"
        c=[T.CMS,"--notextw","-E","10","--incE","10","--noali","-A",sto,"--tblout",tb,"--cpu","1"]+(["-g"] if glocal else [])+(["--max"] if maxs else [])+[cm,fa]
        r=subprocess.run(c,capture_output=True,text=True,timeout=900)
        if r.returncode: raise RuntimeError("cmsearch "+(r.stderr or r.stdout)[-150:])
        if not os.path.exists(sto): return []
        tbl={}
        for l in open(tb):
            if l[0]!="#": f=l.split(); tbl[f"{f[0]}/{f[7]}-{f[8]}"]=(int(f[7]),int(f[8]),f[9],float(f[14]),float(f[15]))
        S,rf=T.parse_sto(sto); m=T.col_map(rf); out=[]
        C33,C34,S5,S3,LOOP=cols
        for nm,al in S.items():
            if nm not in tbl: continue
            s,e,st,sc,ev=tbl[nm]; g=lambda k: al[m[k]].upper().replace("T","U") if k in m else "-"
            a,b=m[LOOP[0]],m[LOOP[-1]]; lc={m[k] for k in LOOP}
            own=[al[i] for i in range(a,b+1) if i not in lc and al[i] not in ".-~"]   # this sequence's own insertions only
            ok=not own and all(al[m[k]] not in ".-~" for k in LOOP); anti="".join(g(C34+i) for i in range(3))
            stem=sum((g(x),g(y)) in T.PAIR for x,y in zip(S5,reversed(list(S3))))
            pre=sum(1 for ch in al[:m[C34]] if ch not in "-.~"); pos=s+pre if st!="-" else s-pre
            if anti=="GUA" and ok and stem>=3: out.append(dict(start=s,end=e,strand=st,score=sc,evalue=ev,g34_pos=pos,N33=g(C33),loop="".join(g(k) for k in LOOP),stem_pairs=stem))
        return sorted(out,key=lambda h:-h["score"])
MITOS=(31,32,range(25,30),range(37,42),range(30,37))
def one(a):
    pid,var,cm,cols=a; fa=f"data/panel/{pid}.fna"
    try:
        for mx in (False,True):
            h=hits(fa,cm,var!="A",mx,cols)
            if h: return dict(panel_id=pid,status="called",reason="sens" if mx else "std",n_ok=len(h),**h[0])
        return dict(panel_id=pid,status="no_hit",reason="no_verified_GUA_arm")
    except Exception as e: return dict(panel_id=pid,status="exception",reason=str(e)[:200])
if __name__=="__main__":
    var,cm,out,n=sys.argv[1:5]; P=[r["panel_id"] for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")]
    ids=[l.strip() for l in open(sys.argv[5])] if len(sys.argv)>5 else P; panel_guard.check(ids)
    cols=MITOS if len(sys.argv)<7 else eval(sys.argv[6])
    with ProcessPoolExecutor(int(n)) as ex,open(out,"w") as o:
        for k,fu in enumerate(as_completed([ex.submit(one,(i,var,cm,cols)) for i in ids])): o.write(json.dumps(fu.result())+"\n"); o.flush()
    print(var,"done",len(ids))
