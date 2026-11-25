"""Cold review 2, check 4 (POST-HOC): second-copy check. Reads that resemble the tRNA-Tyr body (no flanks, relaxed identity; src/tyscan.c) are aligned to tRNA-Tyr and, over the same read segment, to every other mitochondrial tRNA of that genome on both strands and to
cytosolic tRNA-Tyr. A read counts only if tRNA-Tyr wins outright, the alignment covers 33-10..33+10 and >=40 nt; the bases aligned to 34-36 are checked. usage: cr2_secondcopy.py scan ACC RUN | classify ACC RUN | simulate ACC | all"""
import sys,os,json,csv,subprocess,collections,random,hashlib
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,"src"); B="env/prefix/bin/"; W="results/work/r2/sc/"; O="results/posthoc2/secondcopy/"; os.makedirs(W,exist_ok=True); os.makedirs(O,exist_ok=True)
RC=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(RC)[::-1]
CYT=("cytosolic_Tyr_GUA","CCTTCGATAGCTCAGCTGGTAGAGCGGAGGACTGTAGATCCTTAGGTCGCTGGTTCAATTCCGGCTCGAAGGA")   # human nuclear tRNA-Tyr(GUA), intron removed
JOBS=[("NC_037360.1","SRR39401262"),("NC_037360.1","SRR4556446"),("CM165726.1","SRR38145622"),("CM165726.1","SRR13036352"),("NC_013627.1","SRR36227421"),("NC_013436.1","SRR22396609"),("NC_013628.2","SRR38094758"),("CM163399.1","SRR37037826"),("OZ318138.1","ERR14835550"),("OZ318138.1","SRR33975184"),("OZ318138.1","ERR14792883"),("OZ407450.1","ERR14842678"),("OZ407450.1","ERR14857142")]
def fasta(acc): return "".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna") if l[0]!=">").upper()
_o={}
def call(acc):
    if not _o:
        want={a for a,_ in JOBS}
        for l in open("results/gate3/run/out.jsonl"):
            x=json.loads(l)
            if x["acc"] in want: _o[x["acc"]]=x
    return _o[acc]
_f={}
def feats(acc):
    if not _f:
        want={a for a,_ in JOBS}
        for l in open("data/sweep/features_fixed.jsonl"):
            x=json.loads(l)
            if x["acc"] in want: _f[x["acc"]]=x["f"]
    return _f.get(acc)
def body(acc):
    s=fasta(acc); t=call(acc)["trny"]; a,b=sorted((int(t["start"]),int(t["end"]))); q=s[a-1:b]; q=q if t["strand"]!="-" else rc(q); g=int(t["g34_pos"]); p34=(g-a) if t["strand"]!="-" else (b-g); assert q[p34:p34+3]=="GTA",(acc,q,p34); return q,p34-1,s,a,b,t["strand"]
def others(acc):
    """other mitochondrial tRNAs: annotated ones plus ARAGORN calls; any interval overlapping tRNA-Tyr by >50% on the same strand is tRNA-Tyr itself and left out"""
    q,p,s,a,b,st=body(acc); f=feats(acc) or {}; L=[(x[0],x[1],x[2],str(x[3])) for x in f.get("trna",[])]; L+=[(x[0],x[1],x[2],"ARAGORN_"+x[3]+"_"+x[4]) for x in call(acc)["trna"]["aragorn"]]; out=[];seen=set()
    for u,v,z,n in L:
        u,v=min(u,v),max(u,v)
        if v-u>120 or u<1: continue
        if min(b,v)-max(a,u)>0.5*min(b-a,v-u) and z==st: continue
        x=s[u-1:v]; x=x if z!="-" else rc(x)
        if x in seen: continue
        seen.add(x); out.append((n,x)); out.append((n+"(antisense)",rc(x)))
    out.append(("Tyr(antisense)",rc(q))); out.append(CYT); out.append((CYT[0]+"(antisense)",rc(CYT[1]))); return out
def al():
    from Bio.Align import PairwiseAligner
    A=PairwiseAligner(); A.mode="local"; A.match_score=2; A.mismatch_score=-2; A.open_gap_score=-4; A.extend_gap_score=-1; return A
def scan(acc,run):
    fa=W+acc+".trnY.fa"; open(fa,"w").write(f">{acc}\n{body(acc)[0]}\n"); out=W+f"{acc}.{run}.win"
    if os.path.exists(out) and os.path.getsize(out)>0: return out
    subprocess.run(["bash","-c",f"{B}fastq-dump -Z --fasta 0 --split-spot results/work/{run}.sra 2>/dev/null | env/tools/tyscan/tyscan {fa} 46 > {out}.tmp 2> {out}.log && mv {out}.tmp {out}"]); return out
def read(fn):
    h=None
    for l in open(fn):
        if l[0]==">": h=l[1:].strip()
        else: yield h,l.strip()
_g={}
def setup(acc):
    if acc not in _g: q,p,s,a,b,st=body(acc); fl=s[max(0,a-1-150):b+150]; fl=fl if st!="-" else rc(fl); _g[acc]=(al(),q,p,others(acc),fl)
    return _g[acc]
def one(j):
    acc,seq=j; A,q,p,oth,fl=setup(acc)
    x=A.align(q,seq)
    try: x=x[0]
    except Exception: return None
    qi,si=x.aligned; m=n=0; mp={}
    for (u,v),(c,d) in zip(qi,si):
        for k in range(v-u): mp[u+k]=c+k; n+=1; m+=q[u+k]==seq[c+k]
    q0,q1=qi[0][0],qi[-1][1]; s0,s1=si[0][0],si[-1][1]; ncol=max(q1-q0,s1-s0); ident=m/ncol
    cover=q0<=p-10 and q1>=p+11
    seg=seq[max(0,s0-3):s1+3]; best=max(((A.score(o,seg),nm) for nm,o in oth),default=(0,"")); 
    b33=seq[mp[p]] if p in mp else "-"; ac="".join(seq[mp[p+k]] if p+k in mp else "-" for k in(1,2,3)); contig=all(p+k in mp for k in range(-1,6)) and all(mp[p+k]-mp[p]==k for k in range(-1,6))
    loop="".join(seq[mp[p+k]] if p+k in mp else "-" for k in range(-1,6))
    # flank agreement with the mitogenome outside the tRNA body
    y=A.align(fl,seq); fm=0; fn=0
    try:
        y=y[0]
        for (u,v),(c,d) in zip(*y.aligned):
            for k in range(v-u):
                if not(s0<=c+k<s1): fn+=1; fm+=fl[u+k]==seq[c+k]
        flank=(round(fm/fn,2) if fn else 0.0,fn)
    except Exception: flank=(0.0,0)
    return dict(score=x.score,ident=ident,n=ncol,cover=cover,win=x.score>best[0],rival=best[1] if best[0]>=x.score else "",b33=b33.replace("T","U"),ac=ac.replace("T","U"),contig=contig,loop=loop.replace("T","U"),flank=flank)
def band(i): return ">=97" if i>=.97 else "90-97" if i>=.90 else "80-90" if i>=.80 else "70-80" if i>=.70 else "<70"
def classify(acc,run,wfile=None,ex=None):
    fn=wfile or W+f"{acc}.{run}.win"; J=list(read(fn)); uniq=collections.Counter(s for _,s in J); U=list(uniq)
    res=list(ex.map(one,[(acc,s) for s in U],chunksize=200)) if ex else [one((acc,s)) for s in U]; D=dict(zip(U,res)); return J,uniq,D
def report(acc,run,ex):
    J,uniq,D=classify(acc,run,ex=ex); rec=call(acc)["trny"]["N33"]; tot=collections.Counter(); rows=[]; rival=collections.Counter(); T=collections.Counter(); oddac=collections.Counter(); susp=[]
    for s,k in uniq.items():
        d=D[s]
        if d is None or not d["cover"] or d["n"]<40: tot["not spanning or <40 nt"]+=k; continue
        if d["ident"]<.70: tot["identity <70%"]+=k; continue
        if not d["win"]: tot["another tRNA wins or ties"]+=k; rival[d["rival"]]+=k; continue
        if not d["contig"]: tot["gap in the loop"]+=k; continue
        if d["ac"]!="GUA": tot["anticodon not GUA"]+=k; oddac[(d["ac"],d["b33"],band(d["ident"]),"mt flank" if d["flank"][1]>=20 and d["flank"][0]>=.9 else "other/short flank")]+=k; continue
        tot["tRNA-Tyr(GUA)"]+=k; fl="mt flank" if d["flank"][1]>=20 and d["flank"][0]>=.9 else "no flank in read" if d["flank"][1]<20 else "other flank"
        T[(d["b33"],band(d["ident"]),fl)]+=k
        if d["b33"]!=rec: susp.append((k,d["b33"],round(d["ident"],3),d["loop"],fl,d["flank"],s))
    with open(O+f"{acc}.{run}.txt","w") as o:
        print("POST-HOC check 4",acc,run,"record N33",rec,"| scanner:",open(W+f"{acc}.{run}.win.log").read().strip(),file=o); print("windows",len(J),dict(tot),file=o); print("rivals",rival.most_common(8),file=o); print("non-GUA anticodons",oddac.most_common(8),file=o)
        for k,v in sorted(T.items(),key=lambda x:-x[1]): print("  N33=%s identity %s %s: %d"%(*k,v),file=o)
        print("reads with another base at 33:",len(susp),file=o); [print("   ",x,file=o) for x in sorted(susp,key=lambda x:-x[0])[:40]]
    by=collections.Counter(); byf=collections.Counter()
    for (b,i,f),v in T.items(): by[b]+=v; byf[(b,i)]+=v
    return dict(acc=acc,run=run,record=rec,windows=len(J),gua=tot["tRNA-Tyr(GUA)"],A=by["A"],C=by["C"],G=by["G"],U=by["U"],U_hi=byf[("U",">=97")],U_mid=byf[("U","90-97")],U_lo=byf[("U","80-90")]+byf[("U","70-80")],rival=tot["another tRNA wins or ties"],rivals=";".join(f"{a}:{b}" for a,b in rival.most_common(3)),nonGUA=tot["anticodon not GUA"],U_mtflank=sum(v for (b,i,f),v in T.items() if b=="U" and f=="mt flank"),U_otherflank=sum(v for (b,i,f),v in T.items() if b=="U" and f=="other flank"))
def simulate(acc,ex):
    q,p,*_=body(acc); rng=random.Random(7); out=[]
    for lev in (1.0,.95,.90,.85,.80,.75,.70):
        fn=W+f"{acc}.sim{int(lev*100)}.fa"; n=2000
        with open(fn,"w") as o:
            for i in range(n):
                s=list(q); sites=[k for k in range(len(q)) if not p<=k<=p+3]; 
                for k in rng.sample(sites,round((1-lev)*len(q))): s[k]=rng.choice([c for c in "ACGT" if c!=s[k]])
                s[p]="T"; s="".join(s); off=rng.randint(0,150-len(s)); r="".join(rng.choice("ACGT") for _ in range(off))+s+"".join(rng.choice("ACGT") for _ in range(150-len(s)-off)); r=r if rng.random()<.5 else rc(r); o.write(f">r{i}\n{r}\n")
        w=fn+".win"; subprocess.run(["bash","-c",f"env/tools/tyscan/tyscan {W}{acc}.trnY.fa 46 < {fn} > {w} 2>/dev/null"]); J,uniq,D=classify(acc,"",w,ex); ids={h.split(":")[0] for h,s in J if D[s] and D[s]["cover"] and D[s]["n"]>=40 and D[s]["ident"]>=.70 and D[s]["win"] and D[s]["contig"] and D[s]["ac"]=="GUA" and D[s]["b33"]=="U"}; out.append((lev,len(ids)/n))
    return out
if __name__=="__main__":
    m=sys.argv[1]
    if m=="scan": print(scan(sys.argv[2],sys.argv[3]))
    elif m=="others": [print(x) for x in others(sys.argv[2])]
    elif m=="all":
        with ProcessPoolExecutor(10) as ex:
            R=[report(a,r,ex) for a,r in JOBS if os.path.exists(W+f"{a}.{r}.win")]; w=csv.DictWriter(open(O+"summary.tsv","w"),fieldnames=list(R[0]),delimiter="\t"); w.writeheader(); w.writerows(R)
            with open(O+"sensitivity.tsv","w") as o:
                o.write("acc\tbody_identity\trecovered_as_tRNA-Tyr_U33\n")
                for a in dict.fromkeys(a for a,_ in JOBS):
                    if os.path.exists(W+a+".trnY.fa"): [o.write(f"{a}\t{l}\t{v:.3f}\n") for l,v in simulate(a,ex)]
