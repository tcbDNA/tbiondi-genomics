"""Cold review 2, check 4 supplement, second step (POST-HOC): for each subsample assembly, every contig >=500 bp with coverage >=50% of the mitogenome's is searched alone with the locked extractor (cmsearch --max, so the E-value is per contig as in the sweep) and ARAGORN. usage: cr2_highcopy2.py LABEL..."""
import sys,os,subprocess,collections,re
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import trny_extract as TE; B="env/prefix/bin/"; M="env/tools/mapenv/bin/samtools"
from concurrent.futures import ProcessPoolExecutor
def one(j):
    W,k,s=j; t=f"{W}one_{os.getpid()}.fa"; open(t,"w").write(f">c\n{s}\n"); c=TE.call(t); r=subprocess.run([B+"aragorn","-mt","-gc5","-w",t],capture_output=True,text=True).stdout; os.remove(t); return k,c,[" ".join(l.split()) for l in r.splitlines() if "(gta)" in l.lower()]
if __name__=="__main__":
    o=open("results/posthoc2/secondcopy/highcopy_summary.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True)); P("POST-HOC check 4 supplement, step 2 (contigs >=500 bp at >=50% of the mitochondrial depth)")
    for lab in sys.argv[1:]:
        W=f"results/work/r2/hc/{lab}/"
        if not os.path.exists(W+"c.bam.bai"): P(lab,"not assembled"); continue
        S={};c=None
        for l in open(W+"all.fa"):
            if l[0]==">": c=l[1:].strip(); S[c]=[]
            else: S[c].append(l.strip())
        S={k:"".join(v) for k,v in S.items()}; cov={}
        for l in subprocess.run([M,"coverage",W+"c.bam"],capture_output=True,text=True).stdout.splitlines():
            if l[0]!="#": f=l.split("\t"); cov[f[0]]=float(f[6])
        d=sorted(int(x.split()[2]) for x in subprocess.run([M,"depth","-aa",W+"d.bam"],capture_output=True,text=True).stdout.splitlines()); md=d[len(d)//2]
        mtc=collections.defaultdict(int)
        for l in subprocess.run([B+"blastn","-query",W+"all.fa","-db",W+"mt","-evalue","1e-10","-dust","no","-outfmt","6 qseqid pident length"],capture_output=True,text=True).stdout.splitlines(): k,i,n=l.split("\t"); mtc[k]+=int(n)*(float(i)>=90)
        hi=[k for k in S if len(S[k])>=500 and cov.get(k,0)>=.5*md]; nm=[k for k in hi if mtc[k]<.5*len(S[k])]; mtbp=sum(len(S[k]) for k in hi if k not in nm)
        with ProcessPoolExecutor(6) as ex: R=list(ex.map(one,[(W,k,S[k]) for k in hi]))
        P(f"\n{lab}: mitochondrial depth {md}x; contigs at >=50% of it: {len(hi)} (mitogenome record {mtbp} bp; other {len(nm)} contigs, {sum(len(S[k]) for k in nm)} bp)")
        n=collections.Counter()
        for k,c,a in R:
            if c.get("status")!="called" and not a: n["no candidate"]+=1; continue
            u,v=sorted((c.get("start",0),c.get("end",0))); near=[z for z in a for q in [re.search(r"\[(-?\d+),(\d+)\]",z)] if q and min(v,int(q[2]))-max(u,int(q[1]))>0.5*min(v-u,int(q[2])-int(q[1]))]
            sc=float(c.get("score",0) or 0) if c.get("status")=="called" else 0; gate=sc>=20 or (sc>=10 and bool(near))
            if c.get("status")=="called" and c.get("N33")!="U": gate=sc>=10 and bool(near) and c.get("stem_pairs",0)>=4
            if k in nm and not gate: n["candidate fails the locked gate (score >=20, or >=10 with ARAGORN at the same place)"]+=1; continue
            P("   ","record" if k not in nm else "OTHER contig passing the locked gate",k[:44],"cov",round(cov[k],1),"| extractor:",(c.get("N33"),c.get("loop"),"score",c.get("score"),"stem",c.get("stem_pairs"),c.get("start"),c.get("end")) if c.get("status")=="called" else "none","| ARAGORN at the same place:",near or "none")
        P("    other contigs:",dict(n))
