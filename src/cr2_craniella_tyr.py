"""Cold review 2, check 4 supplement (POST-HOC): what are the high-copy C. infrequens windows whose Tyr(GUA) U33 candidates pass the locked gate? tRNAscan-SE in eukaryotic and organellar mode, identity of the candidate to cytosolic tRNA-Tyr and to the mitochondrial tRNA-Tyr, number of copies, other mitochondrial content (the 13 protein profiles, blastn to the mitogenome). writes results/posthoc2/secondcopy/craniella_tyr_contigs.txt"""
import sys,os,subprocess,re,collections
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import trny_extract as TE, uaa_caller as U, cr2_secondcopy as SC; B="env/prefix/bin/"; W="results/work/r2/cran/"; o=open("results/posthoc2/secondcopy/craniella_tyr_contigs.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True)); P("POST-HOC")
S={};c=None
for l in open("results/work/numt/GCA_974939545.1.fna"):
    if l[0]==">": c=l[1:].split()[0]; S[c]=[]
    else: S[c].append(l.strip())
S={k:"".join(v).upper() for k,v in S.items()}; mt=SC.body("OZ318138.1")[0]; cyt=SC.CYT[1]; A=SC.al(); rc=SC.rc; seen=collections.Counter()
def ident(a,b):
    x=A.align(a,b)[0]; m=n=0
    for (u,v),(p,q) in zip(*x.aligned): m+=sum(a[u+i]==b[p+i] for i in range(v-u)); n+=v-u
    return f"{100*m/n:.0f}% over {n} nt"
for ct,a,b in (("CDFYFP010000042.1",70001,85000),("CDFYFP010000083.1",65001,85000),("CDFYFP010000180.1",20001,30000),("CDFYFP010000180.1",35001,50000),("CDFYFP010000288.1",25001,30000),("CDFYFP010000419.1",15001,20000),("CDFYFP010000449.1",1,24249),("CDFYFP010000849.1",1,15000)):
    w=S[ct][a-1:b]; t=W+"tyr.fa"; open(t,"w").write(f">w\n{w}\n"); h,e=TE.run(t,True); ok=[x for x in h if x["anticodon"]=="GUA" and x["loop_ok"] and x["stem_pairs"]>=3 and x["score"]>=10]
    euk=[l.split() for l in subprocess.run([B+"tRNAscan-SE","-E","-Q","-q","--brief",t],capture_output=True,text=True,stdin=subprocess.DEVNULL).stdout.splitlines()]; ey=[(x[2],x[3],x[4],x[5],x[8]) for x in euk if len(x)>=9 and x[5].upper()=="GTA"]
    try: k=U.analyse(t,False)[2]; g=sorted({x["gene"] for x in k})
    except Exception as ex_: g=str(ex_)[:60]
    bl=subprocess.run([B+"blastn","-task","blastn","-query",t,"-subject",W+"OZ318138.1.q","-evalue","1e-5","-dust","no","-outfmt","6 pident length"],capture_output=True,text=True).stdout.splitlines()
    P(f"\n{ct} {a}-{b} (contig {len(S[ct])} bp): CM Tyr(GUA) candidates score>=10: {len(ok)}, loops {dict(collections.Counter(x['loop'] for x in ok))}; tRNAscan-SE eukaryotic mode Tyr(GTA): {len(ey)} {ey[:4]}; all tRNAs in eukaryotic mode {len(euk)}; mitochondrial protein domains {g}; blastn to the mitogenome (E<=1e-5): {len(bl)} {bl[:2]}")
    for x in ok[:3]:
        u,v=sorted((x["start"],x["end"])); q=w[u-1:v]; q=rc(q) if x["strand"]=="-" else q; P("    ",q,"| to cytosolic tRNA-Tyr",ident(q,cyt),"| to the mitochondrial tRNA-Tyr",ident(q,mt),"score",x["score"])
