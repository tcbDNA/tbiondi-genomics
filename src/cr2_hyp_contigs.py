"""Cold review 2, check 4 supplement (POST-HOC): identify high-copy contigs whose U33 Tyr(GUA) candidates pass the locked gate. For each: k-mer coverage against the mitogenome contig, GC, read divergence, the 13 mitochondrial protein profiles, blastn to the mitogenome record, tRNAscan-SE (bacterial and organellar mode), Pfam domains on six-frame ORFs (table 11).
usage: cr2_hyp_contigs.py LABEL NODE_x[,NODE_y...] -> results/posthoc2/secondcopy/contig_id_LABEL.txt"""
import sys,os,subprocess,re,collections,statistics as st
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import uaa_caller as U, pyhmmer; B="env/prefix/bin/"; M="env/tools/mapenv/bin/samtools"
lab,nodes=sys.argv[1],sys.argv[2].split(","); W=f"results/work/r2/hc/{lab}/"; o=open(f"results/posthoc2/secondcopy/contig_id_{lab}.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True)); P("POST-HOC",lab); S={};c=None
for l in open(W+"all.fa"):
    if l[0]==">": c=l[1:].strip(); S[c]=[]
    else: S[c].append(l.strip())
S={k:"".join(v) for k,v in S.items()}; rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
Bs="TCAG"; T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; AA={a+b+c:T[16*Bs.index(a)+4*Bs.index(b)+Bs.index(c)] for a in Bs for b in Bs for c in Bs}
def orfs(s):
    for st_,x in (("+",s),("-",rc(s))):
        for f in range(3):
            p="".join(AA.get(x[i:i+3],"X") for i in range(f,len(x)-2,3))
            for m in re.finditer(r"[^*]{60,}",p): yield m[0]
abc=pyhmmer.easel.Alphabet.amino(); hm=[f for f in os.listdir("data/pfam") if f.endswith(".hmm") or f.endswith(".h3m")]; hmm="data/pfam/"+([f for f in hm if "Pfam-A" in f] or hm or [""])[0]
for n in nodes:
    k=[x for x in S if x.startswith(n+"_length")][0]; s=S[k]; t=W+"id.fa"; open(t,"w").write(f">c\n{s}\n")
    v=[]
    for l in subprocess.run([M,"view","-F","2304",W+"c.bam",k],capture_output=True,text=True).stdout.splitlines()[:4000]:
        f=l.split("\t"); nm=[int(z[5:]) for z in f[11:] if z.startswith("NM:i:")]; a=sum(int(q) for q,y in re.findall(r"(\d+)([M=X])",f[5])); v.append(nm[0]/max(1,a) if nm else 0)
    bl=subprocess.run([B+"blastn","-task","blastn","-query",t,"-db",W+"mt","-evalue","1e-3","-dust","no","-outfmt","6 pident length"],capture_output=True,text=True).stdout.splitlines()
    try: g=sorted({x["gene"] for x in U.analyse(t,False)[2]})
    except Exception as e: g=str(e)[:60]
    bt=[l.split() for l in subprocess.run([B+"tRNAscan-SE","-B","-Q","-q","--brief",t],capture_output=True,text=True,stdin=subprocess.DEVNULL).stdout.splitlines()]
    P(f"\n{k}: GC {100*sum(s.count(x) for x in 'GC')/len(s):.0f}% | median read divergence {st.median(v) if v else None} (n={len(v)}) | mitochondrial protein domains {g} | blastn to the mitogenome: {bl[:2]} | tRNAscan-SE bacterial mode: {[(x[4],x[5],x[2],x[3],x[8]) for x in bt if len(x)>=9]}")
    O=[pyhmmer.easel.TextSequence(name=f"o{i}".encode(),sequence=p).digitize(abc) for i,p in enumerate(orfs(s))]; hits=collections.Counter()
    try:
        with pyhmmer.plan7.HMMFile(hmm) as H:
            for r in pyhmmer.hmmsearch(H,O,cpus=4,E=1e-5):
                for h in r:
                    if h.included: hits[str(r.query.name.decode() if hasattr(r.query.name,'decode') else r.query.name)]=max(hits[str(r.query.name.decode() if hasattr(r.query.name,'decode') else r.query.name)],round(h.score))
    except Exception as e: P("   Pfam error",str(e)[:120])
    P("   ORFs >=60 aa:",len(O),"| Pfam domains (score):",hits.most_common(14))
