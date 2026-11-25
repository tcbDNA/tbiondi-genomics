"""Cold review 2, check 4 supplement (POST-HOC): a second mitochondrial tRNA-Tyr that is too divergent for a body search (the Globodera case: the two genes share about 19 of 26 aligned bases) but sits on a high-copy molecule.
Reads are subsampled so the mitogenome is at ~TARGET x; at that depth nuclear single-copy sequence is mostly below 3x. SPAdes assembles the subsample; the tRNA search on high-copy contigs is step 2 (src/cr2_highcopy2.py).
Only informative where nuclear depth << mitochondrial depth (genome skims). M. kubotai and A. glossodonta use the contigs from before repeat resolution (--disable-rr; for M. kubotai K77/before_rr.fasta, because SPAdes stalled in distance estimation). usage: cr2_highcopy.py LABEL RUN MTFASTA(accession or path) TARGET"""
import sys,os,re,subprocess,random,collections
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import trny_extract as TE
lab,run,ref,target=sys.argv[1:5]; target=float(target); B="env/prefix/bin/"; M="env/tools/mapenv/bin/"; W=f"results/work/r2/hc/{lab}/"; os.makedirs(W,exist_ok=True); out=open(f"results/posthoc2/secondcopy/highcopy_{lab}.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=out,flush=True))
mt=ref if os.path.exists(ref) else "data/sweep/fa/"+str(sum(map(ord,ref))%200).zfill(3)+f"/{ref}.fna"; seq="".join(l.strip() for l in open(mt) if l[0]!=">").upper(); open(W+"mt.fa","w").write(">mt\n"+seq+"\n")
def depth(fq):
    subprocess.run(["bash","-c",f"{M}minimap2 -t 6 -ax sr {W}mt.fa {' '.join(fq)} 2>/dev/null | {M}samtools view -b -F 4 - | {M}samtools sort -@ 2 -o {W}d.bam - && {M}samtools index {W}d.bam"]); d=sorted(int(x.split()[2]) for x in subprocess.run([M+"samtools","depth","-aa",W+"d.bam"],capture_output=True,text=True).stdout.splitlines()); return d[len(d)//2]
if len(sys.argv)>5 and sys.argv[5]=="have":   # FASTQ already in the folder (streamed subset)
    fq=[f"{W}{run}_{m}.fastq" for m in (1,2)]; P("POST-HOC check 4 supplement",lab,run,"spots used",sum(1 for _ in open(fq[0]))//4,"(streamed subset, first spots of the run)")
else:
    # pilot: first 400k spots -> mitochondrial depth per spot
    subprocess.run([B+"fastq-dump","-X","400000","--split-3","-O",W,f"results/work/{run}.sra"],capture_output=True); fq=[f"{W}{run}_{m}.fastq" for m in (1,2) if os.path.exists(f"{W}{run}_{m}.fastq")] or [f"{W}{run}.fastq"]; pd=depth(fq); nsp=sum(1 for _ in open(fq[0]))//4
    info=subprocess.run([B+"vdb-dump","--info",f"results/work/{run}.sra"],capture_output=True,text=True).stdout; tot=int(re.search(r"SEQ\s*:\s*([\d,]+)",info)[1].replace(",",""))
    need=int(target/max(pd,0.1)*nsp); need=min(need,tot); P("POST-HOC check 4 supplement",lab,run,"spots",tot,"pilot",nsp,"spots -> mitochondrial median depth",pd,"| spots used",need,f"({100*need/tot:.1f}%) | estimated full-run mitochondrial depth",round(pd*tot/nsp))
    for f in fq: os.remove(f)
    subprocess.run([B+"fastq-dump","-N",str(max(1,(tot-need)//2)),"-X",str(max(1,(tot-need)//2)+need-1),"--split-3","-O",W,f"results/work/{run}.sra"],capture_output=True); fq=[f"{W}{run}_{m}.fastq" for m in (1,2) if os.path.exists(f"{W}{run}_{m}.fastq")] or [f"{W}{run}.fastq"]
md=depth(fq); P("  mitochondrial median depth in the subsample",md)
a=["-1",fq[0],"-2",fq[1]] if len(fq)==2 else ["-s",fq[0]]
if not os.path.exists(W+"sp/contigs.fasta"): subprocess.run([B+"spades.py","--only-assembler"]+(["--disable-rr"] if lab in("A_glossodonta","M_kubotai") else [])+["-k","33,55,77"]+a+["-o",W+"sp","-t","10","-m","20"],stdout=open(W+"sp.log","w"),stderr=subprocess.STDOUT)
cf=W+"sp/contigs.fasta"
if not os.path.exists(cf): P("  SPAdes failed; see",W+"sp.log"); sys.exit(1)
S={};c=None
for l in open(cf):
    if l[0]==">": c=l[1:].strip(); S[c]=[]
    else: S[c].append(l.strip())
S={k:"".join(v) for k,v in S.items()}; open(W+"all.fa","w").write("".join(f">{k}\n{v}\n" for k,v in S.items()))
subprocess.run(["bash","-c",f"{M}minimap2 -t 6 -ax sr {W}all.fa {' '.join(fq)} 2>/dev/null | {M}samtools view -b -F 2308 - | {M}samtools sort -@ 2 -o {W}c.bam - && {M}samtools index {W}c.bam"])
cov={}
for l in subprocess.run([M+"samtools","coverage",W+"c.bam"],capture_output=True,text=True).stdout.splitlines():
    if l[0]!="#": f=l.split("\t"); cov[f[0]]=float(f[6])
subprocess.run([B+"makeblastdb","-in",W+"mt.fa","-dbtype","nucl","-out",W+"mt"],capture_output=True); mtc=collections.defaultdict(int)
for l in subprocess.run([B+"blastn","-query",W+"all.fa","-db",W+"mt","-evalue","1e-10","-dust","no","-outfmt","6 qseqid pident length"],capture_output=True,text=True).stdout.splitlines(): k,i,n=l.split("\t"); mtc[k]+=int(n)*(float(i)>=90)
hi=[k for k in S if cov.get(k,0)>=0.2*md and len(S[k])>=200]; ism=[k for k in hi if mtc[k]>=0.5*len(S[k])]; P("  contigs",len(S),"| coverage >=20% of the mitogenome and >=200 bp:",len(hi),"total bp",sum(len(S[k]) for k in hi),"| of them mitochondrial record:",len(ism),sum(len(S[k]) for k in ism),"bp")
P("DONE step 1 (assembly and coverage); the tRNA search is step 2, src/cr2_highcopy2.py")
