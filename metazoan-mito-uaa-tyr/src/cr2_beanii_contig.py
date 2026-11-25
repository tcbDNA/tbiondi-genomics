"""Cold review 2, check 4 supplement (POST-HOC): what is NODE_189 in the S. beanii subsample assembly (a U33 GUA arm passing the locked gate)? k-mer coverage, divergence of the reads piled on it, similarity to mitochondrial genes. writes results/posthoc2/secondcopy/beanii_contig.txt"""
import subprocess,re,statistics as st
W="results/work/r2/hc/S_beanii/"; M="env/tools/mapenv/bin/samtools"; B="env/prefix/bin/"; o=open("results/posthoc2/secondcopy/beanii_contig.txt","w"); P=lambda *x:(print(*x),print(*x,file=o)); S={};c=None
for l in open(W+"all.fa"):
    if l[0]==">": c=l[1:].strip(); S[c]=""
    else: S[c]+=l.strip()
k=[x for x in S if x.startswith("NODE_189_")][0]; rec=[x for x in S if x.startswith("NODE_1_")][0]; P("POST-HOC",k,"sequence:",S[k])
def div(x):
    v=[];d=[];clip=[]
    for l in subprocess.run([M,"view","-F","2304",W+"c.bam",x],capture_output=True,text=True).stdout.splitlines():
        f=l.split("\t"); n=[int(z[5:]) for z in f[11:] if z.startswith("NM:i:")][0]; a=sum(int(q) for q,t in re.findall(r"(\d+)([M=X])",f[5])); s=sum(int(q) for q,t in re.findall(r"(\d+)([SH])",f[5])); v.append(n/max(1,a)); clip.append(s/(a+s)); d.append(int(f[4]))
    return len(v),round(st.median(v),4),round(st.mean(clip),2),sum(q>=20 for q in d)/len(d)
for x in (rec,k): P("  ",x,"reads",div(x)[0],"median divergence of reads from the contig",div(x)[1],"mean clipped fraction",div(x)[2],"share MAPQ>=20",round(div(x)[3],2),"| SPAdes k-mer coverage",x.split("cov_")[1])
open(W+"n189.fa","w").write(">n\n"+S[k]+"\n")
r=subprocess.run([B+"blastn","-task","blastn","-word_size","7","-query",W+"n189.fa","-db",W+"mt","-evalue","10","-dust","no","-outfmt","6 pident length qstart qend sstart send evalue"],capture_output=True,text=True).stdout; P("  blastn to the S. beanii mitogenome (E<=10):",[l.replace("\t"," ") for l in r.splitlines()[:4]])
tr=subprocess.run(["env/prefix/bin/python","-c",f"import sys;sys.path.insert(0,'src');import os;os.environ['GATE3_UNLOCK']='1';import uaa_caller as U;n,L,k,*_=U.analyse('{W}n189.fa',False);print([(x['gene'],x['ns'],x['ne']) for x in k])"],capture_output=True,text=True); P("  mitochondrial protein-gene domains on it (the sweep's profiles):",tr.stdout.strip() or tr.stderr[-200:])
d=[l.split() for l in subprocess.run([M,"depth","-aa","-r",k,W+"c.bam"],capture_output=True,text=True).stdout.splitlines()]; v=[int(x[2]) for x in d]; P("  read depth along it in 125-nt bins:",[round(st.mean(v[i:i+125])) for i in range(0,len(v),125)],"; candidate arm at 546-627")
for nm,l in (("ARAGORN other tRNAs on the contig",subprocess.run([B+"aragorn","-mt","-gc2","-w",W+"n189.fa"],capture_output=True,text=True).stdout),): P("  ",nm,[" ".join(x.split()) for x in l.splitlines() if "mtRNA" in x])
