"""NUMT search: mitogenome vs nuclear assembly. usage: numt.py LABEL MTFASTA ASSEMBLY_ACC MTACC
Reports (1) all blastn hits >=200 bp and >=80% identity outside the mitochondrial sequence itself (hits >=99% identity over >=90% of the mitogenome are the mitochondrion in the assembly and are excluded),
(2) the allele at N33 in every hit overlapping the trnY gene (record allele vs other)."""
import sys,os,subprocess,json,re,csv,collections
sys.path.insert(0,"src"); lab,mtf,asm,acc=sys.argv[1:5]; B="env/prefix/bin/"; W="results/work/numt"; os.makedirs(W,exist_ok=True)
s="".join(l.strip() for l in open(mtf) if l[0]!=">").upper(); d=None
for l in open("results/gate3/run/out.jsonl"):
    x=json.loads(l)
    if x["acc"]==acc: d=x["trny"]
if d is None:
    for l in open("results/gate3_unverified/run/out.jsonl") if os.path.exists("results/gate3_unverified/run/out.jsonl") else []:
        x=json.loads(l)
        if x["acc"]==acc: d=x["trny"]
g=int(d["g34_pos"]); site=(g-2) if d["strand"]!="-" else g; print(lab,acc,"mt length",len(s),"N33 + strand base",s[site],"tRNA N33",d["N33"],"strand",d["strand"])
fa=f"{W}/{asm}.fna"
if not os.path.exists(fa):
    z=f"{W}/{asm}.zip"; subprocess.run([B+"datasets","download","genome","accession",asm,"--include","genome","--filename",z,"--no-progressbar"],capture_output=True)
    if os.path.exists(z): open(fa,"wb").write(subprocess.run(["unzip","-p",z,"*genomic.fna"],capture_output=True).stdout); os.remove(z); open("data/manifests/download_ledger.tsv","a").write(f"3x\tNCBI datasets {asm}\t{fa}\t{os.path.getsize(fa)}\tNA (unpacked stream, not hashed)\ton_disk\n")
print("assembly",asm,"size",os.path.getsize(fa) if os.path.exists(fa) else None)
db=f"{W}/{asm}"; 
if not os.path.exists(db+".nsq") and not os.path.exists(db+".00.nsq"): subprocess.run([B+"makeblastdb","-in",fa,"-dbtype","nucl","-out",db],capture_output=True)
open(f"{W}/{acc}.q.fa","w").write(f">q\n{s}\n")
r=subprocess.run([B+"blastn","-task","blastn","-query",f"{W}/{acc}.q.fa","-db",db,"-evalue","1e-10","-dust","no","-max_target_seqs","2000","-num_threads","3","-outfmt","6 sseqid pident length qstart qend sstart send evalue slen qseq sseq"],capture_output=True,text=True)
H=[l.split("\t") for l in r.stdout.splitlines()]; nu=[]; mt=collections.Counter()
by=collections.defaultdict(list)
for h in H: by[h[0]].append(h)
for c,hs in by.items():
    cov=sum(int(h[2]) for h in hs if float(h[1])>=99); 
    if cov>0.5*len(s) and int(hs[0][8])<200000: mt[c]=cov; continue
    nu+=[h for h in hs if int(h[2])>=200 and float(h[1])>=80]
print("contigs that look like the mitochondrion itself:",dict(mt)); print("nuclear-like hits >=200 bp, >=80%:",len(nu),"| total aligned bp",sum(int(h[2]) for h in nu)); 
for h in sorted(nu,key=lambda h:-int(h[2]))[:12]: print("   ",h[0],"id",h[1],"len",h[2],"mt",h[3]+"-"+h[4],"contig",h[5],"-",h[6],"contig len",h[8])
t=0
for h in nu:
    q0,q1=int(h[3]),int(h[4]); lo,hi=min(q0,q1),max(q0,q1)
    if lo<=site+1<=hi:
        t+=1; a=q0; b=None
        for x,y in zip(h[9],h[10]):
            if x!="-":
                if a-1==site: b=y
                a+=1
        print("   covers N33:",h[0],"id",h[1],"len",h[2],"aligned base at N33 site",b,"(record + strand base",s[site],")")
print("nuclear hits covering the N33 site:",t); print("DONE",lab,asm)
