"""Cold review 2, check 2 (POST-HOC): are the Craniella G33 reads mitochondrial? Uses the competitive alignments from src/cr2_craniella.sh (each run against its own assembly, mitogenome and nuclear contigs together). writes results/posthoc2/craniella_depth.txt"""
import subprocess,re,json,collections,sys,statistics as st
M="env/tools/mapenv/bin/samtools"; W="results/work/r2/cran/"; o=open("results/posthoc2/craniella_depth.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True))
def run(*a): return subprocess.run([M]+list(a),capture_output=True,text=True).stdout
def cig(c): return [(int(n),k) for n,k in re.findall(r"(\d+)([MIDNSHP=X])",c)]
def base_at(f,pos):
    p=int(f[3]); i=0
    for n,k in cig(f[5]):
        if k in"M=X":
            if p<=pos<p+n: return f[9][i+pos-p]
            p+=n;i+=n
        elif k in"IS": i+=n
        elif k in"DN":
            if p<=pos<p+n: return "-"
            p+=n
T={}
for l in open("results/gate3/run/out.jsonl"):
    x=json.loads(l)
    if x["acc"] in("OZ318138.1","OZ407450.1"): T[x["acc"]]=x
J=[("infr_hifi","ERR14835550","OZ318138.1","GCA_974939545.1","DNA HiFi, own specimen"),("zetl_hifi","ERR14842678","OZ407450.1","GCA_978397455.1","DNA HiFi, own specimen"),("infr_rna","SRR33975184","OZ318138.1","GCA_974939545.1","RNA-seq, independent specimen"),("zetl_rna","ERR14857142","OZ407450.1","GCA_978397455.1","RNA-seq, own specimen")]
P("POST-HOC check 2.")
import os
for lab,rn,acc,asm,what in J:
    if not os.path.exists(W+lab+".bam.bai"): P("\n== ",lab,"not mapped yet"); continue
    bam=W+lab+".bam"; t=T[acc]["trny"]; g=int(t["g34_pos"]); n33=g-1; a,b=sorted((int(t["start"]),int(t["end"]))); L=len("".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna") if l[0]!=">"))
    P(f"\n== {lab} {rn} ({what}); {acc} {L} bp; trnY {a}-{b} (+), N33 at {n33}")
    dep=[int(x.split("\t")[2]) for x in run("depth","-aa","-G","2304","-r",acc,bam).splitlines()]; d20=[int(x.split("\t")[2]) for x in run("depth","-aa","-G","2304","-Q","20","-r",acc,bam).splitlines()]
    if not dep: P("  no alignments on the mitogenome"); continue
    P("  mitochondrial depth (primary): whole-genome median",st.median(dep),"range",min(dep),max(dep),"| MAPQ>=20 median",st.median(d20),"| at N33: all",dep[n33-1],"MAPQ>=20",d20[n33-1],"| trnY gene median",st.median(dep[a-1:b]),"| 1 kb each side median",st.median(dep[a-1001:a-1]),st.median(dep[b:b+1000]))
    G=sorted([(x[0],x[1],x[3],x[2]) for x in T[acc]["pcg"]]); near=sorted(G,key=lambda x:min(abs(x[0]-b),abs(a-x[1])))[:4]; P("  neighbouring genes (profile hits):",[(x[2],f"{x[0]}-{x[1]}",'median',st.median(dep[x[0]-1:x[1]]),'MAPQ>=20',st.median(d20[x[0]-1:x[1]])) for x in sorted(near)])
    cv=[]
    for l in run("coverage","-q","0","--ff","2304",bam).splitlines():
        if l[0]=="#": continue
        f=l.split("\t")
        if f[0]!=acc and int(f[2])>1000000: cv.append((f[0],int(f[2]),float(f[6])))
    P("  nuclear mean depth over contigs >1 Mb (length-weighted):",round(sum(x[1]*x[2] for x in cv)/max(1,sum(x[1] for x in cv)),2),"contigs",len(cv))
    R=[l.split("\t") for l in run("view","-F","2304",bam,f"{acc}:{n33}-{n33}").splitlines()]; by=collections.Counter(); full=0; clip=[]; nm=[]; ln=[]
    for f in R:
        bs=base_at(f,n33); q="MAPQ>=20" if int(f[4])>=20 else "MAPQ 1-19" if int(f[4]) else "MAPQ 0"; by[(q,(bs or "?").replace("T","U"))]+=1; c=cig(f[5]); al=sum(n for n,k in c if k in"M=XI"); sc=sum(n for n,k in c if k in"SH"); e=sum(n for n,k in c if k in"M=XD"); ln.append(al+sc); clip.append(sc); x=[z for z in f[11:] if z.startswith("NM:i:")]; nm.append(int(x[0][5:])/max(1,e) if x else 0); full+=sc<=50
    P("  reads with a primary alignment over N33 on the mitogenome:",len(R),{f"{k[0]} base {k[1]}":v for k,v in by.items()}); P("    read length median",st.median(ln) if ln else None,"| aligned end to end (<=50 nt clipped):",full,"| divergence from the record median",round(st.median(nm),5) if nm else None,"max",round(max(nm),4) if nm else None)
    if "hifi" in lab:
        ids={f[0] for f in R}; sup=collections.Counter()
        open(W+lab+".ids","w").write("\n".join(ids)+"\n")
        for l in run("view","-f","2048","-N",W+lab+".ids",bam).splitlines(): sup["mitogenome" if l.split("\t")[2]==acc else "nuclear"]+=1
        P("    supplementary parts of those reads:",dict(sup),"(a part on the mitogenome = the read crosses the record's origin)")
    # nuclear copies covering N33
    N=[]
    for l in open("results/gate3/numt_craniella_%s.txt"%("infr" if acc.startswith("OZ3") else "zetl")):
        m=re.match(r"\s+covers N33: (\S+) id (\S+) len (\d+)",l)
        if m: N.append((m[1],float(m[2]),int(m[3])))
    B="env/prefix/bin/"; s="".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna") if l[0]!=">").upper(); q=W+acc+".q"; open(q,"w").write(f">q\n{s}\n"); r=subprocess.run([B+"blastn","-task","blastn","-query",q,"-db",f"results/work/numt/{asm}","-evalue","1e-10","-dust","no","-max_target_seqs","2000","-num_threads","4","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout; sites={}
    for l in r.splitlines():
        f=l.split("\t")
        if f[0]==acc or int(f[2])<200: continue
        qs,ss,se=int(f[3]),int(f[5]),int(f[6]); step=1 if se>ss else -1; pos=ss; k=qs
        for x,y in zip(f[7],f[8]):
            if k==n33 and x!="-" and y!="-": sites[(f[0],pos)]=(float(f[1]),int(f[2]),step,y)
            if x!="-": k+=1
            if y!="-": pos+=step
    tot=collections.Counter(); rows=[]
    for (c,p),(idn,ln_,step,by_) in sites.items():
        z=collections.Counter(); hq=0
        for f in (l.split("\t") for l in run("view","-F","2304",bam,f"{c}:{p}-{p}").splitlines()):
            bb=base_at(f,p)
            if bb and bb!="-": bb=bb if step>0 else {"A":"T","C":"G","G":"C","T":"A"}.get(bb,bb); z[bb.replace("T","U")]+=1; hq+=int(f[4])>=20
        rows.append((c,p,idn,ln_,sum(z.values()),hq,dict(z))); tot.update(z)
    P("  nuclear copies whose sequence covers N33:",len(rows),"| reads placed on them (primary):",sum(x[4] for x in rows),"MAPQ>=20:",sum(x[5] for x in rows),"| bases",dict(tot),"| depth per copy: median",st.median([x[4] for x in rows]) if rows else None)
    for x in sorted(rows,key=lambda x:-x[2])[:5]: P("     ",x)
    tt=len(R)+sum(x[4] for x in rows); P(f"  => share of all N33-covering reads placed on the mitogenome: {len(R)}/{tt} = {len(R)/max(1,tt):.2f}; confidently (MAPQ>=20): {sum(v for k,v in by.items() if k[0]=='MAPQ>=20')}")
P("\nDONE")
