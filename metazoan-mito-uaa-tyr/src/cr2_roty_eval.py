"""Cold review 2, check 6 (POST-HOC): evaluation against results/posthoc2/roty_criteria.md (written before mapping). usage: cr2_roty_eval.py recruit | pile | eval"""
import sys,os,re,json,subprocess,collections
W="results/work/r2/roty/"; M="env/tools/mapenv/bin/"; B="env/prefix/bin/"; RC=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(RC)[::-1]; O="results/posthoc2/roty_eval.txt"
mt="".join(l.strip() for l in open(W+"mt.fa") if l[0]!=">").upper(); AR=(22643,22704); CMI=(22632,22684); G34=22668; N33=G34+1   # minus-strand gene: N33 is the + strand base 1 above G34
def sam(args): return subprocess.run([M+"samtools"]+args,capture_output=True,text=True).stdout
def cig(c): return [(int(n),o) for n,o in re.findall(r"(\d+)([MIDNSHP=X])",c)]
def base_at(l,pos):
    f=l.split("\t"); p=int(f[3]); i=0
    for n,o in cig(f[5]):
        if o in"M=X":
            if p<=pos<p+n: return f[9][i+pos-p],ord(f[10][i+pos-p])-33 if f[10]!="*" else 40
            p+=n;i+=n
        elif o in"IS": i+=n
        elif o in"DN":
            if p<=pos<p+n: return "-",40
            p+=n
    return None,0
if sys.argv[1]=="pile":
    o=open(O,"w"); P=lambda *x:(print(*x),print(*x,file=o)); P("POST-HOC check 6. Criteria: results/posthoc2/roty_criteria.md (sha256 in its timestamp file), written before mapping.")
    L=[l for l in sam(["view","-F","2304",W+"all.bam","CM003310.1:%d-%d"%(N33,N33)]).splitlines()]; mq=collections.Counter(); c=collections.Counter(); q20=collections.Counter()
    for l in L:
        f=l.split("\t"); mq["MAPQ>=20" if int(f[4])>=20 else "MAPQ 1-19" if int(f[4])>0 else "MAPQ 0"]+=1; b,q=base_at(l,N33)
        if b: c[b]+=1; q20[b]+=(q>=20)
    cp={"A":"U","C":"G","G":"C","T":"A","-":"-"}; P("\nC1/C2. primary alignments covering N33 on CM003310.1:",len(L),dict(mq)); P("  tRNA base at N33 (all):",{cp[k]:v for k,v in c.items()},"(base quality >=20):",{cp[k]:v for k,v in q20.items()})
    # alternative alleles over the ARAGORN interval
    bad=[]; cov=[]
    for pos in range(AR[0],AR[1]+1):
        d=collections.Counter()
        for l in sam(["view","-F","2304","-q","1",W+"all.bam","CM003310.1:%d-%d"%(pos,pos)]).splitlines():
            b,q=base_at(l,pos)
            if b and q>=20: d[b]+=1
        t=sum(d.values()); cov.append(t); top=d.most_common(2)
        if len(top)>1 and top[1][1]/t>=.05: bad.append((pos,top))
        if top and top[0][0]!=mt[pos-1]: bad.append((pos,"major differs from record",top,mt[pos-1]))
    P("  depth (bq>=20) over the 62-nt interval: min",min(cov),"median",sorted(cov)[len(cov)//2],"max",max(cov),"| positions with a second base >=5% or a major base unlike the record:",bad or "none")
    dm=[int(x.split("\t")[2]) for x in sam(["depth","-aa","-G","2304","-r","CM003310.1",W+"all.bam"]).splitlines()]; P("  mitochondrial depth: median",sorted(dm)[len(dm)//2],"| trnY interval median",sorted(dm[AR[0]-1:AR[1]])[31],"| 500 nt each side median",sorted(dm[AR[0]-501:AR[0]-1]+dm[AR[1]:AR[1]+500])[500])
    # nuclear places that resemble the tRNA
    open(W+"trny.fa","w").write(">t\n"+rc(mt[AR[0]-1:AR[1]])+"\n"); subprocess.run([B+"makeblastdb","-in",W+"ref.fa","-dbtype","nucl","-out",W+"refdb"],capture_output=True)
    r=subprocess.run([B+"blastn","-task","blastn","-query",W+"trny.fa","-db",W+"refdb","-word_size","7","-evalue","10","-dust","no","-outfmt","6 sseqid pident length qstart qend sstart send evalue sseq"],capture_output=True,text=True).stdout.splitlines(); P("  places in the reference resembling the tRNA (blastn word 7, E<=10):",len(r))
    for l in r[:12]:
        f=l.split("\t"); u,v=sorted((int(f[5]),int(f[6]))); D=[int(x.split("\t")[2]) for x in sam(["depth","-aa","-r",f"{f[0]}:{u}-{v}",W+"all.bam"]).splitlines()]; P("    ",f[0],"id",f[1],"len",f[2],"q",f[3]+"-"+f[4],"s",f[5]+"-"+f[6],"E",f[7],"median depth",sorted(D)[len(D)//2] if D else 0,"criterion (>=80% over >=50 nt):",float(f[1])>=80 and int(f[2])>=50)
    # secondary / supplementary placements of the reads that cover N33
    ids={l.split("\t")[0] for l in L}; sec=collections.Counter(); 
    open(W+"n33.ids","w").write("\n".join(ids)+"\n")
    for l in sam(["view","-N",W+"n33.ids",W+"all.bam"]).splitlines():
        f=l.split("\t"); sec[("secondary" if int(f[1])&256 else "supplementary" if int(f[1])&2048 else "primary",f[2])]+=1
    P("  all placements of those read pairs:",sec.most_common(8)); o.close()
if sys.argv[1]=="recruit":
    # read names with a primary alignment in the window trnY +-3 kb, and separately anywhere on the mitochondrion; mates are taken with them
    for tag,reg in (("win","CM003310.1:%d-%d"%(AR[0]-3000,min(len(mt),AR[1]+3000))),("mt","CM003310.1")):
        ids={l.split("\t")[0] for l in sam(["view","-F","2304",W+"all.bam",reg]).splitlines()}; open(W+tag+".ids","w").write("\n".join(ids)+"\n"); print(tag,"pairs",len(ids),flush=True)
        for m in "12":
            keep=0
            with open(W+f"SRR26324223_{m}.fastq") as f,open(W+f"{tag}_{m}.fq","w") as o:
                while True:
                    h=f.readline()
                    if not h: break
                    s,p,q=f.readline(),f.readline(),f.readline()
                    if h[1:].split()[0] in ids: o.write(h+s+p+q); keep+=1
            print(tag,m,"reads",keep,flush=True)
        subprocess.run([B+"spades.py","--isolate","-1",W+tag+"_1.fq","-2",W+tag+"_2.fq","-o",W+"spades_"+tag,"-t","6","-m","12"],stdout=open(W+f"spades_{tag}.log","w"),stderr=subprocess.STDOUT); print(tag,"SPAdes contigs",sum(l[0]=='>' for l in open(W+f"spades_{tag}/contigs.fasta")) if os.path.exists(W+f"spades_{tag}/contigs.fasta") else "FAILED",flush=True)
if sys.argv[1]=="eval":
    sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import trny_extract as TE, trny_decide as TD, gate3_classify as GC
    o=open(O,"a"); P=lambda *x:(print(*x),print(*x,file=o)); PAIR={("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
    rna=lambda s:s.replace("T","U")
    def fa(p):
        S={};c=None
        for l in open(p):
            if l[0]==">": c=l[1:].strip(); S[c]=[]
            else: S[c].append(l.strip())
        return {k:"".join(v).upper() for k,v in S.items()}
    ref62=rc(mt[AR[0]-1:AR[1]])
    for tag in ("win","mt"):
        S=fa(W+f"spades_{tag}/contigs.fasta"); P(f"\nC3. local assembly '{tag}': contigs",len(S),"longest",max(len(v) for v in S.values()))
        hit=None
        for k,v in S.items():
            for st,x in (("+",v),("-",rc(v))):
                i=x.find(ref62)
                if i>=0: hit=(k,st,x,i)
        if not hit:
            subprocess.run([B+"makeblastdb","-in",W+f"spades_{tag}/contigs.fasta","-dbtype","nucl","-out",W+"sdb"],capture_output=True); r=subprocess.run([B+"blastn","-task","blastn","-query",W+"trny.fa","-db",W+"sdb","-word_size","7","-evalue","1e-5","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout; P("   the exact 62-nt sequence is not in the assembly; blastn:",r[:600]); continue
        k,st,x,i=hit; P("   contig",k,"holds the 62-nt ARAGORN interval exactly (tRNA sense",st+");",i,"nt before it and",len(x)-i-62,"after it on the contig")
        # flanks against the record (tRNA sense)
        full=rc(mt); j=full.find(ref62)
        for nm,u in (("5' flank",min(i,j,1500)),("3' flank",min(len(x)-i-62,len(full)-j-62,1500))):
            a=x[i-u:i] if nm[0]=="5" else x[i+62:i+62+u]; b=full[j-u:j] if nm[0]=="5" else full[j+62:j+62+u]; same=sum(p==q for p,q in zip(a,b)); P("   ",nm,u,"nt identical positions",same,f"({100*same/max(1,u):.2f}%, ungapped)",[(n if nm[0]=="3" else n-u,p,q) for n,(p,q) in enumerate(zip(a,b)) if p!=q][:6])
        subprocess.run([B+"makeblastdb","-in",W+"mt.fa","-dbtype","nucl","-out",W+"mtdb"],capture_output=True); open(W+"ctg.fa","w").write(">c\n"+x+"\n"); r=subprocess.run([B+"blastn","-query",W+"ctg.fa","-db",W+"mtdb","-dust","no","-outfmt","6 pident length mismatch gapopen qstart qend sstart send"],capture_output=True,text=True).stdout.splitlines(); P("   contig against CM003310.1 (top HSPs: identity, length, mismatches, gaps, contig, record):",[l.replace("\t"," ") for l in r[:3]])
        # C5: the locked extractor and callers on the contig (a 12-kb piece centred on the locus)
        a0=max(0,i-6000); piece=x[a0:i+6062]; open(W+f"locus_{tag}.fa","w").write(">locus\n"+piece+"\n"); c=TE.call(W+f"locus_{tag}.fa"); import gate3_run as GR
        trn,sec=GR.trnas(W+f"locus_{tag}.fa"); P("C5. locked extractor on the contig:",{q:c.get(q) for q in "status reason start end strand score anticodon N33 loop stem_pairs g34_pos".split()}); P("    second callers:",sec,"| ARAGORN tRNAs within 250 nt:",[t for t in trn["aragorn"] if abs(t[0]-(i-a0))<250]); z=GC.n33(dict(trny=c,second=sec,trna=trn,pcg=[]),None); P("    registered status (gate and audit functions unchanged; no annotation, tRNA calls used for overlap):",z["status"],z["reason"],"N33",z["N33"],"stem at model columns",z["stem"],"score",z["score"],"callers",z["n_agree"])
        # C4
        g=c.get("g34_pos") or (i-a0)+ref62.find("GTA",25)+1
        def struct(nm,s0,e0):
            t=rna(piece[s0:e0]); p34=(g-1)-s0; L=len(t); arm5=t[p34-7:p34-2]; arm3=t[p34+5:p34+10]; loop=t[p34-2:p34+5]; ac=sum((a,b) in PAIR for a,b in zip(arm5,arm3[::-1]))
            # acceptor: first 7 against the 7 before the discriminator (last base)
            acc=sum((a,b) in PAIR for a,b in zip(t[:7],t[L-8:L-1][::-1])); acc2=max((sum((a,b) in PAIR for a,b in zip(t[:7],t[L-7-d:L-d][::-1])),d) for d in (0,1,2))
            nb=[y for y in trn["aragorn"] if not(y[0]<=g<=y[1])]; ov=max([min(e0,y[1])-max(s0+1,y[0])+1 for y in nb]+[-999]); P(f"C4. {nm}: {t} ({L} nt); 5' of the arm {p34-7} nt, 3' of the arm {L-(p34+10)} nt; acceptor pairs {acc}/7 (best of 0-2 unpaired 3' bases: {acc2[0]}/7 with {acc2[1]}); anticodon arm {arm5}-[{loop}]-{arm3} {ac}/5, GUA {'yes' if loop[2:5]=='GUA' else 'no'}; largest overlap with a neighbouring ARAGORN tRNA {max(ov,0)} nt -> {'supported' if acc>=5 and ac>=4 and ov<=10 else 'not supported'}")
        s=i-a0; struct("ARAGORN interval",s,s+62)
        if c.get("start"): u,v=sorted((int(c["start"]),int(c["end"]))); struct("Infernal interval (locked call on the contig)",u-1,v)
    P("\nDONE"); o.close()
if sys.argv[1]=="note":
    import random
    o=open(O,"a"); P=lambda *x:(print(*x),print(*x,file=o)); PAIR={("A","T"),("T","A"),("G","C"),("C","G"),("G","T"),("T","G")}; full=rc(mt); rng=random.Random(1); P("\nNOTES WRITTEN AFTER THE RESULT (not part of the criteria)")
    for L in (53,62):
        n=20000; k5=k6=0
        for _ in range(n):
            i=rng.randrange(0,len(full)-L); t=full[i:i+L]; p=sum((a,b) in PAIR for a,b in zip(t[:7],t[L-8:L-1][::-1])); k5+=p>=5; k6+=p>=6
        P(f"  acceptor criterion base rate: random {L}-nt stretches of this mitogenome reach >=5/7 in {100*k5/n:.1f}% and >=6/7 in {100*k6/n:.1f}% (AT-rich; G.U allowed)")
    r=subprocess.run([B+"blastn","-task","blastn","-query",W+"trny.fa","-db",W+"refdb","-word_size","7","-evalue","1e-6","-dust","no","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout.splitlines()
    for l in r:
        f=l.split("\t")
        if f[0]=="CM003310.1": continue
        q=int(f[3]); b=None; k=q; p33=ref_p=rc(mt[AR[0]-1:AR[1]]).find("GTA",25)  # 0-based index of G34
        for x,y in zip(f[7],f[8]):
            if x!="-":
                if k-1==p33-1: b=y
                k+=1
        P("  nuclear copy",f[0],"id",f[1],"len",f[2],"base at N33:",(b or "-").replace("T","U"),"anticodon+loop region",f[8])
        d=[int(x.split("\t")[2]) for x in sam(["depth","-aa","-r",f"{f[0]}:{min(int(f[5]),int(f[6]))-1000}-{max(int(f[5]),int(f[6]))+1000}",W+"all.bam"]).splitlines()]; P("    local nuclear depth (median, +-1 kb)",sorted(d)[len(d)//2])
        # long mitochondrial segment?
    for l in open(W+"mt_vs_nuc.tsv"):
        f=l.split("\t")
        if f[0]=="LDKF01000008.1": P("  mitochondrial similarity on LDKF01000008.1:",l.strip())
    o.close()
