"""[P] Post-hoc (2026-11-25): poly(A) 3' ends of mitochondrial transcripts for the Tier 1 C33 fish that have RNA-seq (A. glossodonta SRR13036352 liver; L. verrillii SRR36303867 muscle, the counted fish).
Rules: results/posthoc2/sixfish_polyA_criteria.md. Reads: all spots of each run streamed from SRA (a first look used the first 20,000,000) and kept when they share a 25-mer with a gene 3' end (src/preprint_sixfish_pats.py, src/kreads.c); results/work/polyA/<run>.fa.
Reads are placed with minimap2 (sr) for the diagonal only; the junction is the leftmost split that maximises matches to the genome before it and A after it
 (a first version took the first acceptable split instead and placed some junctions 2 to 4 bases too far left, and a second let reads with fewer than 10 tail bases through; both corrected 2026-11-25 to follow the written rule).
usage: preprint_sixfish_polya.py   writes results/posthoc2/sixfish_polya_junctions.tsv (one row per gene) and sixfish_polya_positions.tsv (all junction positions with >=3 reads)"""
import sys,os,csv,re,subprocess,collections
sys.path.insert(0,"src"); import uaa_caller as U
M="env/tools/mapenv/bin/"; W="results/work/polyA/"; RUNS=[("Albula glossodonta","CM165726.1","SRR13036352"),("Lycenchelys verrillii","CM163399.1","SRR36303867")]
G=list(csv.DictReader(open("results/posthoc2/sixfish_terminal_genes.tsv"),delimiter="\t")); BIC={"ATP8":"ATP6","ND4L":"ND4"}
def path(a): return f"data/sweep/fa/{str(sum(map(ord,a))%200).zfill(3)}/{a}.fna"
out=[]; pos=[]
for sp,acc,run in RUNS:
    s=U.load(path(acc))[1]; L=len(s); ref=W+acc+".lin.fa"; open(ref,"w").write(f">mt\n{s+s[:2000]}\n")
    sam=subprocess.Popen(f"{M}minimap2 -t 6 -ax sr -k 15 -w 5 -s 30 --secondary=no {ref} {W}{run}.fa 2>{W}{run}.mm.log",shell=True,stdout=subprocess.PIPE,text=True)
    s2=s+s[:2000]; RC=U.rc(s2); n=0; J={"+":collections.defaultdict(set),"-":collections.defaultdict(set)}; TL={"+":collections.defaultdict(list),"-":collections.defaultdict(list)}; nm=[]; span=collections.Counter(); ST={"+":collections.defaultdict(set),"-":collections.defaultdict(set)}
    ends={}   # gene -> strand, end (strand coords exclusive), class
    for r in G:
        if r["acc"]==acc and r["end_1based"]:
            e=int(r["end_1based"]); ends.setdefault(r["gene"],(r["strand"],e if r["strand"]=="+" else L-e+1,r["stop_class"],r["flag"]))
    for l in sam.stdout:
        if l[0]=="@": continue
        f=l.split("\t"); fl=int(f[1])
        if fl&4 or fl&256 or fl&2048: continue
        seq=f[9]; p=int(f[3])-1; cig=[(int(a),b) for a,b in re.findall(r"(\d+)([MIDSH])",f[5])]; n+=1
        t=re.search(r"NM:i:(\d+)",l); nm.append((int(t[1]) if t else 0,sum(a for a,b in cig if b=="M")))
        # last M block diagonal (for + genes) and first M block diagonal (for - genes)
        q=0; g=p; bl=[]
        for a,b in cig:
            if b=="M": bl.append((q,g,a)); q+=a; g+=a
            elif b in "IS": q+=a
            elif b=="D": g+=a
        for st in "+-":
            if st=="+": rd,rf,(qq,gg,aa)=seq,s2,bl[-1]
            else:
                b0=bl[0]; rd=U.rc(seq); rf=RC; qq=len(seq)-(b0[0]+b0[2]); gg=len(s2)-(b0[1]+b0[2]); aa=b0[2]
            lo=qq; N=len(rd)
            # candidate splits inside or after the block
            best=(-1,None); m=0; pref=[0]; k0=None
            for i in range(lo,N):
                gi=gg+(i-qq); m+=(gi<len(rf) and rd[i]==rf[gi]); pref.append(m)
            for k in range(lo+30,N-9):                                  # first place where a tail of >=10 bases, >=90% A, could start
                if rd[k]=="A" and rd[k:k+10].count("A")>=9 and pref[k-lo]>=0.95*(k-lo): k0=k; break
            if k0 is not None:
                E=min(N,k0+25)                                          # same read interval scored for every split: genome matches before it, A after it; leftmost on ties
                for k in range(k0,min(k0+13,N)):
                    sc=pref[k-lo]+rd[k:E].count("A")
                    if sc>best[0]: best=(sc,(gg+(k-qq),len(re.match("A*",rd[k:])[0]),k))
                k=best[1][2]                                            # the rule (>=30 matched bases, >=10 tail bases, >=90% A) is applied at the best split
                if N-k<10 or rd[k:k+10].count("A")<9 or k-lo<30 or pref[k-lo]<0.95*(k-lo): best=(-1,None)
                if best[1] and rf[best[1][0]:best[1][0]+10].count("A")>4: best=(-1,None)
            if best[1]:
                gj,tl,_=best[1]; e=gj if st=="+" else gj-(len(s2)-L)      # end in strand coordinates of the L-long genome
                if e<0: e+=L
                e%=L; sp1=f[0].rsplit(".",1)[0]; J[st][e].add(sp1); TL[st][e].append(tl); ST[st][e].add(gg-qq)     # reads counted per spot (a mate pair counts once)
        # spans over a stop codon (bicistronic check): ungapped read covering 20 nt on both sides
        for gn,(st,e,cl,fg) in ends.items():
            a,b=(e-3,e) if st=="+" else (L-e,L-e+3)
            for qq,gg,aa in bl:
                for sh in (0,L):
                    if gg<=a+sh-20 and gg+aa>=b+sh+20:
                        cd=seq[qq+a+sh-gg:qq+b+sh-gg]; cd=cd if st=="+" else U.rc(cd); span[(gn,cd)]+=1
    ident=1-sum(a for a,_ in nm)/max(1,sum(b for _,b in nm)); print(sp,run,"reads placed",n,"| identity over aligned bases",round(100*ident,2),"| junction positions",sum(len(v) for v in J.values()),"| tail reads",sum(len(x) for v in J.values() for x in v.values()))
    for st in "+-":
        for e,x in sorted(J[st].items()):
            if len(x)>=3: pos.append(dict(acc=acc,run=run,strand=st,end_strand_coord=e,genome_1based=e if st=="+" else L-e+1,reads=len(x),median_A=sorted(TL[st][e])[len(x)//2] if x else "",max_A=max(TL[st][e])))
    for gn,(st,e,cl,fg) in ends.items():
        x=s if st=="+" else U.rc(s); k={"T--":1,"TA-":2}.get(cl,3); c0=e-k                        # start of the stop codon
        near={j:len(v) for j,v in J[st].items() if e-60<=j<=e+150}; tot=sum(near.values())
        # transcript codon at the stop for a junction j: genome up to j, then A
        def codon(j): return (x[c0:j]+"AAA")[:3] if j>c0 else "junction before the stop codon"
        at={j:v for j,v in near.items() if (abs(j-e)<=2 and j>c0) or (j>e and k==3 and j<=e+150)}; sup=collections.Counter()
        dist=set()
        for j,v in at.items():
            c=codon(j).replace("T","U"); sup[c]+=v
            if c=="UAA": dist|={(j,q) for q in ST[st][j]}
        uaa=sup.get("UAA",0); sp_=[(c,v) for (g2,c),v in span.items() if g2==gn]; sp_.sort(key=lambda z:-z[1]); ntot=sum(v for _,v in sp_)
        top=sorted(near.items(),key=lambda z:-z[1])[:3]
        if uaa>=5 and uaa>=0.5*max(1,tot): v=f"transcript-supported UAA ({'encoded' if k==3 and not fg else 'TA- plus poly(A)' if fg else cl+' plus poly(A)'})"
        elif gn in BIC and cl=="UAA" and ntot>=5 and sp_[0][0]=="TAA" and tot<5: v="UAA in bicistronic transcript (no 3' end expected here)"
        elif uaa>=5: v=f"transcript-supported UAA at {uaa} of {tot} tails ({'encoded' if k==3 else cl+' plus poly(A)'}); other tails nearby"
        elif sup and max(sup.values())>=5: c,vv=sup.most_common(1)[0]; v=f"transcript ends with {c}"
        elif tot>=5: v="poly(A) ends elsewhere near the gene end"
        else: v=f"too few tail reads ({tot})"
        out.append(dict(species=sp,acc=acc,run=run,gene=gn,strand=st,annotated_end_1based=e if st=="+" else L-e+1,annotated_class=cl,flag=fg,tail_reads_near_end=tot,tail_reads_UAA=uaa,UAA_tail_distinct_read_starts=len(dist),codon_counts=";".join(f"{a}:{b}" for a,b in sup.most_common()),top_junctions=";".join(f"{j-e:+d}:{v}" for j,v in top),reads_spanning_stop=ntot,spanning_codons=";".join(f"{c}:{v}" for c,v in sp_[:4]),verdict=v))
def w(p,R):
    with open(p,"w",newline="") as f: o=csv.DictWriter(f,list(R[0]),delimiter="\t",lineterminator="\n"); o.writeheader(); o.writerows(R)
w("results/posthoc2/sixfish_polya_junctions.tsv",out)
if pos: w("results/posthoc2/sixfish_polya_positions.tsv",pos)
for r in out: print(r["species"][:10],r["gene"],r["annotated_class"],r["tail_reads_near_end"],r["codon_counts"],r["top_junctions"],"| span",r["reads_spanning_stop"],r["spanning_codons"],"|",r["verdict"])
