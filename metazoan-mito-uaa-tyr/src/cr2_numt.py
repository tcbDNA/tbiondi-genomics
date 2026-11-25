"""Cold review 2, check 8a / item 12 (POST-HOC): src/numt.py's exclusion rule against its stated rule. The script sums HSP lengths at >=99% identity and excludes a contig <200 kb when the sum exceeds 50% of the mitogenome; its header says >=99% over >=90% of the mitogenome.
Here: true coverage = merged mitogenome intervals of >=99% HSPs; contig excluded only at >=90%; and the length gate for hits that cover N33 is lowered from 200 to 50 bp (blastn E<=1e-3 around trnY). Reports what changes. writes results/posthoc2/numt_truecov.txt"""
import sys,os,json,csv,subprocess,collections
B="env/prefix/bin/"; W="results/work/numt"; o=open("results/posthoc2/numt_truecov.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True))
P("POST-HOC. NUMT checks exist for 6 assemblies of 5 of the 17 hits (two for A. glossodonta); the other hits have no nuclear assembly on disk.")
J=[("CM165726.1","GCA_019788955.1"),("CM165726.1","GCA_056783665.1"),("CM163399.1","GCA_056747025.1"),("OZ318138.1","GCA_974939545.1"),("OZ407450.1","GCA_978397455.1"),("KU726971.1","GCA_001723225.1")]
T={}
for l in open("results/gate3/run/out.jsonl"):
    x=json.loads(l)
    if x["acc"] in {a for a,_ in J}: T[x["acc"]]=x["trny"]
for acc,asm in J:
    d=T[acc]; g=int(d["g34_pos"]); site=g+1 if d["strand"]=="-" else g-1   # 1-based + strand position of N33
    q=f"{W}/{acc}.q.fa"; s="".join(l.strip() for l in open(q) if l[0]!=">"); L=len(s); rb=s[site-1]
    r=subprocess.run([B+"blastn","-task","blastn","-query",q,"-db",f"{W}/{asm}","-evalue","1e-10","-dust","no","-max_target_seqs","2000","-num_threads","4","-outfmt","6 sseqid pident length qstart qend sstart send evalue slen qseq sseq"],capture_output=True,text=True).stdout
    by=collections.defaultdict(list); [by[l.split("\t")[0]].append(l.split("\t")) for l in r.splitlines()]
    old=set();new=set();cov={}
    for c,h in by.items():
        x=[y for y in h if float(y[1])>=99]; su=sum(int(y[2]) for y in x); iv=sorted((int(y[3]),int(y[4])) for y in x); m=0;e=0
        for a,b in iv:
            if b>e: m+=b-max(a,e+1)+1; e=b
        slen=int(h[0][8]); cov[c]=(su,m/L,slen)
        if su>0.5*L and slen<200000: old.add(c)
        if m>=0.9*L and slen<200000: new.add(c)
    P(f"\n{acc} vs {asm} (mitogenome {L} bp, N33 + strand base {rb})"); P("  excluded by the script (sum >50%, contig <200 kb):",{c:(cov[c][0],f"true coverage {100*cov[c][1]:.1f}%","contig",cov[c][2]) for c in old}); P("  excluded under true coverage >=90%:",sorted(new)); P("  contigs that change class:",sorted(old^new))
    P("  contigs >=200 kb with >=99% HSPs covering >50% of the mitogenome:",[(c,f"{100*v[1]:.0f}%") for c,v in cov.items() if v[1]>0.5 and v[2]>=200000])
    def n33(H):
        out=[]
        for h in H:
            q0=int(h[3]); b=None; a=q0
            for x,y in zip(h[9],h[10]):
                if x!="-":
                    if a==site: b=y
                    a+=1
            if b is not None: out.append((h[0],float(h[1]),int(h[2]),b,int(h[5]),int(h[6])))
        return out
    hits=[h for c,H in by.items() if c not in new for h in H if int(h[2])>=200 and float(h[1])>=80]; x=n33(hits); P("  nuclear-like HSPs (>=200 bp, >=80%):",len(hits),"| covering N33:",len(x),"| base at N33:",dict(collections.Counter(y[3] for y in x)))
    xo=n33([h for c,H in by.items() if c not in old for h in H if int(h[2])>=200 and float(h[1])>=80]); P("  same with the script's exclusion:",len(xo),dict(collections.Counter(y[3] for y in xo)),"| changed:",len(x)!=len(xo))
    if old^new: P("    N33-covering HSPs on the contigs that change class:",n33([h for c in old^new for h in by[c]]))
    # short copies: query = trnY +-40 nt
    w0=max(1,site-110); w1=min(L,site+110); qs=f"{W}/{acc}.trny.fa"; open(qs,"w").write(f">t\n{s[w0-1:w1]}\n")
    r2=subprocess.run([B+"blastn","-task","blastn","-query",qs,"-db",f"{W}/{asm}","-evalue","1e-3","-dust","no","-max_target_seqs","5000","-word_size","7","-num_threads","4","-outfmt","6 sseqid pident length qstart qend sstart send evalue slen qseq sseq"],capture_output=True,text=True).stdout
    H=[l.split("\t") for l in r2.splitlines()]; off=w0-1
    sh=[]
    for h in H:
        if int(h[2])<50: continue
        a=int(h[3])+off; b=None
        for X,Y in zip(h[9],h[10]):
            if X!="-":
                if a==site: b=Y
                a+=1
        if b is not None and h[0] not in new and h[0]!=acc: sh.append((h[0],float(h[1]),int(h[2]),b))
    P("  short-copy search (>=50 bp, E<=1e-3, trnY +-110 nt) covering N33 outside the mitochondrial contig(s):",len(sh),"| base:",dict(collections.Counter(x[3] for x in sh)),"| other than the record base:",[x for x in sh if x[3]!=rb][:6])
P("\nDONE")
