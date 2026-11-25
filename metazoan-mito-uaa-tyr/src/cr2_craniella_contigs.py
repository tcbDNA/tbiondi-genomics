"""Cold review 2, check 4 supplement for Craniella (POST-HOC): the HiFi libraries cannot be assembled cheaply, so the deposited assembly is used. Depth of the own-specimen HiFi run (check 2 alignments) in 5-kb windows of every contig; windows at >=50% of the mitochondrial depth (record written twice, craniella_depth3.txt) outside the mitogenome record are merged and searched with the locked extractor and ARAGORN.
Limit: a molecule that the assembly does not contain is not tested. usage: cr2_craniella_contigs.py -> appends to results/posthoc2/secondcopy/highcopy_craniella.txt"""
import sys,os,subprocess,collections,re
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import trny_extract as TE; M="env/tools/mapenv/bin/samtools"; B="env/prefix/bin/"; W="results/work/r2/cran/"; o=open("results/posthoc2/secondcopy/highcopy_craniella.txt","w"); P=lambda *x:(print(*x,flush=True),print(*x,file=o,flush=True)); P("POST-HOC check 4 supplement, Craniella (deposited assembly, own-specimen HiFi depth)")
for lab,acc,asm,md in (("infr_hifi","OZ318138.1","GCA_974939545.1",46),("zetl_hifi","OZ407450.1","GCA_978397455.1",49)):
    S={};c=None
    for l in open(f"results/work/numt/{asm}.fna"):
        if l[0]==">": c=l[1:].split()[0]; S[c]=[]
        else: S[c].append(l.strip())
    S={k:"".join(v).upper() for k,v in S.items()}; win=collections.defaultdict(list)
    p=subprocess.Popen([M,"depth","-a","-G","2304",W+lab+".bam"],stdout=subprocess.PIPE,text=True); cur=None; acc_=[0,0]
    for l in p.stdout:
        k,x,d=l.split("\t"); i=(int(x)-1)//5000
        if cur!=(k,i):
            if cur and acc_[1]>=500 and acc_[0]/acc_[1]>=0.5*md: win[cur[0]].append((cur[1],round(acc_[0]/acc_[1],1)))
            cur=(k,i); acc_=[0,0]
        acc_[0]+=int(d); acc_[1]+=1
    if cur and acc_[1]>=500 and acc_[0]/acc_[1]>=0.5*md: win[cur[0]].append((cur[1],round(acc_[0]/acc_[1],1)))
    n=sum(len(v) for k,v in win.items() if k!=acc); P(f"\n{acc} ({asm}, {len(S)} contigs): mitochondrial depth {md}x; 5-kb windows at >=50% of it outside the record: {n} on {len([k for k in win if k!=acc])} contigs")
    for k,v in win.items():
        if k==acc: P("   the record itself:",len(v),"windows"); continue
        i=0
        while i<len(v):
            j=i
            while j+1<len(v) and v[j+1][0]==v[j][0]+1: j+=1
            a,b=v[i][0]*5000,min(len(S[k]),(v[j][0]+1)*5000); t=W+"hc.fa"; open(t,"w").write(f">w\n{S[k][a:b]}\n"); c=TE.call(t); r=subprocess.run([B+"aragorn","-mt","-gc5","-w",t],capture_output=True,text=True).stdout; ar=[" ".join(x.split()) for x in r.splitlines() if "(gta)" in x.lower()]
            bl=subprocess.run([B+"blastn","-query",t,"-subject",W+acc+".q","-evalue","1e-10","-dust","no","-outfmt","6 pident length"],capture_output=True,text=True).stdout; mt=sum(int(z.split()[1]) for z in bl.splitlines() if float(z.split()[0])>=90)
            P("   ",k,f"{a+1}-{b}","depth",[x[1] for x in v[i:j+1]][:4],"| bp similar to the mitogenome (>=90%):",mt,"| extractor:",(c.get("N33"),c.get("loop"),"score",c.get("score"),"stem",c.get("stem_pairs")) if c.get("status")=="called" else "none","| ARAGORN Tyr(GTA):",ar or "none"); i=j+1
