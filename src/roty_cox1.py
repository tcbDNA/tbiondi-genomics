"""R. reniformis: do independent COX1 fragments (nuccore partials, many labs) carry in-frame TAA at the positions where the CM003310 genome has TAA?"""
import sys,os,re,collections,subprocess; sys.path.insert(0,"src"); import ncbi, uaa_caller as U
os.environ["GATE3_UNLOCK"]="1"
fa="data/sweep/fa/"+str(sum(map(ord,"CM003310.1"))%200).zfill(3)+"/CM003310.1.fna"
n,L0,kept,recs,me,idx,H=U.analyse(fa,True); seq=U.load(fa)[1]; rc=U.rc
site={}
for r in recs:
    if r[0]=="COX1" and r[1]=="TAA": site[r[6]]=(r[2],r[3],r[5])   # + strand 0-based codon start (npos), col
nn,ids,_=ncbi.search("nuccore",'"Rotylenchulus reniformis"[Organism] AND (COX1 OR "cytochrome c oxidase subunit" OR COI) AND 300:800[SLEN]',100); t=ncbi.call("efetch",db="nuccore",id=",".join(ids),rettype="gb",retmode="text").decode(errors="replace")
recsgb=[b for b in t.split("\n//") if b.strip()]; print("records",len(recsgb))
fas=open("results/work/roty_cox.fa","w"); meta={}
for b in recsgb:
    m=re.search(r"^ACCESSION\s+(\S+)",b,re.M); o=re.search(r"^ORIGIN(.*)$",b,re.M|re.S); 
    s=re.sub(r"[^a-z]","",b.split("ORIGIN",1)[1].lower()) if "ORIGIN" in b else ""; tt=re.findall(r'/transl_table=(\d+)',b); tr=re.findall(r'/translation="([A-Z\s*]*)',b); au=re.search(r'AUTHORS\s+(.*)',b); ti=re.search(r'TITLE\s+(.*)',b); geo=re.search(r'geo_loc_name="([^"]*)',b) or re.search(r'country="([^"]*)',b)
    if not m or not s: continue
    meta[m[1]]=dict(tt=tt,ntrans=tr[0].replace("\n","").replace(" ","")[:60] if tr else "",author=(au[1].split(",")[0] if au else "")[:30],title=(ti[1] if ti else "")[:60],geo=geo[1] if geo else "",len=len(s)); fas.write(f">{m[1]}\n{s.upper()}\n")
fas.close(); subprocess.run(["env/prefix/bin/makeblastdb","-in",fa,"-dbtype","nucl","-out","results/work/db_rot"],capture_output=True)
r=subprocess.run(["env/prefix/bin/blastn","-query","results/work/roty_cox.fa","-db","results/work/db_rot","-task","blastn","-evalue","1e-20","-outfmt","6 qseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True)
cnt=collections.Counter(); ident=[]; hits={}
for l in r.stdout.splitlines():
    q,p,ln,qs,qe,ss,se,qq,sq=l.split("\t"); ss,se=int(ss),int(se)
    if q in hits: continue
    hits[q]=1; ident.append(float(p)); st=-1 if ss>se else 1
    # map each subject position (0-based +) to query base
    mp={}; a=ss; b=0
    for x,y in zip(qq,sq):
        if y!="-": 
            if x!="-": mp[a-1]=x
            a+=st
    for npos,(col,py,msk) in site.items():
        cod=[mp.get(npos+i) for i in range(3)]
        if None in cod: continue
        cod="".join(cod)[::-1].translate(str.maketrans("ACGTN","TGCAN")).translate(str.maketrans("ACGTN","TGCAN")) if st==-1 else "".join(cod)   # minus-strand hit: codon is the reversed aligned bases
        # site codon is read on '-' strand for COX1 (strand - of genome)
        cnt[(npos,col,"TAA" if cod=="TAA" else cod)]+=1
print("identity to CM003310",min(ident),max(ident),"n",len(ident)); from collections import defaultdict
agg=defaultdict(collections.Counter)
for (np_,col,c),k in cnt.items(): agg[(np_,col)][c]+=k
print("COX1 genome TAA codons covered by partial records (position, HMM col, p(Tyr)):")
for (np_,col),c in sorted(agg.items()): print("  ",np_,"col",col,"pY",round(site[np_][1],2),"masked",site[np_][2],dict(c))
print("records by table:",collections.Counter(m["tt"][0] if m["tt"] else "?" for m in meta.values())); print("countries:",collections.Counter(m["geo"].split(":")[0] for m in meta.values()).most_common(8)); print("authors:",collections.Counter(m["author"] for m in meta.values()).most_common(8))
