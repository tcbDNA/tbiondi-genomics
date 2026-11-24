"""Craniella audit: divergence, identical tRNA/flank, GenBank project metadata, trnY copies and tRNA complement, Tetillidae relatives in the sweep/NCBI."""
import sys,os,json,re,subprocess,csv,xml.etree.ElementTree as ET; sys.path.insert(0,"src"); import ncbi
R={}
for l in open("results/gate3/run/out.jsonl"):
    x=json.loads(l)
    if x["acc"] in("OZ318138.1","OZ407450.1","OZ465198.1"): R[x["acc"]]=x
load=lambda a:"".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,a))%200).zfill(3)+f"/{a}.fna") if l[0]!=">").upper(); rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
S={a:load(a) for a in R}
def trn(a,fl=0):
    d=R[a]["trny"]; x,y=sorted((d["start"],d["end"])); s=S[a][max(0,x-1-fl):y+fl]; return rc(s) if d["strand"]=="-" else s
A,B="OZ407450.1","OZ318138.1"
print("lengths",{a:len(S[a]) for a in S}); print("trnY",{a:(R[a]["trny"]["N33"],R[a]["trny"]["loop"],R[a]["trny"]["stem_pairs"],trn(a).replace('T','U')) for a in S}); print("trnY zetlandica == infrequens:",trn(A)==trn(B))
for fl in (100,300,1000): x,y=trn(A,fl),trn(B,fl); n=sum(p!=q for p,q in zip(x,y)); print(f"flank +/-{fl}: len {len(x)},{len(y)} equal-length mismatches {n if len(x)==len(y) else 'length differs'}")
for a in S: open(f"results/work/cran_{a}.fa","w").write(f">{a}\n{S[a]}\n"); subprocess.run(["env/prefix/bin/makeblastdb","-in",f"results/work/cran_{a}.fa","-dbtype","nucl","-out",f"results/work/cran_{a}"],capture_output=True)
for q,s in ((A,B),(A,"OZ465198.1"),(B,"OZ465198.1")):
    r=subprocess.run(["env/prefix/bin/blastn","-task","blastn","-query",f"results/work/cran_{q}.fa","-db",f"results/work/cran_{s}","-evalue","1e-20","-outfmt","6 pident length"],capture_output=True,text=True); L=[l.split("\t") for l in r.stdout.splitlines()]; al=sum(int(l[1]) for l in L); idt=sum(float(l[0])*int(l[1]) for l in L)/al if al else 0; print(q,"vs",s,"aligned",al,"identity",round(idt,2))
for a in (A,B): d=R[a]; print(a,"GUA loci ARAGORN",len(d["second"]["aragorn"]),"tRNAscan",len(d["second"]["trnascan"]),"verified GUA arms",d["trny"].get("n_ok"),"tRNAs: ARAGORN",len(d["trna"]["aragorn"]),"tRNAscan",len(d["trna"]["trnascan"]),"anticodons",sorted({x[4] for x in d["trna"]["aragorn"]}),"|dup GTA",[x for x in d["trna"]["aragorn"] if x[4]=="GTA"])
n,ids,_=ncbi.search("nuccore",'(Tetillidae[Organism] OR Cinachyrella[Organism] OR Tetilla[Organism]) AND mitochondrion[filter] AND 5000:60000[SLEN]',40); d=json.loads(ncbi.call("esummary",db="nuccore",id=",".join(ids),retmode="json"))["result"]; print("Tetillidae mitogenomes in nuccore:",[(d[i]["accessionversion"],d[i]["organism"],d[i]["slen"],d[i]["createdate"],d[i]["title"][:60]) for i in ids])
for a in (A,B):
    root=ET.fromstring(ncbi.call("efetch",db="nuccore",id=a,rettype="gb",retmode="xml"))
    for s in root.iter("GBSeq"): print(a,[x.findtext("GBXref_id") for x in s.iter("GBXref")],"project",s.findtext("GBSeq_project"),(s.findtext("GBSeq_comment") or "")[:300].replace("\n"," "))
