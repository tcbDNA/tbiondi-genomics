"""Plicatula muricata mtRF1a: rebuild CDS from genome + GFF, locate loop codons, compare to the deposited protein. (panel species)"""
import sys,re; sys.path.insert(0,"src"); import mtrf_guard; mtrf_guard.check(["FIG_Pmu"])
gid="PMU02G001039"; cds=[]; chrom=None
for l in open("data/figshare/Pmu.chr.gff3"):
    if gid in l:
        f=l.rstrip("\n").split("\t")
        if f[2]=="CDS": cds.append((int(f[3]),int(f[4]),f[6],f[7])); chrom=f[0]
        if f[2] in("gene","mRNA"): print(f[2],f[0],f[3],f[4],f[6],f[8][:90])
st=cds[0][2]; cds.sort(reverse=(st=="-")); print("exons",len(cds),cds)
lo,hi=min(a for a,*_ in cds)-2000,max(b for _,b,*_ in cds)+2000
seq=[];on=False;pos=0;buf=""
with open("data/figshare/plic/Pmu.chr.fna") as fh:
    for l in fh:
        if l[0]==">": on=(l[1:].split()[0]==chrom); pos=0; continue
        if on:
            n=len(l)-1; s,e=pos+1,pos+n
            if e>=lo and s<=hi: buf+=l.strip()[max(0,lo-s):hi-s+1]
            pos+=n
            if pos>hi: break
g=buf.upper(); rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
B="TCAG"; T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; CODE={a+b+c:T[16*B.index(a)+4*B.index(b)+B.index(c)] for a in B for b in B for c in B}
ex=[g[a-lo:b-lo+1] for a,b,*_ in cds]; ex=[rc(x) for x in ex] if st=="-" else ex
c="".join(ex); p="".join(CODE.get(c[i:i+3],"X") for i in range(0,len(c)-2,3)); print("CDS",len(c),"mod3",len(c)%3,"prot",len(p),"internal stops",p[:-1].count("*"))
D="";f=0
for l in open("data/figshare/Pmu.chr.pep"):
    if l[0]==">": f=gid in l
    elif f: D+=l.strip()
print("genomic translation == deposited protein:",p.rstrip("*")==D.rstrip("*"))
i=p.find("GGVHRVQRKPVTD"); L=p[i-1:i+26]; print("loop",L,"at aa",i+1)
for nm,off in (("Q203 eq",7),("T216 eq",20)):
    k=i-1+off; print(nm,p[k],"codon",c[3*k:3*k+3],"context",p[k-3:k+4])
b=0
for a,x in zip(cds,ex):
    if b<=3*(i-1)<b+len(x) or b<=3*(i+25)<b+len(x): print("loop lies in exon",a[:2],"CDS nt",b,b+len(x))
    b+=len(x)
