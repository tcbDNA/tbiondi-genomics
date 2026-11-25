"""R. reniformis mtRF1a codon 216 in reads: exact 31-mers centred on the codon (record = Ser codon; alternatives = the four Thr codons and other Ser codons), from independent projects. usage: gate4_s216_reads.py RUN..."""
import sys,os,collections; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import kmer_support as KS
B="TCAG"; T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; CODE={a+b+c:T[16*B.index(a)+4*B.index(b)+B.index(c)] for a in B for b in B for c in B}; rc=KS.rc
def fa(p):
    S={};c=None
    for l in open(p):
        if l[0]==">": c=l[1:].split()[0]; S[c]=[]
        else: S[c].append(l.strip())
    return {k:"".join(v) for k,v in S.items()}
import glob
P="GRVHSSTV"; hits=collections.Counter(); src=fa("results/work/g4/TSA_R_reniformis.fna")
for want in ("GGVP01154292.1","GGVQ01123836.1","GGVO01139257.1","GGVR01090710.1","GGVV01154875.1"):
    s=src.get(want,"").upper()
    for st,x in(("+",s),("-",rc(s))):
        for f in range(3):
            pr="".join(CODE.get(x[i:i+3],"X") for i in range(f,len(x)-2,3)); j=pr.find(P)
            if j>=0:
                codon_pos=f+3*(j+4); w=x[codon_pos-14:codon_pos+17]; print(want,st,"codon",x[codon_pos:codon_pos+3],"31-mer",w); hits[(w,x[codon_pos:codon_pos+3])]+=1
lab={}; (w,cod),n=hits.most_common(1)[0]
for c in [a+b+d for a in B for b in B for d in B]:
    if CODE[c] in("S","T"): lab[f"{CODE[c]}_{c}"]=[w[:14]+c+w[17:]]
for i,k in enumerate(("GGVP01154292.1",)):
    s=src[k]; m=len(s)//2
    for q,off in(("ctrlA",40),("ctrlB",len(s)-70)): lab[q]=[s[off:off+31]]
print("test 31-mer (record codon",cod+"):",w,"aa",CODE[cod])
for run in sys.argv[1:]:
    c,e=KS.support(run,lab,maxgb=25); print("run",run,{k:v for k,v in (c or {}).items()} if c is not None else e,flush=True)
