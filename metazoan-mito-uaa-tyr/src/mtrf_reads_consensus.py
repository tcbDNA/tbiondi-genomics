"""Strict loop consensus from reads: E<=1e-8, aligned >=25 aa, read spans human 200-222 (loop), translated in the aligned frame.
Reports read count and residue counts at human 203 and 216 plus distinct loop peptides."""
import sys,re,collections,json
import sys as _s; _s.path.insert(0,"src"); import mtrf_guard
mtrf_guard.check(sys.argv[1:])
A={};cur=None
for l in open("data/mtrf/loop.aln"):
    if l[0]==">": cur=int(l[1:]); A[cur]=""
    else: A[cur]+=l.strip()
h=A[0]; hp={};p=0
for i,ch in enumerate(h):
    if ch!="-": p+=1; hp[i]=p
def smap(sid):
    m={};q=0
    for i,ch in enumerate(A[sid]):
        if ch!="-": q+=1; m[q]=hp.get(i)
    return m
B="TCAG"; T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE={a+b+c:T[16*B.index(a)+4*B.index(b)+B.index(c)] for a in B for b in B for c in B}
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def run(r):
    c203=collections.Counter();c216=collections.Counter();pep=collections.Counter();n=0;ids=[]
    for l in open(f"data/mtrf/reads/{r}.hits.tsv"):
        f=l.rstrip("\n").split("\t"); sid,ln,ev,qs,qe,ss,se,seq=int(f[1]),int(f[3]),float(f[4]),int(f[5]),int(f[6]),int(f[7]),int(f[8]),f[10]
        if ev>1e-8 or ln<25: continue
        m=smap(sid); 
        # nt codons: query coordinate of each aligned residue (ungapped HSP assumed; diamond reports gapped, so use conserved anchor search in the translated frame)
        s=seq if qs<=qe else rc(seq); fr=(min(qs,qe)-1)%3 if qs<=qe else (len(seq)-max(qs,qe))%3
        pep_full="".join(CODE.get(s[i:i+3],"X") for i in range(fr,len(s)-2,3))
        a=[x for x in (m.get(y) for y in range(ss,se+1)) if x]
        if not a or min(a)>205 or max(a)<220: continue
        n+=1; ids.append(f[0]); pep[pep_full[:80]]+=1
    return n,pep
for r in sys.argv[1:]:
    n,pep=run(r); print("##",r,"strict loop-spanning reads:",n)
    for x,c in pep.most_common(6):
        print("   ",c,x)
