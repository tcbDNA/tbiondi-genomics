"""Map read hits to human NP_061914.3 columns via loop.aln; list reads covering human 196-222, with their translated loop segment."""
import sys,json,re
A={};cur=None
for l in open("data/mtrf/loop.aln"):
    if l[0]==">": cur=int(l[1:]); A[cur]=""
    else: A[cur]+=l.strip()
h=A[0]; hp={}  # alignment column -> human position
p=0
for i,ch in enumerate(h):
    if ch!="-": p+=1; hp[i]=p
def col2human(sid):  # subject residue (1-based) -> human pos
    m={};q=0
    for i,ch in enumerate(A[sid]):
        if ch!="-": q+=1; m[q]=hp.get(i)
    return m
GC={}
B="TCAG"; T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE={a+b+c:T[16*B.index(a)+4*B.index(b)+B.index(c)] for a in B for b in B for c in B}
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def trans(s): return "".join(CODE.get(s[i:i+3],"X") for i in range(0,len(s)-2,3))
for run in sys.argv[1:]:
    n=0; print("##",run)
    for l in open(f"data/mtrf/reads/{run}.hits.tsv"):
        f=l.rstrip("\n").split("\t"); sid,pid,ln,ev,qs,qe,ss,se,b,seq=int(f[1]),f[2],f[3],float(f[4]),int(f[5]),int(f[6]),int(f[7]),int(f[8]),f[9],f[10]
        m=col2human(sid); a=[m.get(x) for x in range(ss,se+1) if m.get(x)]; 
        if not a or min(a)>222 or max(a)<196: continue
        s=seq if qs<=qe else rc(seq); st=(min(qs,qe)-1) if qs<=qe else len(seq)-max(qs,qe)
        pep=trans(s[st%3 if False else st:][:]) ; n+=1
        print(f[0][-16:],"sid",sid,"id",pid,"human",min(a),"-",max(a),"e",ev,"aln_pep",trans(s[st:st+3*int(ln)+6])[:60])
    print("reads covering loop:",n)
