"""Cold review 2, check 8 / item 13 (POST-HOC): the two +14 P. vulnus mtRF1a models (results/gate4/loop.tsv). Are the 15 extra residues an isoform or an unspliced intron? Reads the miniprot GFF and the assembly; writes results/posthoc2/pvulnus_models.txt"""
import re,subprocess,os
o=open("results/posthoc2/pvulnus_models.txt","w"); P=lambda *x:(print(*x),print(*x,file=o)); RC=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(RC)[::-1]
def fa(p):
    S={};c=None
    for l in open(p):
        if l[0]==">": c=l[1:].split()[0]; S[c]=[]
        else: S[c].append(l.strip())
    return {k:"".join(v).upper() for k,v in S.items()}
B="TCAG"; T5="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; AA={a+b+c:T5[i] for i,(a,b,c) in enumerate((a,b,c) for a in B for b in B for c in B)}
tr=lambda s:"".join(AA.get(s[i:i+3],"X") for i in range(0,len(s)-2,3))
G=fa("data/mtrf/g4/Pratylenchus_vulnus.fna"); g=G["CM134649.1"]; lo,hi=11494700,11496400; reg=rc(g[lo:hi]); P("POST-HOC. region CM134649.1:%d-%d, minus strand"%(lo+1,hi))
# locate the pieces in the three frames of the coding strand
for nm,pep in (("inserted residues","VILIFKQHLFSPNKK"),("downstream exon","GGVHRVQRTPVTDKSRIHTSTASLA"),("read-through E+insert+exon","EVILIFKQHLFSPNKKGGVHRVQRTPVTD")):
    for f in range(3):
        i=tr(reg[f:]).find(pep)
        if i>=0: s=f+3*i; P(f"  {nm} {pep}: frame {f}, coding-strand offset {s}-{s+3*len(pep)} genomic {hi-s}..{hi-s-3*len(pep)+1}")
C={}
for l in open("data/mtrf/g4/cand.faa"):
    if l[0]==">": c=l[1:].strip(); C[c]=""
    else: C[c]+=l.strip()
for k,v in C.items():
    if k.startswith("GEN_P_vulnus") and "CM134649.1:114949" in k:
        j=v.find("GGVHR"); P("   model",k.split("|",1)[1],"len",len(v),"... ",v[max(0,j-45):j],"|",v[j:j+26])
for pep in ("EVILIFK","VILIFK","LIFKQHL","QHLFSPNKK","FSPNKKGGVHR","SPNKK","NKKGGVHRVQ","KGGVHRVQ"):
    for f in range(3):
        for m in re.finditer(pep,tr(reg[f:])): s_=f+3*m.start(); P(f"    piece {pep}: frame {f} offset {s_}-{s_+3*len(pep)} (coding strand of the region)")
P("  coding strand 480-760 of the region, translated in three frames:")
for st in range(480,740,60):
    P("   ",st,reg[st:st+60].lower()); 
    for f in range(3): q=st+((f-st)%3); P("        f%d "%f+" "*(q-st)+"  ".join(tr(reg[q:st+62])))
P("\n  miniprot models overlapping the locus (target, score, CIGAR-like cs):")
cur=None
for l in open("results/work/g4/GEN_P_vulnus.gff"):
    if l.startswith("##PAF"):
        f_=l.rstrip("\n").split("\t")
        if f_[6]=="CM134649.1" and 11490000<int(f_[8])<11500000:
            t={x[:2]:x[5:] for x in f_[13:]}; P("   ",f_[1][-22:],"q",f_[3],"-",f_[4],"of",f_[2],"t",f_[8],"-",f_[9],f_[5],"AS",t.get("AS"),"fs",t.get("fs"),"st",t.get("st"),"cg",t.get("cg"))
# the homologous intron in other tylenchids
P("\n  same junction in other genomes (exon structure from miniprot cg of the best model):")
for src,fn in (("GEN_R_reniformis","Rotylenchulus_reniformis"),("GEN_R_robustus","Rotylenchus_robustus")):
    best=None
    for l in open(f"results/work/g4/{src}.gff"):
        if l.startswith("##PAF"):
            f_=l.rstrip("\n").split("\t"); t={x[:2]:x[5:] for x in f_[13:]}
            if f_[1].startswith("HUMAN_MTRF1L") or "MTRF1L" in f_[1]: P("   ",src,f_[1][-22:],f_[6],f_[8],f_[9],f_[5],"AS",t.get("AS"),"cg",t.get("cg"))
S=fa("data/mtrf/g4/Rotylenchulus_reniformis.fna"); RR=S["LDKF01008021.1"]
for f in range(3):
    for strand,x in (("+",RR),("-",rc(RR))):
        t=tr(x[f:]); j=t.find("GVHRVQRLPFTEK")
        if j>=0: s=f+3*j; P("    R. reniformis LDKF01008021.1",strand,"around the junction:",x[s-70:s].lower(),"|",x[s:s+12],"upstream translation",t[j-25:j])
# protein evidence: does any annotated or transcript protein of Pratylenchus carry the insert?
P("\n  SUMMARY")
it=reg[613:660]; P("  intron, region 613-660 of the coding strand:",it,"length",len(it),"mod 3 =",len(it)%3,"| donor",it[:6],"acceptor",it[-6:],"| upstream exon ends",reg[604:613],tr(reg[604:613]),"| downstream exon starts",reg[660:672],tr(reg[660:672]))
import collections
n=collections.Counter()
for l in open("results/work/g4/GEN_P_vulnus.gff"):
    if l.startswith("##PAF"):
        f_=l.rstrip("\n").split("\t")
        if f_[6]=="CM134649.1" and 11490000<int(f_[8])<11500000:
            t={x[:2]:x[5:] for x in f_[13:]}; cg=t.get("cg",""); n["models"]+=1; n["47-nt intron at the junction"]+=bool(re.search(r"47[NUV]",cg)); n["2-nt frameshift + 15 codons (2F15D) instead"]+="2F15D" in cg; n["2F15D with frameshift flag fs>=1"]+=("2F15D" in cg and t.get("fs")!="0")
P("  miniprot models at the locus:",dict(n)); P("  2 + 15x3 = 47: the +14 models translate the 47-nt GT..AG intron after skipping 2 nt. Retaining a 47-nt intron would shift the frame (47 mod 3 = 2), so it cannot give the in-frame 'isoform'.")
