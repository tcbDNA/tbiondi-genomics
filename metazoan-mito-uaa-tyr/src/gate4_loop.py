"""Gate 4: align MTRF1L-like candidates and genotype human NP_061914.3 positions 203 and 216 and the loop 196-222, as mtrf_loop.py. Adds a bacterial RF1 check (best blastp identity to E. coli RF1/RF2 vs human)."""
import csv,subprocess,sys,os,collections,json
sys.path.insert(0,"src"); import mtrf_guard
B="env/prefix/bin/"; G="data/mtrf/g4/"; W="results/work/g4/"
def fa(p):
    S={};c=None
    for l in open(p):
        if l[0]==">": c=l[1:].strip(); S[c]=""
        else: S[c]+=l.strip()
    return S
C=fa(G+"cand.faa"); L={r["protein"]:r for r in csv.DictReader(open("results/gate4/labels.tsv"),delimiter="\t")}; mtrf_guard.check([k.split("|")[0] for k in C])
H=fa("data/mtrf/human_refs.faa"); H={k.split()[0]:v for k,v in H.items()}; hs=H["HUMAN_MTRF1L"]; assert hs[202]=="Q" and hs[215]=="T"
K=[k for k in C if L.get(k) and L[k]["best_human"]=="HUMAN_MTRF1L" and 100<=len(C[k])<=700]
K2=[k for k in C if L.get(k) and L[k]["best_human"]=="HUMAN_MTRF1" and 100<=len(C[k])<=700]
recs=[("HUMAN_MTRF1L",hs),("HUMAN_MTRF1",H["HUMAN_MTRF1"])]+[(k,C[k]) for k in K+K2]
open(W+"loop_in.faa","w").write("".join(f">{i}\n{s}\n" for i,(k,s) in enumerate(recs)))
r=subprocess.run([B+"mafft","--quiet","--localpair","--maxiterate","200","--thread","6",W+"loop_in.faa"],capture_output=True,text=True); A={};c=None
for l in r.stdout.splitlines():
    if l[0]==">": c=int(l[1:]); A[c]=""
    else: A[c]+=l.strip()
open("results/gate4/loop.aln","w").write(r.stdout)
h=A[0]; col={};p=0
for i,ch in enumerate(h):
    if ch!="-": p+=1; col[p]=i
# bacterial check
open(W+"cand_q.faa","w").write("".join(f">q{i}\n{C[k]}\n" for i,k in enumerate(K+K2)))
subprocess.run([B+"makeblastdb","-in","data/mtrf/bact_refs.faa","-dbtype","prot","-out",W+"db_bact"],capture_output=True)
b=subprocess.run([B+"blastp","-query",W+"cand_q.faa","-db",W+"db_bact","-evalue","1e-3","-outfmt","6 qseqid sseqid pident length bitscore"],capture_output=True,text=True).stdout; bac={}
for l in b.splitlines():
    q,s,pi,ln,bs=l.split("\t"); bac[q]=max(bac.get(q,(0,))[0:1]+(0,)) if False else max(bac.get(q,(0,"")),(float(bs),f"{s}:{pi}%"))
out=[]
for i,(k,s) in enumerate(recs):
    a=A[i]; seg=a[col[196]:col[222]+1]; hseg=h[col[196]:col[222]+1]; al=[(x,y) for x,y in zip(seg,hseg) if x!="-" and y!="-"]
    loopu=seg.replace("-",""); covered=a[col[203]]!="-" and a[col[216]]!="-"
    q=i-2; bb=bac.get(f"q{q}",(0,"")) if i>=2 else (0,"")
    out.append(dict(id=k,source=k.split("|")[0],label=("HUMAN" if i<2 else L[k]["best_human"].replace("HUMAN_","")),length=len(s),p203=a[col[203]],p216=a[col[216]],covered=covered,loop=loopu,loop_len_vs_human=len(loopu)-len(hseg.replace("-","")),ident_to_human=round(100*sum(x==y for x,y in al)/max(1,len(al))),bact_best=bb[1],bact_bits=bb[0]))
w=csv.DictWriter(open("results/gate4/loop.tsv","w"),fieldnames=list(out[0]),delimiter="\t"); w.writeheader(); w.writerows(out)
for o in out:
    if o["label"] in("MTRF1L","HUMAN"): print(f"{o['source']:26s} {o['id'][-32:]:34s} {o['label']:7s} len {o['length']:4d} 203 {o['p203']} 216 {o['p216']} cov {o['covered']} loopΔ {o['loop_len_vs_human']:+d} id {o['ident_to_human']}% bact {o['bact_best']} | {o['loop']}")
