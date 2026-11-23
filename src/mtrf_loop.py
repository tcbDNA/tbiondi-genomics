"""Align RF1-type candidates (best human hit MTRF1L/MTRF1, genome models minus bacterial contaminants) and genotype human NP_061914.3 positions 203 and 216,
plus codon-recognition-loop indels and divergence. Anchors are columns of human MTRF1L motif residues: HRVQRVPKTE (200-209) and GRVHTST (212-218)."""
import csv,subprocess,sys,re,json,collections
B="env/prefix/bin/"
L={r["protein"]:r for r in csv.DictReader(open("results/work/mtrf_labels.tsv"),delimiter="\t")}
S={};cur=None
for l in open("data/mtrf/candidates.faa"):
    if l[0]==">": cur=l[1:].strip(); S[cur]=""
    else: S[cur]+=l.strip()
H={};cur=None
for l in open("data/mtrf/human_refs.faa"):
    if l[0]==">": cur=l[1:].split()[0]; H[cur]=""
    else: H[cur]+=l.strip()
hs=H["HUMAN_MTRF1L"]; assert hs[202]=="Q" and hs[215]=="T", (hs[200:210],hs[210:220])
IDX={p:p-1 for p in range(1,381)}
CONT={"EniAen_m1","EniAen_m2","EniAen_m3","RadSim_m1"}
recs=[("HUMAN_MTRF1L",hs),("HUMAN_MTRF1",H["HUMAN_MTRF1"])]; meta={}
for k,s in S.items():
    r=L.get(k)
    if r and r["best_human"] in("HUMAN_MTRF1L","HUMAN_MTRF1") and 100<=len(s)<=700: recs.append((k,s)); meta[k]=dict(src=k.split("|")[0],cls=r["best_human"])
for l in open("data/mtrf/genome_models.faa"):
    pass
cur=None;G={}
for l in open("data/mtrf/genome_models.faa"):
    if l[0]==">": cur=l[1:].strip(); G[cur]=""
    else: G[cur]+=l.strip()
for k,s in G.items():
    m=k.split("|")[1]
    if m not in CONT: recs.append((k,s)); meta[k]=dict(src="GENOME_"+m.split("_")[0],cls="genome_model")
open("data/mtrf/loop_in.faa","w").write("".join(f">{i}\n{s}\n" for i,(k,s) in enumerate(recs)))
r=subprocess.run([B+"mafft","--quiet","--localpair","--maxiterate","200","--thread","8","data/mtrf/loop_in.faa"],capture_output=True,text=True); open("data/mtrf/loop.aln","w").write(r.stdout)
A={};cur=None
for l in r.stdout.splitlines():
    if l[0]==">": cur=int(l[1:]); A[cur]=""
    else: A[cur]+=l.strip()
h=A[0]; colmap={};p=0
for i,ch in enumerate(h):
    if ch!="-": p+=1; colmap[p]=i
out=[]
for i,(k,s) in enumerate(recs):
    a=A[i]; res={203:a[colmap[203]],216:a[colmap[216]]}
    loop=(colmap[196],colmap[222]); seg=a[loop[0]:loop[1]+1]; hseg=h[loop[0]:loop[1]+1]
    ins=sum(1 for x,y in zip(seg,hseg) if y=="-" and x!="-"); dele=sum(1 for x,y in zip(seg,hseg) if x=="-" and y!="-")
    al=[(x,y) for x,y in zip(seg,hseg) if x!="-" and y!="-"]; ident=sum(x==y for x,y in al)/max(1,len(al))
    out.append(dict(id=k,source=meta.get(k,{}).get("src","HUMAN"),cls=meta.get(k,{}).get("cls",k),p203=res[203],p216=res[216],loop=a[loop[0]:loop[1]+1].replace("-",""),loop_len_vs_human=len(seg.replace("-",""))-len(hseg.replace("-","")),ins=ins,dele=dele,ident=round(ident,3),covers=res[203]!="-" and res[216]!="-"))
json.dump(out,open("results/work/mtrf_loop.json","w"),indent=0)
w=csv.writer(open("results/work/mtrf_loop.tsv","w"),delimiter="\t"); w.writerow(list(out[0])); [w.writerow(list(o.values())) for o in out]
print(len(out),collections.Counter(o["p216"] for o in out if o["covers"]))
for o in out:
    if o["source"].startswith(("FIG","GENOME","HUMAN")) or o["p216"]!="T": print(o["source"],o["id"][-28:],"203",o["p203"],"216",o["p216"],"ident",o["ident"],"ins",o["ins"],"del",o["dele"],o["loop"])
