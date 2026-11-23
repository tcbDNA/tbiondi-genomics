import csv
import sys as _s; _s.path.insert(0,"src"); import mtrf_guard
mtrf_guard.check([l[1:].split("|")[0].strip() for fn in ("data/mtrf/candidates.faa","data/mtrf/genome_models.faa") for l in open(fn) if l[0]==">"])
L={r["protein"]:r for r in csv.DictReader(open("results/work/mtrf_labels.tsv"),delimiter="\t")}
S={};cur=None
for l in open("data/mtrf/candidates.faa"):
    if l[0]==">": cur=l[1:].strip(); S[cur]=""
    else: S[cur]+=l.strip()
o=open("data/mtrf/tree_in.faa","w")
for fn in ("data/mtrf/human_refs.faa","data/mtrf/bact_refs.faa"):
    for l in open(fn): o.write(l if l[0]!=">" else ">"+l[1:].split()[0]+"\n")
n={}
for k,s in S.items():
    r=L.get(k)
    if not r: continue
    cls=r["best_human"].replace("HUMAN_",""); src=k.split("|")[0]; key=(src,cls)
    if cls in("MRPL58","MTRFR") and n.get(key,0)>=1: continue
    if cls in("MTRF1L","MTRF1") and n.get(key,0)>=3: continue
    n[key]=n.get(key,0)+1
    if len(s)<60 or len(s)>700: continue
    o.write(f">{cls}~{k.replace('|','_')}\n{s}\n")
for l in open("data/mtrf/genome_models.faa"): o.write(l if l[0]!=">" else ">GENOME~"+l[1:].strip().replace("|","_")+"\n")
o.close()
