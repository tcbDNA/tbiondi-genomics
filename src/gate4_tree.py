"""Gate 4 tree: human MTRF1L/MTRF1/MRPL58/MTRFR, E. coli RF1/RF2, Gate 2 MTRF1L proteins and one Gate 4 model per source (best covered). FastTree support; list the nearest neighbours of each Gate 4 sequence."""
import csv,subprocess,sys,re
sys.path.insert(0,"src"); import mtrf_guard
B="env/prefix/bin/"; W="results/work/g4/"
def fa(p):
    S={};c=None
    for l in open(p):
        if l[0]==">": c=l[1:].strip(); S[c]=""
        else: S[c]+=l.strip()
    return S
C=fa("data/mtrf/g4/cand.faa"); rep={}
for r in csv.DictReader(open("results/gate4/loop.tsv"),delimiter="\t"):
    if r["label"]=="MTRF1L" and r["covered"]=="True": rep.setdefault(r["source"],[]).append((int(r["ident_to_human"]),r["id"]))
pick={s:max(v)[1] for s,v in rep.items()}; mtrf_guard.check(list(pick))
H=fa("data/mtrf/human_refs.faa"); E=fa("data/mtrf/bact_refs.faa"); G2=fa("data/mtrf/candidates.faa"); L={r["protein"]:r["best_human"] for r in csv.DictReader(open("results/work/mtrf_labels.tsv"),delimiter="\t")}
o=open(W+"tree_in.faa","w")
for k,s in H.items(): o.write(f">{k.split()[0]}\n{s}\n")
for k,s in E.items(): o.write(f">{k.split()[0]}\n{s}\n")
n={}
for k,s in G2.items():
    if L.get(k)=="HUMAN_MTRF1L" and 150<=len(s)<=600:
        src=k.split("|")[0]
        if n.get(src,0)>=1: continue
        n[src]=1; o.write(f">G2_{src}\n{s}\n")
for s,k in pick.items(): o.write(f">G4_{s}\n{C[k]}\n")
o.close()
subprocess.run(f"{B}mafft --quiet --auto --thread 4 {W}tree_in.faa > {W}tree.aln",shell=True); subprocess.run(f"{B}FastTree -quiet -lg {W}tree.aln > results/gate4/tree.nwk 2>/dev/null",shell=True)
from Bio import Phylo; t=Phylo.read("results/gate4/tree.nwk","newick"); t.root_at_midpoint() if not any("ECOLI" in x.name for x in t.get_terminals()) else t.root_with_outgroup([x for x in t.get_terminals() if "ECOLI" in x.name][0])
par={c:p for p in t.find_clades() for c in p.clades}
for x in t.get_terminals():
    if not x.name.startswith("G4_"): continue
    p=par.get(x); sib=[y.name for y in p.get_terminals() if y is not x][:4]; print(x.name,"| nearest:",sib,"| support",p.confidence)
