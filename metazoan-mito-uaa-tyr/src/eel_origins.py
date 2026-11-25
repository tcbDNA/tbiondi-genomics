"""Eels (Anguilliformes + other Elopomorpha): independent-origin estimate for C33. (a) NCBI taxonomy parsimony; (b) FastTree on HMM-column codon alignment of 12 OXPHOS genes (all Elopomorpha with resolved N33 in the sweep), Hartigan changes and a permutation-free count on the tree."""
import sys,os,csv,json,subprocess,collections; sys.path.insert(0,"src"); import uaa_caller as U, origins
R=[r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included" and r["order"] in("Albuliformes","Anguilliformes","Elopiformes","Notacanthiformes") and r["n33"]=="resolved"]
G=["COX1","CYTB","ND4","ND5","ND1","ND2","COX2","COX3","ATP6"]; MS={str(h.name):h.M for h in U.hmms()}
def path(r): return "data/sweep/fa/"+str(sum(map(ord,r["acc"]))%200).zfill(3)+f"/{r['acc']}.fna"
def row(r):
    n,L0,kept,recs,*_=U.analyse(path(r),r["topology"]!="linear"); d=collections.defaultdict(dict)
    for x in recs:
        if x[0] in G: d[x[0]].setdefault(x[2],x[1])
    return r["acc"],"".join("".join((d[g].get(c,"---") if set(d[g].get(c,"---"))<=set("ACGT-") else "---") for c in range(1,MS[g]+1)) for g in G)
from concurrent.futures import ProcessPoolExecutor
def main():
    print("resolved Elopomorpha",len(R),collections.Counter((r["order"],r["base"]) for r in R)); print("C33:",[(r["acc"],r["species"],r["family"]) for r in R if r["base"]=="C"])
    par=origins.load_nodes(); print("taxonomy origins of C33 (Anguilliformes only):",origins.count({r["taxid"]:int(r["base"]=="C") for r in R if r["order"]=="Anguilliformes"},par),"| all Elopomorpha:",origins.count({r["taxid"]:int(r["base"]=="C") for r in R},par))
    fam=collections.defaultdict(lambda:[0,0]); [fam.__getitem__(r["family"]).__setitem__(0,fam[r["family"]][0]+1) or (fam[r["family"]].__setitem__(1,fam[r["family"]][1]+(r["base"]=="C"))) for r in R if r["order"]=="Anguilliformes"]; print("families with C33 (C/all):",{k:f"{v[1]}/{v[0]}" for k,v in fam.items() if v[1]})
    global t,B
    with ProcessPoolExecutor(3) as ex: A=dict(ex.map(row,R))
    open("results/work/eel.fa","w").write("".join(f">{a}\n{s}\n" for a,s in A.items() if s.count("-")<0.6*len(s))); print("alignment",len(A),"len",len(next(iter(A.values()))))
    subprocess.run("env/prefix/bin/FastTree -quiet -gtr -nt results/work/eel.fa > results/gate3/eel_tree.nwk 2>/dev/null",shell=True)
    from Bio import Phylo; t=Phylo.read("results/gate3/eel_tree.nwk","newick"); B={r["acc"]:r["base"] for r in R}; tips=[x.name for x in t.get_terminals()]
    def root_at(n):
        t.root_with_outgroup(n)
    og=[r["acc"] for r in R if r["order"] in("Elopiformes",) and r["acc"] in tips]
    if og: root_at(og[0])
    def fitch(c):
        if c.is_terminal(): return {B[c.name]=="C"},0
        s=[fitch(x) for x in c.clades]; n=sum(x[1] for x in s); cnt=collections.Counter(v for a,_ in s for v in a); m=max(cnt.values()); return {k for k,v in cnt.items() if v==m},n+len(s)-m
    print("tree: minimum C<->U changes",fitch(t.root)[1],"tips",len(tips))
    def nodes(c):
        for x in c.find_clades(order="preorder"): yield x
    C33=[x for x in tips if B[x]=="C"]; mono=0; grp=[]
    for x in t.find_clades():
        n=[y.name for y in x.get_terminals()]
        if all(B[k]=="C" for k in n) and len(n)>=2: grp.append(sorted(n))
    print("all-C33 clades (>=2 tips):",[(len(g),[next(r['species'] for r in R if r['acc']==a) for a in g][:4]) for g in sorted(grp,key=len,reverse=True)[:8]]); print("C33 tips",len(C33),"support values not evaluated")

if __name__=='__main__': main()
