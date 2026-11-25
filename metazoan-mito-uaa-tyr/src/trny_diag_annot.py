"""Diagnostic: for PANEL nematodes, list the annotated trnY sequence and what the anticodon arm looks like (panel allowlist enforced)."""
import csv,gzip,re,json,sys
sys.path.insert(0,"src"); import panel_guard
from Bio import SeqIO
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
want={k for k in P if P[k]["role"]!="pos" and P[k]["lineage"]!="Chordata" and k.startswith("NC_")}; panel_guard.check(want)
PAIR={("A","T"),("T","A"),("G","C"),("C","G"),("G","T"),("T","G")}; out={}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id not in want: continue
        L=[]
        for f in r.features:
            if f.type=="tRNA" and re.search(r"Tyr|trnY",str(f.qualifiers.get("product",""))+str(f.qualifiers.get("gene","")),re.I):
                s=str(f.extract(r.seq)).upper(); arms=[]
                for m in re.finditer("(?=GTA)",s):
                    i=m.start(); a,b=s[i-7:i-2],s[i+5:i+10]
                    if i>=7 and len(b)==5: arms.append((sum((x,y) in PAIR for x,y in zip(a,b[::-1])),i,s[i-1],s[i-2:i+5]))
                L.append(dict(start=int(f.location.start)+1,end=int(f.location.end),strand=f.location.strand,seq=s,arms=sorted(arms,reverse=True)))
        out[r.id]=L
json.dump(out,open("results/gate2/trny_diag/annotated_trny_panel.json","w"))
fail="NC_033868.1 NC_033869.1 NC_071371.1 NC_031516.1 NC_008640.1 NC_008828.1 NC_005928.1 NC_008046.1 NC_008047.1 NC_008692.1 NC_017750.1".split()
for k in fail:
    for x in out.get(k,[]): print(k,P[k]["species"],len(x["seq"]),x["start"],x["strand"],x["seq"],"best arm (pairs,pos,N33,loop):",x["arms"][:1])
