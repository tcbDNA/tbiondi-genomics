"""Annotated CDS and tRNA intervals for panel records (annotation only)."""
import csv,gzip,json,re,sys,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"src"); import panel_guard, ncbi
from Bio import SeqIO; from io import StringIO
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; panel_guard.check(P); out={}
def rec(pid,r):
    d={"cds":[],"trna":[]}
    for f in r.features:
        PT=[(int(x.start)+1,int(x.end),"+" if x.strand==1 else "-") for x in f.location.parts]   # post-lock fix: parts, not min..max
        if f.type=="CDS": d["cds"]+=[[u,v,w,(f.qualifiers.get("gene") or ["?"])[0]] for u,v,w in PT]
        elif f.type=="tRNA":
            t=str(f.qualifiers.get("product",""))+str(f.qualifiers.get("gene","")); d["trna"]+=[[u,v,w,(f.qualifiers.get("product") or ["?"])[0],bool(re.search(r"Tyr|trnY",t,re.I))] for u,v,w in PT]
    out[pid]=d if (d["cds"] or d["trna"]) else None
want={k for k in P if k.startswith("NC_")}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id in want: rec(r.id,r)
gb=[k for k in P if re.match(r"^(MF|FN|OZ|ON|LC)\d",k)]
for r in SeqIO.parse(StringIO(ncbi.call("efetch",db="nuccore",id=",".join(gb),rettype="gbwithparts",retmode="text").decode()),"genbank"): rec(r.id,r)
for k in P: out.setdefault(k,None)
json.dump(out,open("results/gate2/panel_features_fixed.json","w")); print(len(out),sum(v is None for v in out.values()),"without annotation")
