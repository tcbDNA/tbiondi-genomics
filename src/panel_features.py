"""Annotated CDS and tRNA intervals for panel records (annotation only)."""
import csv,gzip,json,re,sys,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"src"); import panel_guard, ncbi
from Bio import SeqIO; from io import StringIO
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; panel_guard.check(P); out={}
def rec(pid,r):
    d={"cds":[],"trna":[]}
    for f in r.features:
        if f.type=="CDS": d["cds"].append([int(f.location.start)+1,int(f.location.end),"+" if f.location.strand==1 else "-",(f.qualifiers.get("gene") or ["?"])[0]])
        elif f.type=="tRNA":
            t=str(f.qualifiers.get("product",""))+str(f.qualifiers.get("gene","")); d["trna"].append([int(f.location.start)+1,int(f.location.end),"+" if f.location.strand==1 else "-",(f.qualifiers.get("product") or ["?"])[0],bool(re.search(r"Tyr|trnY",t,re.I))])
    out[pid]=d if (d["cds"] or d["trna"]) else None
want={k for k in P if k.startswith("NC_")}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id in want: rec(r.id,r)
gb=[k for k in P if re.match(r"^(MF|FN|OZ|ON|LC)\d",k)]
for r in SeqIO.parse(StringIO(ncbi.call("efetch",db="nuccore",id=",".join(gb),rettype="gbwithparts",retmode="text").decode()),"genbank"): rec(r.id,r)
for k in P: out.setdefault(k,None)
json.dump(out,open("results/gate2/panel_features.json","w")); print(len(out),sum(v is None for v in out.values()),"without annotation")
