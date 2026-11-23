"""Annotated terminal codons for panel genomes (RefSeq gbff, GenBank records, Lin .gb). Annotation only."""
import csv,gzip,json,re,sys,glob,os,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,"src"); import panel_guard, ncbi
from Bio import SeqIO; from io import StringIO
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; panel_guard.check(P)
NM=[("COX1",r"^cox1$|^coxi$|^co1$|^coi$|cytochrome c oxidase subunit (1|i)$"),("COX2",r"^cox2$|^coii$|^co2$|^coxii$|oxidase subunit (2|ii)$"),("COX3",r"^cox3$|^coiii$|^co3$|^coxiii$|oxidase subunit (3|iii)$"),("CYTB",r"^cob$|^cytb$|^cyt ?b$|^cytochrome b$"),("ATP6",r"^atp6$|^atpase ?6$|synthase (f0 )?subunit 6"),("ATP8",r"^atp8$|^atpase ?8$|synthase (f0 )?subunit 8")]+[("ND"+n,rf"^n(a?d|adh) ?{n}$|dehydrogenase subunit {n}$") for n in("1","2","3","4","5","6","4L")]
def gname(q):
    for k in("gene","label","product","standard_name"):
        for v in q.get(k,[]):
            for g,p in NM:
                if re.search(p,v.strip(),re.I): return g
    return None
def term(f,s):
    try: c=str(f.extract(s)).upper()
    except Exception: return None
    c=c[int((f.qualifiers.get("codon_start") or ["1"])[0])-1:]; r=len(c)%3
    if r: return {"T":"T--","TA":"TA-"}.get(c[-r:],"none")
    return c[-3:] if c[-3:] in("TAA","TAG") else "none"
out={}
def rec(pid,r):
    d={}
    for f in r.features:
        if f.type=="CDS":
            g=gname(f.qualifiers); t=term(f,r.seq)
            if g and t and g not in d: d[g]=t
    out[pid]=d
want={k for k in P if k.startswith("NC_")}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id in want: rec(r.id,r)
gb=[k for k in P if re.match(r"^(MF|FN|OZ|ON|LC)\d",k)]
for r in SeqIO.parse(StringIO(ncbi.call("efetch",db="nuccore",id=",".join(gb),rettype="gbwithparts",retmode="text").decode()),"genbank"): rec(r.id,r)
for f in glob.glob("data/figshare/mito/*.gb"):
    pid="LIN_"+re.sub(r"\W+","_",os.path.basename(f)[:-3])
    if pid in P:
        try: rec(pid,next(SeqIO.parse(f,"genbank")))
        except Exception as e: print("no annotation parsed",pid,str(e)[:60])
json.dump(out,open("results/gate2/annot_terminals.json","w")); print(len(out),"genomes;",sum(len(v) for v in out.values()),"genes")
