"""Compare extractor calls on the broad panel with RefSeq-annotated trnY coordinates (independent of the CM score).
Outputs results/work/trny_concordance.tsv"""
import csv,gzip,json,re,collections
from Bio import SeqIO
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
D={}
for l in open("results/work/trny_panel.jsonl"): d=json.loads(l); D[d["panel_id"]]=d
want={k for k,v in P.items() if v["role"]=="neg_broad"}; ann={}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id in want:
            L=[]
            for f in r.features:
                if f.type=="tRNA" and re.search(r"Tyr|trnY",str(f.qualifiers.get("product",""))+str(f.qualifiers.get("gene","")),re.I):
                    L.append((int(f.location.start)+1,int(f.location.end),"+" if f.location.strand==1 else "-"))
            ann[r.id]=L
w=csv.writer(open("results/work/trny_concordance.tsv","w"),delimiter="\t"); w.writerow("panel_id lineage status reason score evalue N33 n_annotated overlap".split()); out=[]
for k in want:
    d=D[k]; A=ann.get(k,[]); ov=""
    if d["status"]=="called":
        s,e=sorted((int(d["start"]),int(d["end"]))); ov="yes" if any(min(e,b)-max(s,a)>0.5*min(e-s,b-a) for a,b,st in A) else ("no" if A else "NA_not_annotated")
    w.writerow([k,P[k]["lineage"],d["status"],d.get("reason"),d.get("score"),d.get("evalue"),d.get("N33"),len(A),ov]); out.append((k,P[k]["lineage"],d,len(A),ov))
called=[x for x in out if x[2]["status"]=="called"]
print("called",len(called),collections.Counter((x[2]["reason"],x[4]) for x in called))
for thr in (0,5,10,12,14,15,16,18,20,22,25,30):
    c=[x for x in called if x[2]["score"]>=thr]; a=[x for x in c if x[3]>0]; 
    print("score>=",thr,"accepted",len(c),"| with annotation",len(a),"overlap_yes",sum(x[4]=="yes" for x in a),"overlap_no",sum(x[4]=="no" for x in a),"| unannotated",sum(x[3]==0 for x in c))
print("\noverlap=no calls:",[(x[0],x[1],x[2]["score"],x[2]["N33"],x[2]["reason"]) for x in called if x[4]=="no"])
print("non-U33 calls:",[(x[0],x[1],x[2]["score"],x[2]["N33"],x[4]) for x in called if x[2]["N33"]!="U"])
print("annotated trnY genomes total",sum(1 for x in out if x[3]>0),"by lineage none-annotated:",collections.Counter(x[1] for x in out if x[3]==0))
