"""Bug fix (post-lock, found after interim look 2): annotation intervals must use each part of a compound location.
Origin-spanning and spliced features were stored as min..max, i.e. as covering the whole genome or the introns. Writes a NEW file; does not touch sweep FASTA or the running sweep.
usage: gate3_features_fix.py out.jsonl [panel]"""
import gzip,glob,json,re,sys,csv,warnings; warnings.filterwarnings("ignore")
from Bio import SeqIO
def parts(f): return [(int(p.start)+1,int(p.end),"+" if p.strand==1 else "-") for p in f.location.parts]
def feats(r):
    d={"cds":[],"trna":[]}
    for f in r.features:
        if f.type=="CDS":
            nm=(f.qualifiers.get("gene") or f.qualifiers.get("product") or ["?"])[0]; d["cds"]+=[[a,b,s,nm] for a,b,s in parts(f)]
        elif f.type=="tRNA":
            t=str(f.qualifiers.get("product",""))+str(f.qualifiers.get("gene","")); nm=(f.qualifiers.get("product") or f.qualifiers.get("gene") or ["?"])[0]; d["trna"]+=[[a,b,s,nm,bool(re.search(r"Tyr|trnY",t,re.I))] for a,b,s in parts(f)]
    return d if d["cds"] or d["trna"] else None
if __name__=="__main__":
    acc={r["acc"] for r in csv.DictReader(open("results/gate3/sweep_audit.tsv"),delimiter="\t") if r["status"]=="included"}; done=set(); o=open(sys.argv[1],"w")
    with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
        for r in SeqIO.parse(fh,"genbank"):
            if r.id in acc: o.write(json.dumps(dict(acc=r.id,f=feats(r)))+"\n"); done.add(r.id)
    for fn in sorted(glob.glob("data/sweep/gb/*.gb")):
        for r in SeqIO.parse(fn,"genbank"):
            if r.id in acc and r.id not in done: o.write(json.dumps(dict(acc=r.id,f=feats(r)))+"\n"); done.add(r.id)
    for a in acc-done: o.write(json.dumps(dict(acc=a,f=None))+"\n")
    print(len(done),"annotated or parsed;",len(acc-done),"without annotation")
