"""Corrected read support for the candidates: record k-mer vs tRNA-U33 k-mer (31-mers centred on N33), over the species' public runs (largest first capped at MAXGB). usage: candidate_reads.py ACC [RUN...]"""
import sys,os,csv,json; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import kmer_support as KS, ncbi
acc=sys.argv[1]; M={r["acc"]:r for r in csv.DictReader(open("results/gate3/sweep_audit.tsv"),delimiter="\t")}[acc]; p="data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna"; s="".join(l.strip() for l in open(p) if l[0]!=">").upper(); d=None
for l in open("results/gate3/run/out.jsonl"):
    x=json.loads(l)
    if x["acc"]==acc: d=x["trny"]
g=int(d["g34_pos"]); site=(g-2) if d["strand"]!="-" else g; a,b=KS.variants(s,site,d["strand"]); print(acc,M["species"],"N33",d["N33"],"record",a,"U33",b)
runs=sys.argv[2:]
for r in runs:
    c,e=KS.support(r,{"record":[a],"U33":[b]},maxgb=25); print("  run",r,dict(c) if c is not None else e,flush=True)
