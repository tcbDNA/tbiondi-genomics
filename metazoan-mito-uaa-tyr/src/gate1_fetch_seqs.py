"""Gate 1: download organelle (non-RefSeq) FASTA and RefSeq FASTA for records lacking inline sequence. Resumable; logs missing."""
import sys, os, csv, hashlib
sys.path.insert(0,"src"); import ncbi
os.makedirs("data/ncbi/organelle_fasta",exist_ok=True)
A=[r[0] for r in list(csv.reader(open("results/organelle_mito_inventory.tsv"),delimiter="\t"))[1:]]
R=[r for r in csv.DictReader(open("results/refseq_mito_inventory.tsv"),delimiter="\t") if r["seq_in_gbff"]=="no"]
A+= [r["accession"] for r in R]
for i in range(0,len(A),300):
    fn=f"data/ncbi/organelle_fasta/b{i:05d}.fna"
    if os.path.exists(fn) and os.path.getsize(fn)>1000: continue
    t=ncbi.call("efetch",db="nuccore",id=",".join(A[i:i+300]),rettype="fasta",retmode="text").decode()
    open(fn,"w").write(t); print(i,t.count(">"),flush=True)
