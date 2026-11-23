"""Gate 3 data: GenBank flat files (sequence + annotation) for INSDC-only eligible species, annotated organelle representatives and the named T4/T6 extras.
Resumable; each batch checked for completeness, hashed and written to the download ledger. No caller is run here."""
import sys,os,csv,hashlib,re,shutil; sys.path.insert(0,"src"); import ncbi
I=[r["accession"] for r in csv.DictReader(open("results/insdc_only_set.tsv"),delimiter="\t") if r["eligible"]=="yes"]
O=[r["accession"] for r in csv.DictReader(open("results/gate1_species.tsv"),delimiter="\t") if r["set"]!="RefSeq" and int(r["n_cds"])>0]
X=["AP002973.1","LR745847.1"]
A=I+[a for a in O if a not in set(I)]+[a for a in X if a not in set(I)|set(O)]
open("data/sweep/accessions.txt","w").write("\n".join(A)+"\n"); L=open("data/manifests/download_ledger.tsv","a"); miss=[]
for i in range(0,len(A),100):
    f=f"data/sweep/gb/b{i//100:04d}.gb"; want=A[i:i+100]
    if os.path.exists(f): continue
    if shutil.disk_usage(".").free<32e9: print("STOP: free disk near 30 GB",flush=True); sys.exit(2)
    for k in range(5):
        try:
            t=ncbi.call("efetch",db="nuccore",id=",".join(want),rettype="gbwithparts",retmode="text").decode(errors="replace")
            got=set(re.findall(r"^VERSION\s+(\S+)",t,re.M))
            if len(got&set(want))==len(want): break
        except Exception as e: print("retry",i,str(e)[:80],flush=True); t="";got=set()
    m=[a for a in want if a not in got]; miss+=m
    open(f+".tmp","w").write(t); os.replace(f+".tmp",f); b=os.path.getsize(f); h=hashlib.sha256(open(f,"rb").read()).hexdigest()
    L.write(f"3\tNCBI efetch nuccore gbwithparts batch {i//100} ({len(want)} accessions, {len(m)} missing)\t{f}\t{b}\t{h}\ton_disk\n"); L.flush()
    print(f"{i+len(want)}/{len(A)} missing {len(m)}",flush=True)
open("results/audit/gate3_fetch_missing.tsv","w").write("accession\n"+"\n".join(miss)+"\n"); print("FETCH DONE",len(A),"missing",len(miss),flush=True)
