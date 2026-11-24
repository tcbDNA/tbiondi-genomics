"""EXPLORATORY / POST-HOC: GenBank flat files for the INSDC-only species excluded as UNVERIFIED. Separate from Gate 3. Resumable; ledgered; 30 GB floor."""
import sys,os,csv,hashlib,re,shutil; sys.path.insert(0,"src"); import ncbi
A=[r["accession"] for r in csv.DictReader(open("results/insdc_only_set.tsv"),delimiter="\t") if r["eligible"]=="no" and "unverified" in r["title"].lower()]
os.makedirs("data/sweep_unverified/gb",exist_ok=True); L=open("data/manifests/download_ledger.tsv","a"); miss=[]
for i in range(0,len(A),100):
    f=f"data/sweep_unverified/gb/b{i//100:04d}.gb"; want=A[i:i+100]
    if os.path.exists(f): continue
    if shutil.disk_usage(".").free<32e9: print("STOP: free disk near 30 GB",flush=True); sys.exit(2)
    t="";got=set()
    for k in range(5):
        try:
            t=ncbi.call("efetch",db="nuccore",id=",".join(want),rettype="gbwithparts",retmode="text").decode(errors="replace"); got=set(re.findall(r"^VERSION\s+(\S+)",t,re.M))
            if len(got&set(want))==len(want): break
        except Exception as e: print("retry",i,str(e)[:80],flush=True)
    m=[a for a in want if a not in got]; miss+=m; open(f+".tmp","w").write(t); os.replace(f+".tmp",f); b=os.path.getsize(f)
    L.write(f"3x\tNCBI efetch gbwithparts UNVERIFIED batch {i//100} ({len(want)} accessions, {len(m)} missing)\t{f}\t{b}\t{hashlib.sha256(open(f,'rb').read()).hexdigest()}\ton_disk\n"); L.flush(); print(f"{i+len(want)}/{len(A)} missing {len(m)}",flush=True)
open("results/audit/gate3u_fetch_missing.tsv","w").write("accession\n"+"\n".join(miss)+"\n"); print("FETCH DONE",len(A),"missing",len(miss),flush=True)
