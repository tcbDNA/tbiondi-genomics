"""Gate 1 extension: all non-RefSeq INSDC metazoan mitochondrial records 10-200 kb (metadata only, resumable, no title filter)."""
import sys,json,csv,os,time; sys.path.insert(0,"src"); import ncbi
Q=[("invert",'Metazoa[Organism] NOT Chordata[Organism] AND mitochondrion[filter] AND 10000:200000[SLEN] NOT refseq[filter] NOT "genome assembly"[Title]'),
   ("chordate",'Chordata[Organism] NOT "Homo sapiens"[Organism] AND mitochondrion[filter] AND 10000:200000[SLEN] NOT refseq[filter] NOT "genome assembly"[Title]'),
   ("invert_gatitle_shorttopup",'Metazoa[Organism] NOT Chordata[Organism] AND mitochondrion[filter] AND 8000:9999[SLEN] AND complete[Title] NOT refseq[filter]')]
os.makedirs("results/work",exist_ok=True)
for name,q in Q:
    out=f"results/work/insdc_{name}.tsv"
    done=sum(1 for _ in open(out))-1 if os.path.exists(out) else 0
    d=json.loads(ncbi.call("esearch",db="nuccore",term=q,retmode="json",retmax=0,usehistory="y"))["esearchresult"]; N=int(d["count"]); print(name,N,"already",done,flush=True)
    f=open(out,"a"); w=csv.writer(f,delimiter="\t")
    if done<=0: w.writerow("accession taxid length topology title create_date".split()); done=max(done,0)
    for s in range((done//500)*500,N,500):
        for k in range(6):
            try: r=json.loads(ncbi.call("esummary",db="nuccore",query_key=d["querykey"],WebEnv=d["webenv"],retstart=s,retmax=500,retmode="json"))["result"]; break
            except Exception as e: print("retry",s,str(e)[:80],flush=True); time.sleep(5)
        else: print("FAILED_BATCH",name,s,flush=True); continue
        for i in r.get("uids",[]):
            x=r[i]; w.writerow([x.get("accessionversion"),x.get("taxid"),x.get("slen"),x.get("topology"),x.get("title"),x.get("createdate")])
        f.flush()
    print(name,"done",flush=True)
