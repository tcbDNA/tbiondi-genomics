"""Gate 1: size the INSDC-only invertebrate mitogenomes absent from RefSeq + organelle set. Metadata only."""
import sys,json,csv; sys.path.insert(0,"src"); import ncbi
base='Metazoa[Organism] NOT Chordata[Organism] AND mitochondrion[filter] AND "complete genome"[Title] AND 8000:200000[SLEN] NOT refseq[filter] NOT "genome assembly"[Title]'
d=json.loads(ncbi.call("esearch",db="nuccore",term=base,retmode="json",retmax=0,usehistory="y"))["esearchresult"]
we,qk,N=d["webenv"],d["querykey"],int(d["count"]); rows=[]
for s in range(0,N,500):
    for k in range(4):
        try: j=json.loads(ncbi.call("esummary",db="nuccore",query_key=qk,WebEnv=we,retstart=s,retmax=500,retmode="json")); break
        except Exception as e: print("retry",s,e,flush=True)
    r=j.get("result",{}); 
    for i in r.get("uids",[]):
        x=r[i]; rows.append((x.get("accessionversion"),x.get("taxid"),x.get("slen"),x.get("topology"),x.get("title")))
w=csv.writer(open("results/insdc_invert_summaries.tsv","w"),delimiter="\t"); w.writerow("accession taxid length topology title".split()); w.writerows(rows)
print("expected",N,"got",len(rows))
