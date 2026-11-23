"""Gate 1: non-RefSeq metazoan organelle assemblies (DToL/ATLASea/etc.). Metadata via esummary + gb records; no caller is run."""
import sys, csv, re, json, xml.etree.ElementTree as ET
sys.path.insert(0,"src"); import ncbi
q='Metazoa[Organism] AND "genome assembly"[Title] AND mitochondrion[filter] NOT refseq[filter]'
d=json.loads(ncbi.call("esearch",db="nuccore",term=q,retmode="json",retmax=0,usehistory="y"))["esearchresult"]
N=int(d["count"]); we,qk=d["webenv"],d["querykey"]; print("records",N,flush=True)
out=csv.writer(open("results/organelle_mito_inventory.tsv","w"),delimiter="\t")
out.writerow("accession organism taxid length topology definition bioproject dtol create_date n_cds comment_flag".split())
for s0 in range(0,N,200):
    x=ncbi.call("efetch",db="nuccore",query_key=qk,WebEnv=we,retstart=s0,retmax=200,rettype="gb",retmode="xml",complexity=1)
    root=ET.fromstring(x)
    for s in root.iter("GBSeq"):
        com=s.findtext("GBSeq_comment") or ""; d_=s.findtext("GBSeq_definition") or ""
        tx=next((q.findtext("GBQualifier_value")[6:] for q in s.iter("GBQualifier") if q.findtext("GBQualifier_name")=="db_xref" and q.findtext("GBQualifier_value").startswith("taxon:")),"")
        nc=sum(1 for f in s.iter("GBFeature") if f.findtext("GBFeature_key")=="CDS")
        prj=s.findtext("GBSeq_project") or ""
        dt="DToL" if "Darwin Tree of Life" in com else "Sanger/ToL-other" if "Tree of Life" in com else "ATLASea" if "ATLASea" in com else "other"
        out.writerow([s.findtext("GBSeq_accession-version"),s.findtext("GBSeq_organism"),tx,s.findtext("GBSeq_length"),s.findtext("GBSeq_topology"),d_,prj,dt,s.findtext("GBSeq_create-date"),nc,"Oatk" if "Oatk" in com else "MBG" if "MBG" in com else "other"])
    print(s0+200,flush=True)
