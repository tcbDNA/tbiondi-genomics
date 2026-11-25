"""Summarise Lin Figshare mitogenome files: length, topology, CDS, tables, stop-type notes. Metadata only (no caller)."""
import glob,os,re
from Bio import SeqIO
import warnings; warnings.filterwarnings("ignore")
for f in sorted(glob.glob("data/figshare/mito/*")):
    b=os.path.basename(f)
    if f.endswith(".gb"):
        try:
            for r in SeqIO.parse(f,"genbank"):
                cds=[x for x in r.features if x.type=="CDS"]; trn=[x for x in r.features if x.type=="tRNA"]
                tt={t for x in cds for t in x.qualifiers.get("transl_table",[])}
                print(b,"|",r.id,"|",r.annotations.get("organism"),"|",len(r.seq),r.annotations.get("topology"),"| CDS",len(cds),sorted((x.qualifiers.get("gene") or x.qualifiers.get("product") or ["?"])[0] for x in cds)[:16],"tRNA",len(trn),"tables",tt)
        except Exception as e: print(b,"PARSE FAIL",e)
    else:
        rs=list(SeqIO.parse(f,"fasta")); print(b,"| FASTA",[(x.id,len(x.seq)) for x in rs])
