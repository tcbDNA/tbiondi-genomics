"""Gate 1: inventory RefSeq + organelle-assembly metazoan mitogenomes. Reads annotation metadata only; runs no caller.
Outputs results/refseq_mito_inventory.tsv. Every excluded record is logged with a reason."""
import gzip, sys, csv, re
from Bio import SeqIO
src="data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz"
out=csv.writer(open("results/refseq_mito_inventory.tsv","w"),delimiter="\t"); aud=csv.writer(open("results/audit/gate1_refseq_excluded.tsv","w"),delimiter="\t")
H="accession organism taxid phylum lineage length topology completeness date refseq_status derived_from n_cds n_trna n_rrna transl_tables n_ambiguous seq_in_gbff".split()
out.writerow(H); aud.writerow(["accession","organism","reason"])
n=0
with gzip.open(src,"rt") as f:
    for r in SeqIO.parse(f,"genbank"):
        n+=1; lin=r.annotations.get("taxonomy",[]); org=r.annotations.get("organism","")
        if "Metazoa" not in lin: aud.writerow([r.id,org,"not_metazoa"]); continue
        tax=""; tt=set(); c=t=rr=0
        for ft in r.features:
            if ft.type=="source":
                for x in ft.qualifiers.get("db_xref",[]):
                    if x.startswith("taxon:"): tax=x[6:]
            elif ft.type=="CDS": c+=1; tt.update(ft.qualifiers.get("transl_table",[]))
            elif ft.type=="tRNA": t+=1
            elif ft.type=="rRNA": rr+=1
        com=r.annotations.get("comment",""); m=re.search(r"(?:derived from|identical to) ([A-Z]{1,6}\d+)",com)
        st=re.match(r"(\w+) REFSEQ",com); d=r.description
        comp="complete" if "complete genome" in d else "partial" if "partial" in d else "unstated"
        try: s=str(r.seq).upper(); sd="yes"
        except Exception: s=""; sd="no"; L0=len(r.seq)
        ph=next((p for p in ("Porifera","Cnidaria","Ctenophora","Placozoa") if p in lin),"") 
        out.writerow([r.id,org,tax,"",";".join(lin),(len(s) or L0),r.annotations.get("topology",""),comp,r.annotations.get("date",""),st.group(1) if st else "",m.group(1) if m else "",c,t,rr,",".join(sorted(tt,key=int)),sum(s.count(x) for x in "NRYKMSWBDHV"),sd])
print("records parsed",n)
