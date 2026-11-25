"""Gate 0: list public mitochondrial records for candidate control taxa (metadata only)."""
import sys, csv, xml.etree.ElementTree as ET
sys.path.insert(0, "src"); import ncbi
TAXA = ["Anomia","Enigmonia","Placuna","Placunidae","Plicatula","Plicatulidae","Pododesmus","Heteranomia","Monia",
 "Cephalodiscus","Rhabdopleura","Radopholus similis","Pseudotanais","Arctotanais alascensis","Tanaidacea"]
out = csv.writer(open("results/audit/gate0_control_accessions.tsv","w"), delimiter="\t")
out.writerow("query accession organism taxid length topology completeness moltype transl_tables definition create_date".split())
for t in TAXA:
    term = f'"{t}"[Organism] AND mitochondrion[filter] AND 5000:200000[SLEN]'
    n, ids, _ = ncbi.search("nuccore", term, 200)
    print(t, n)
    if not ids: continue
    root = ET.fromstring(ncbi.call("efetch", db="nuccore", id=",".join(ids), rettype="gb", retmode="xml"))
    for s in root.iter("GBSeq"):
        tt = sorted({q.findtext("GBQualifier_value") for q in s.iter("GBQualifier") if q.findtext("GBQualifier_name")=="transl_table"})
        tx = next((q.findtext("GBQualifier_value")[6:] for q in s.iter("GBQualifier") if q.findtext("GBQualifier_name")=="db_xref" and q.findtext("GBQualifier_value").startswith("taxon:")),"")
        d = s.findtext("GBSeq_definition")
        comp = "complete" if "complete genome" in d else "partial" if "partial" in d else "unstated"
        row=[t,s.findtext("GBSeq_accession-version"),s.findtext("GBSeq_organism"),tx,s.findtext("GBSeq_length"),s.findtext("GBSeq_topology"),comp,s.findtext("GBSeq_moltype"),",".join(tt),d,s.findtext("GBSeq_create-date")]
        out.writerow(row); print("  ",*row[1:7],row[8],row[10],sep=" | ")
