"""Gate 0: print references, comments and per-CDS transl_table/transl_except for candidate control records."""
import sys, xml.etree.ElementTree as ET
sys.path.insert(0,"src"); import ncbi
accs = sys.argv[1:]
root = ET.fromstring(ncbi.call("efetch", db="nuccore", id=",".join(accs), rettype="gb", retmode="xml"))
for s in root.iter("GBSeq"):
    print("==", s.findtext("GBSeq_accession-version"), s.findtext("GBSeq_organism"), s.findtext("GBSeq_length"), s.findtext("GBSeq_topology"), "|", s.findtext("GBSeq_definition"))
    print("  project:", s.findtext("GBSeq_project"), "| keywords:", [k.text for k in s.iter("GBKeyword")])
    for r in s.iter("GBReference"): print("  ref:", r.findtext("GBReference_title"), "|", r.findtext("GBReference_journal"), "|", (r.findtext("GBReference_pubmed") or ""))
    if s.findtext("GBSeq_comment"): print("  comment:", s.findtext("GBSeq_comment")[:400])
    n=0
    for f in s.iter("GBFeature"):
        k=f.findtext("GBFeature_key"); q={x.findtext("GBQualifier_name"):x.findtext("GBQualifier_value") for x in f.iter("GBQualifier")}
        if k=="CDS":
            n+=1
            if any(y in q for y in ("transl_except","exception","note","pseudo")): print("  CDS",q.get("gene") or q.get("product"),{y:q[y] for y in ("transl_table","transl_except","exception","note") if y in q})
    print("  n_CDS",n)
