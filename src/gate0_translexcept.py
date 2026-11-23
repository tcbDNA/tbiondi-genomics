"""Gate 0: count submitter-annotated transl_except codons (metadata only, named controls)."""
import sys, re, collections, xml.etree.ElementTree as ET
sys.path.insert(0,"src"); import ncbi
A=["OZ281197.1","OZ281198.1","OZ281199.1","OZ281200.1","FN313571.1","MF374802.1","LC597489.1","LC923628.1"]
root=ET.fromstring(ncbi.call("efetch",db="nuccore",id=",".join(A),rettype="gb",retmode="xml"))
rc=lambda s:s.translate(str.maketrans("ACGT","TGCA"))[::-1]
for s in root.iter("GBSeq"):
    seq=s.findtext("GBSeq_sequence").upper(); per={}; cod=collections.Counter()
    for f in s.iter("GBFeature"):
        if f.findtext("GBFeature_key")!="CDS": continue
        g=""; n=0; comp="complement" in f.findtext("GBFeature_location")
        for q in f.iter("GBQualifier"):
            k,v=q.findtext("GBQualifier_name"),q.findtext("GBQualifier_value") or ""
            if k=="gene": g=v
            if k=="transl_except" and "Tyr" in v:
                n+=1; m=re.search(r"(\d+)\.\.(\d+)",v)
                if m: c=seq[int(m[1])-1:int(m[2])]; cod[rc(c) if "complement" in v else c]+=1
        per[g]=n
    print(s.findtext("GBSeq_accession-version"),s.findtext("GBSeq_organism"),"Tyr transl_except total",sum(per.values()),"codons",dict(cod),per)
    for r in s.iter("GBReference"): print("    ",r.findtext("GBReference_title"),"|",r.findtext("GBReference_journal"),"|",[a.text for a in r.iter("GBAuthor")][:6])
