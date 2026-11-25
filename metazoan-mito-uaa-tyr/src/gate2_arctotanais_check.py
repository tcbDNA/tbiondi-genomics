"""Decision 4: check the documented code of LC597489 (Arctotanais) from annotation only."""
import sys,xml.etree.ElementTree as ET; sys.path.insert(0,"src"); import ncbi
for acc in ("LC597489.1","LC923628.1","LC923627.1","LC923626.1","LC923625.1","LC923624.1"):
    root=ET.fromstring(ncbi.call("efetch",db="nuccore",id=acc,rettype="gb",retmode="xml")); s=root.find("GBSeq"); seq=s.findtext("GBSeq_sequence").upper()
    rc=lambda x:x.translate(str.maketrans("ACGT","TGCA"))[::-1]; stops=[];internal=0;n=0
    for f in s.iter("GBFeature"):
        if f.findtext("GBFeature_key")!="CDS": continue
        q={x.findtext("GBQualifier_name"):x.findtext("GBQualifier_value") for x in f.iter("GBQualifier")}
        iv=[(int(i.findtext("GBInterval_from")),int(i.findtext("GBInterval_to"))) for i in f.iter("GBInterval")]
        c="".join(seq[min(a,b)-1:max(a,b)] for a,b in iv) if "complement" not in f.findtext("GBFeature_location") else "".join(rc(seq[min(a,b)-1:max(a,b)]) for a,b in iv)
        off=int(q.get("codon_start","1"))-1; c=c[off:]; n+=1
        cod=[c[i:i+3] for i in range(0,len(c)-len(c)%3,3)]; internal+=sum(x in("TAA","TAG") for x in cod[:-1]); stops.append((q.get("gene"),cod[-1] if len(c)%3==0 else cod[-1]+"+%d"%(len(c)%3), "transl_except" if "transl_except" in q else "",q.get("transl_table")))
    print(acc,s.findtext("GBSeq_organism"),"CDS",n,"internal annotated TAA/TAG:",internal,"\n   terminal:",stops)
