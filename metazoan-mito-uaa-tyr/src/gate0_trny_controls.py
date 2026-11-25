"""Gate 0 only: locate the GUA anticodon in trnY of NAMED CONTROL accessions and report local vs canonical index.
Annotated tRNA features come from GenBank; unannotated controls use tRNAscan-SE -O and ARAGORN.
The allowlist enforces the pre-registration firewall."""
import sys, re, subprocess, os, xml.etree.ElementTree as ET
sys.path.insert(0,"src"); import ncbi
CONTROLS = {"PP708078.1":"Anomia chinensis (+)","OZ576062.1":"Anomia ephippium (+)","PP708091.1":"Placuna placenta (+)",
 "OZ181722.1":"Heteranomia squamula (-)","MF374802.1":"Cephalodiscus hodgsoni (+)","FN908482.1":"Rhabdopleura compacta (-)",
 "FN313571.1":"Radopholus similis (+)","OZ281197.1":"Pseudotanais sp. (+)","LC597489.1":"Arctotanais alascensis (-)","LC923628.1":"Sinelobus kisui (-)"}
os.makedirs("data/ncbi/controls",exist_ok=True)
root=ET.fromstring(ncbi.call("efetch",db="nuccore",id=",".join(CONTROLS),rettype="gb",retmode="xml"))
rc=lambda s:s.translate(str.maketrans("ACGT","TGCA"))[::-1]
for s in root.iter("GBSeq"):
    acc=s.findtext("GBSeq_accession-version"); seq=s.findtext("GBSeq_sequence").upper()
    open(f"data/ncbi/controls/{acc}.fna","w").write(f">{acc}\n{seq}\n")
    print("##",acc,CONTROLS[acc],len(seq)); hit=0
    for f in s.iter("GBFeature"):
        if f.findtext("GBFeature_key")!="tRNA": continue
        q={x.findtext("GBQualifier_name"):(x.findtext("GBQualifier_value") or "") for x in f.iter("GBQualifier")}
        if not re.search(r"tyr|trnY",q.get("product","")+q.get("gene","")+q.get("note",""),re.I): continue
        loc=f.findtext("GBFeature_location"); ivs=list(f.iter("GBInterval"))
        t="".join(seq[min(int(i.findtext("GBInterval_from")),int(i.findtext("GBInterval_to")))-1:max(int(i.findtext("GBInterval_from")),int(i.findtext("GBInterval_to")))] for i in ivs)
        if "complement" in loc: t=rc(t)
        hit=1; sites=[m.start() for m in re.finditer("(?=GTA)",t)]
        print("  annotated",loc,len(t),t,"anticodon_q=",q.get("anticodon",""),"| GTA starts (1-based)",[(i+1,t[i-1] if i else "-",t[max(0,i-5):i]+"["+t[i:i+3]+"]"+t[i+3:i+8]) for i in sites])
    if not hit: print("  no annotated trnY")
for a in CONTROLS:
    fa=f"data/ncbi/controls/{a}.fna"
    r=subprocess.run(["tRNAscan-SE","-O","-q","--brief",fa],capture_output=True,text=True); 
    print("tRNAscan",a,[l.split()[:9] for l in r.stdout.splitlines() if "Tyr" in l or "GTA" in l], r.stderr.strip()[:200])
    r=subprocess.run(["aragorn","-mt","-gc5","-w",fa],capture_output=True,text=True)
    print("aragorn ",a,[l for l in r.stdout.splitlines() if re.search("Tyr",l)])
