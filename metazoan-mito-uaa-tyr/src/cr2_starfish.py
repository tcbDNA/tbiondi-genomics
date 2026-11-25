"""Cold review 2, check 3 (POST-HOC): the starfish and hemichordate accessions in Li et al. 2019 (supplementary table S2; figures S2, S3) against our calls; and the base 33 of every annotated tRNA in the two starfish the paper cites, to look for a source of the 'starfish C' statement.
writes results/posthoc2/starfish_li2019.txt and results/posthoc2/starfish_asteroids.tsv"""
import csv,gzip,json,re,warnings; warnings.filterwarnings("ignore")
from Bio import SeqIO
M={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")}; o=open("results/posthoc2/starfish_li2019.txt","w"); W=lambda *x:(print(*x),print(*x,file=o))
G={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_trny_calls.tsv"),delimiter="\t")}; L={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/trny_panel_final.jsonl")}
W("POST-HOC. Li et al. 2019, Genome Biol Evol 11:154, PMC6319601; supplement evy254_supp.docx (PMC open-data copy).")
W("Figure S2 (alignment of tRNA-Tyr) rows: G. marginatus, St. canadensis, C. hodgsoni, B. carnosus, Sch. brasiliense, B. clavigerus, R. compacta. Figure S3 (structures): the same plus Sa. kowalevskii and a hemichordate consensus. Both captions say 'in hemichordate mitochondrial genomes'. NEITHER CONTAINS A STARFISH.")
W("\nAccessions in table S2 with a full mitogenome, and our call:")
for a,sp in (("NC_007788","Acanthaster planci (Asteroidea; the only starfish)"),("NC_001453","Strongylocentrotus purpuratus (Echinoidea)"),("NC_005929","Cucumaria miniata (Holothuroidea)"),("NC_001887","Balanoglossus carnosus"),("MG652439","Glossobalanus marginatus"),("NC_040109","Glossobalanus marginatus RefSeq"),("NC_013877","Balanoglossus clavigerus"),("NC_007438","Saccoglossus kowalevskii"),("MH841935","Stereobalanus canadensis"),("MH841936","Schizocardium brasiliense"),("MF374802","Cephalodiscus hodgsoni"),("NC_015649","Rhabdopleura compacta"),("FN908482","Rhabdopleura compacta INSDC")):
    x=[m for k,m in M.items() if k.split(".")[0]==a]; y=[k for k in L if k.split(".")[0]==a]
    if x: m=x[0]; W("  ",m["acc"],sp,"role",m["role"],"UAA",m["uaa"],"tRNA-Tyr",m["n33"],m["base"]+"33","loop",m["loop"])
    elif y: d=L[y[0]]; W("  ",y[0],sp,"panel; gate",G[y[0]]["final"],"tRNA-Tyr",d["N33"]+"33","loop",d["loop"])
    else: W("  ",a,sp,"not in sweep or panel (duplicate of a RefSeq/INSDC record, or absent)")
W("\nOther starfish sources cited in the paper's text: Asakawa et al. 1995 (Asterina pectinifera = Patiria pectinifera, D16387 = NC_001627); Hyouta et al. 1987 and Matsuyama et al. 1998 (starfish tRNA-Ser(GCU)); Tomita et al. 1999 (echinoderm tRNA-Asn(GUU) 'with C33').")
A=[m for m in M.values() if m["class"]=="Asteroidea" and m["status"]=="included"]; w=csv.writer(open("results/posthoc2/starfish_asteroids.tsv","w"),delimiter="\t"); w.writerow("acc species family provenance uaa n33 base loop score stem n_agree reason".split())
for m in sorted(A,key=lambda m:m["species"]): w.writerow([m[k] for k in "acc species family provenance uaa n33 base loop score stem n_agree n33_reason".split()])
W("\nAll asteroids in the sweep:",len(A),"resolved U33",sum(m["n33"]=="resolved" and m["base"]=="U" for m in A),"non-U33",sum(m["base"] not in("U","") for m in A),"unresolved:",[(m["acc"],m["species"],m["base"],m["n33_reason"]) for m in A if m["n33"]!="resolved"])
for k in ("Acanthaster planci","Patiria pectinifera","Asterias amurensis","Asterias rubens"): [W("  ",m["acc"],m["species"],m["n33"],m["base"],m["loop"],"score",m["score"]) for m in A if m["species"]==k]
W("\nBase 5' of the annotated/located anticodon in every annotated tRNA (position 33), starfish records with annotation:")
want={"NC_007788.1","NC_001627.1","NC_006665.1"}; AC={"Asn":"GTT","Tyr":"GTA","Ser":None,"His":"GTG","Asp":"GTC","Cys":"GCA","Phe":"GAA"}
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    for r in SeqIO.parse(fh,"genbank"):
        if r.id not in want: continue
        res=[]
        for f in r.features:
            if f.type!="tRNA": continue
            pr=(f.qualifiers.get("product") or [""])[0]; aa=pr.replace("tRNA-",""); s=str(f.extract(r.seq)).upper(); n=f.qualifiers.get("note",[""])[0]
            for ac in ([AC[aa]] if AC.get(aa) else ["GCT","TGA"] if aa=="Ser" else []):
                c=[m.start() for m in re.finditer(f"(?={ac})",s[24:42])]; c=[24+i for i in c if 24+i>=28]
                if c: i=min(c,key=lambda i:abs(i-33)); res.append(f"{aa}({ac}): pos33={s[i-1].replace('T','U')} loop={s[i-2:i+5].replace('T','U')}")
        W("  ",r.id,r.annotations["organism"],"; ".join(res))
