"""EXPLORATORY: weak-tier audit (locked gate3_post.weak_audit) for one UNVERIFIED accession, relatives = Gate 3 Stop species. usage: gate3u_weak.py ACC"""
import sys,json,csv; sys.path.insert(0,"src"); import gate3_classify as K, gate3_post as P
acc=sys.argv[1]; M={r["acc"]:r for r in csv.DictReader(open("results/gate3_unverified/sweep_audit.tsv"),delimiter="\t")}; PA={r["acc"]:r["path"] for r in csv.DictReader(open("results/gate3_unverified/manifest.tsv"),delimiter="\t")}
d=[x for x in map(json.loads,open("results/gate3_unverified/run/out.jsonl")) if x["acc"]==acc][0]; r=dict(M[acc],path=PA[acc],uaa="weak")
S=[q for q in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if q["status"]=="included" and q["provenance"]!="named_extra" and q["uaa"]=="Stop"]
u=d["uaa"]; print(acc,r["species"],r["phylum"],r["class"],r["order"],r["family"],"len",r["length"],"class",K.uaa_class(d),"| TAA",u["n_taa_internal"],"genes",u["n_taa_genes"],"tyr",u["n_taa_tyrcol"],"LLR",u["llr_Y"],"TAG",u["n_tag_internal"],"codetta",d["codetta"].get("TAA"),"aligned",u["aligned_codons"],"genes_hit",u["n_genes_hit"])
G=[q for q in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if q["status"]=="included" and q["uaa"] in("Tyr","weak")]; print("audit (others = Gate 3 Tyr/weak incl. panel, as in gate3_post)",P.weak_audit(r,S,G+[r])); m,tag,fr=P.cols(r); rel=P.closest(r,S); print("relatives",[(k,x["species"]) for k,x in rel]); print("genes with >1 frame",[g for g,f in fr.items() if len(f)>1],"TAG",tag)
for k in [k for k,c in m.items() if c=="TAA"]: print("site",k,[P.cols(x)[0].get(k) for _,x in rel])
t=d["trny"]; print("trnY",{k:t.get(k) for k in("status","reason","score","N33","loop","stem_pairs")},K.n33(d,None)["reason"])
