"""Cold review 2, check 7 / item 17 (POST-HOC): which included records are not whole mitogenomes by title, topology and gene content. Reads titles from results/insdc_only_set.tsv, organelle inventory, RefSeq gbff definitions.
Classes (first that applies):  partial_title (partial / nearly complete / near-complete); gene_cds_title (names genes or 'cds' and has no whole-genome phrase); segment_title (chromosome/minicircle/minichromosome/segment/circle of a multipartite genome);
 no_complete_claim (circular or assembly record without 'complete' or a whole-genome phrase); whole_title.  writes results/posthoc2/title_audit.tsv, title_audit_summary.txt"""
import csv,gzip,re,collections,json,warnings; warnings.filterwarnings("ignore")
M=[r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included"]
T={}
for r in csv.DictReader(open("results/insdc_only_set.tsv"),delimiter="\t"): T[r["accession"]]=r["title"]
for r in csv.DictReader(open("results/organelle_mito_inventory.tsv"),delimiter="\t"): T[r["accession"]]=r["definition"]
with gzip.open("data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz","rt") as fh:
    d="";v="";ind=False
    for l in fh:
        if l.startswith("DEFINITION"): d=l[12:].strip(); ind=True
        elif ind and l.startswith("            "): d+=" "+l.strip()
        elif l.startswith("ACCESSION"): ind=False
        elif l.startswith("VERSION"): T[l.split()[1]]=d
FE={}
for l in open("data/sweep/features_fixed.jsonl"):
    x=json.loads(l); FE[x["acc"]]=x["f"]
WH=re.compile(r"complete (mitochondrial |mitochondrion )?genome|complete mitochondrial (dna|sequence)|mitochondrion,? complete sequence|complete sequence|mitochondrial genome|mitochondrion genome|genome assembly, organelle|organelle: mitochondrion|whole genome|mitochondrion,? genome|mitogenome",re.I)
SEG=re.compile(r"chromosome|minicircle|mini-?chromosome|segment\b|\bcircle\b|mtdna[ -]?(i|ii|1|2)\b|contig: ",re.I)
def cls(t):
    tl=t.lower()
    if re.search(r"partial|nearly complete|near-complete|near complete",tl): return "partial_title"
    if re.search(r"\bcds\b|\bgenes?[,;]|\bgenes for\b",tl) and not re.search(r"complete genome|mitochondrial genome|mitochondrion genome|genome assembly, organelle|whole genome",tl): return "gene_cds_title"
    if SEG.search(t) and not re.search(r"genome assembly, organelle",tl): return "segment_title"
    if WH.search(t): return "whole_title"
    if "complete" in tl: return "complete_other"
    return "no_complete_claim"
rows=[];miss=[]
for r in M:
    t=T.get(r["acc"])
    if t is None: miss.append(r["acc"]); t=""
    f=FE.get(r["acc"]) or {}; ncds=len({x[3] for x in f.get("cds",[])}); ntr=len(f.get("trna",[]))
    rows.append(dict(acc=r["acc"],species=r["species"],phylum=r["phylum"],role=r["role"],provenance=r["provenance"],topology=r["topology"],length=int(r["length"]),uaa=r["uaa"],n33=r["n33"],base=r["base"],genes=r["genes"],ann_cds=ncds,ann_trna=ntr,cls=cls(t) if t else "no_title",title=t))
w=csv.DictWriter(open("results/posthoc2/title_audit.tsv","w"),fieldnames=list(rows[0]),delimiter="\t"); w.writeheader(); w.writerows(rows)
o=open("results/posthoc2/title_audit_summary.txt","w"); W=lambda *x:(print(*x),print(*x,file=o))
W("POST-HOC. included records",len(rows),"without a title found:",len(miss),miss[:6])
C=[r for r in rows if r["role"]=="confirmatory" and r["provenance"]!="named_extra"]; S=[r for r in C if r["uaa"]=="Stop" and r["n33"]=="resolved"]; H=[r for r in S if r["base"]!="U"]
W("class\tincluded(32,982)\tconfirmatory\tresolved Stop(28,583)\t17 hits")
for k,_ in collections.Counter(r["cls"] for r in rows).most_common(): W(f"{k}\t{sum(r['cls']==k for r in rows)}\t{sum(r['cls']==k for r in C)}\t{sum(r['cls']==k for r in S)}\t{sum(r['cls']==k for r in H)}")
W("\nby provenance (included):"); [W("  ",k,v) for k,v in collections.Counter((r["provenance"],r["cls"]) for r in rows).most_common()]
for k in ("gene_cds_title","partial_title","segment_title","complete_other","no_complete_claim"):
    x=[r for r in rows if r["cls"]==k]; W(f"\n{k}: {len(x)} examples"); [W("   ",r["acc"],r["topology"],r["length"],r["genes"],r["uaa"],r["n33"],"|",r["title"][:110]) for r in x[:12]]
W("\nprofile genes found (of 13) among resolved Stop:",sorted(collections.Counter(int(r["genes"] or 0) for r in S).items()))
W("linear topology among resolved Stop:",sum(r["topology"]=="linear" for r in S)," linear & <10 profile genes:",sum(r["topology"]=="linear" and int(r["genes"] or 0)<10 for r in S))
W("shorter than 12 kb among resolved Stop:",sum(r["length"]<12000 for r in S))
W("OQ946609.1:",[r for r in rows if r["acc"]=="OQ946609.1"])
W("\nthe 17:"); [W("   ",r["acc"],r["species"],r["topology"],r["length"],"genes",r["genes"],r["cls"],"|",r["title"][:100]) for r in H]
