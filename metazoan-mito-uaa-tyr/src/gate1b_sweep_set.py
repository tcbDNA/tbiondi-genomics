"""Gate 1 extension: one INSDC-only record per species absent from RefSeq+organelle set. Provenance recorded. Metadata only (no caller)."""
import csv,collections,re
TAX={}
for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
    f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; TAX[f[0]]=f
OPEN=re.compile(r"\b(sp|cf|aff|nr|gen)\.|\bsp\b|uncultured|environmental|hybrid| x ")
S=list(csv.DictReader(open("results/gate1_species.tsv"),delimiter="\t")); have={r["species"] for r in S}
rows=[]
for fn in ("results/work/insdc_invert.tsv","results/work/insdc_chordate.tsv","results/work/insdc_invert_gatitle_shorttopup.tsv"):
    for r in csv.DictReader(open(fn),delimiter="\t"): rows.append(r)
print("INSDC rows",len(rows),"unique accessions",len({r["accession"] for r in rows}))
aud=csv.writer(open("results/audit/gate1b_insdc.tsv","w"),delimiter="\t"); aud.writerow(["accession","taxid","decision","reason"])
best={}
for r in rows:
    t=(r["title"] or ""); f=TAX.get(r["taxid"]); 
    if not f: aud.writerow([r["accession"],r["taxid"],"drop","no_taxonomy"]); continue
    sp=f[2] or f[1]; L=int(r["length"]); ph=f[7] or "(no phylum)"
    if sp in have and not OPEN.search(sp): aud.writerow([r["accession"],r["taxid"],"drop","species_in_RefSeq_or_organelle_set"]); continue
    q=0
    if "UNVERIFIED" in t: q-=3
    if "complete" in t.lower(): q+=2
    if "circular" in (r["topology"] or ""): q+=1
    if "partial" in t.lower() or "nearly complete" in t.lower(): q-=1
    key=sp+("|"+r["accession"] if OPEN.search(sp) else "")
    cand=(q,L,r,ph,sp)
    if key not in best or cand[:2]>best[key][:2]: best[key]=cand
n=collections.Counter()
w=csv.writer(open("results/insdc_only_set.tsv","w"),delimiter="\t"); w.writerow("accession taxid species phylum length topology quality provenance title".split())
for k,(q,L,r,ph,sp) in best.items(): w.writerow([r["accession"],r["taxid"],sp,ph,L,r["topology"],q,"INSDC_only",r["title"]]); n[ph]+=1
print("INSDC-only species added:",len(best),n.most_common(12)); print("quality counts",collections.Counter(x[0] for x in best.values()))
