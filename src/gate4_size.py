"""Estimate the nuclear data needed for Gate 4 (targeted mtRF1a): flagged species with assemblies + up to 2 nearest Stop sisters with assemblies per flagged family. Sizes from NCBI datasets (total_sequence_length)."""
import csv,collections,subprocess,json,random
R=list(csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")); inc=[r for r in R if r["status"]=="included"]; C=[r for r in inc if r["provenance"]!="named_extra"]
A=collections.defaultdict(list)
for a in csv.DictReader(open("results/assemblies_metazoa.tsv"),delimiter="\t"): A[a["organism"]].append(a)
best=lambda sp: sorted(A.get(sp,[]),key=lambda a:(not a["accession"].startswith("GCF"),{"Chromosome":0,"Complete Genome":0,"Scaffold":1,"Contig":2}.get(a["level"],3),-float(a["contig_n50"] or 0)))[:1]
flag=[r for r in C if r["uaa"]=="Tyr" or (r["uaa"]=="Stop" and r["n33"]=="resolved" and r["base"]!="U" and r["role"]=="confirmatory")]
need={}
for r in flag:
    for a in best(r["species"]): need[a["accession"]]=(r["species"],"flagged",r["family"])
import re
for r in flag:
    fam=r["family"] or r["genus"]
    pool=[x for x in C if x["uaa"]=="Stop" and x["acc"]!=r["acc"] and best(x["species"])]
    for rank in("family","order"):
        s=[x for x in pool if r.get(rank) and x.get(rank)==r[rank] and x["species"]!=r["species"]]; random.Random(1).shuffle(s)
        if len(s)>=2 or rank=="order":
            for x in s[:2]:
                for a in best(x["species"]): need.setdefault(a["accession"],(x["species"],f"sister ({rank})",x["family"]))
            if len(s)>=2: break
acc=list(need); size={}
for i in range(0,len(acc),40):
    out=subprocess.run(["env/prefix/bin/datasets","summary","genome","accession",*acc[i:i+40]],capture_output=True,text=True).stdout
    for x in json.loads(out).get("reports",[]): size[x["accession"]]=int(x["assembly_stats"]["total_sequence_length"])
tot=sum(size.values()); print("assemblies needed",len(acc),"with size",len(size),"total bp %.1f Gb"%(tot/1e9),"=> compressed download ~%.1f GB, unpacked ~%.1f GB, peak with indices ~%.1f GB"%(tot*0.28/1e9,tot*1.02/1e9,tot*1.3/1e9))
for a in sorted(acc,key=lambda a:-size.get(a,0)): print(f"  {a} {need[a][0]} [{need[a][1]}] {size.get(a,0)/1e9:.2f} Gb")
