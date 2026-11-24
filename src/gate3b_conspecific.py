"""Post-hoc: conspecific/congeneric trnY calls for the remaining candidates (score gate as in conspecific_filter.py: >=20, or >=10 with stem>=4). Records >=1 kb. Lists same-sequence duplicates (RefSeq copies) by sequence hash."""
import sys,re,collections,os,hashlib; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import ncbi, trny_extract as T
SETS=[("Serrivomer",["beanii","sector","lanceolatoides"]),("Stemonidium",["hypomelas"]),("Coloconger",["meadi"]),("Ophichthus",["asakusae"]),("Xyrias",["revulsus"]),("Lycenchelys",["verrillii"]),("Globodera",["ellingtonae"]),("Toxocara",["apodemi"]),("Antigone",["rubicunda"]),("Grus",["rubicunda"]),("Craniella",["zetlandica","infrequens"])]
for gen,sps in SETS:
    try:
        n,ids,_=ncbi.search("nuccore",f'"{gen}"[Organism] AND mitochondrion[filter] AND 1000:40000[SLEN]',300); rows=[]
        for i in range(0,len(ids),50):
            t=ncbi.call("efetch",db="nuccore",id=",".join(ids[i:i+50]),rettype="fasta",retmode="text").decode(errors="replace")
            for blk in t.split(">")[1:]:
                h,_,s=blk.partition("\n"); a=h.split()[0]; s=re.sub(r"\s","",s).upper(); p=f"results/work/cons/{a}.fna"
                if not os.path.exists(p): open(p,"w").write(f">{a}\n{s}\n")
                d=T.call(p)
                if d.get("status")=="called" and (d["score"]>=20 or (d["score"]>=10 and d["stem_pairs"]>=4)): rows.append((a,h.split(" ",1)[-1][:70],len(s),d["N33"],d["loop"],d["score"],d["stem_pairs"],hashlib.md5(s.encode()).hexdigest()[:8]))
        print(f"== {gen}: records >=1kb {len(ids)}, with confident trnY {len(rows)}, N33 {dict(collections.Counter(r[3] for r in rows))}",flush=True)
        for r in rows:
            if any(sp in r[1] for sp in sps): print("   target:",*r,flush=True)
            elif r[3]!="U": print("   other non-U:",*r,flush=True)
    except Exception as e: print("==",gen,"ERROR",str(e)[:150],flush=True)
print("DONE",flush=True)
