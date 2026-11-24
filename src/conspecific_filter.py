"""Re-count the conspecific/congeneric trnY calls with the score gate (>=20, or >=10 with stem>=4 and 7-nt GUA loop as already enforced). Uses cached per-record files; re-runs T.call only for records with cached fasta."""
import sys,re,collections,glob,os; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import ncbi, trny_extract as T
for gen,sp in (("Microdevario","kubotai"),("Albula","vulpes"),("Bufotes","pewzowi")):
    n,ids,_=ncbi.search("nuccore",f'"{gen}"[Organism] AND mitochondrion[filter] AND 500:30000[SLEN]',400); rows=[]
    for i in range(0,len(ids),50):
        t=ncbi.call("efetch",db="nuccore",id=",".join(ids[i:i+50]),rettype="fasta",retmode="text").decode(errors="replace")
        for blk in t.split(">")[1:]:
            h,_,s=blk.partition("\n"); a=h.split()[0]; s=re.sub(r"\s","",s).upper()
            p=f"results/work/cons/{a}.fna"
            if not os.path.exists(p): open(p,"w").write(f">{a}\n{s}\n")
            if len(s)<1000: continue
            d=T.call(p)
            if d.get("status")=="called" and (d["score"]>=20 or (d["score"]>=10 and d["stem_pairs"]>=4)): rows.append((a,h.split(" ",1)[-1],len(s),d["N33"],d["loop"],d["score"],d["stem_pairs"]))
    sp_rows=[r for r in rows if sp in r[1]]
    print(gen,"records",len(ids),"| confident calls",len(rows),collections.Counter(r[3] for r in rows),f"| {sp}:",[(r[0],r[2],r[3],r[4],r[5],r[6]) for r in sp_rows])
    for r in rows:
        if r[3]!="U" and sp not in r[1]: print("   other non-U:",r[0],r[1][:60],r[2],r[3],r[4],r[5])
