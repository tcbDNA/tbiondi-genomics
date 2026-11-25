"""POST-HOC, EXPLORATORY: anticodon-stem pairs counted directly from the 5 nt flanking the called 7-nt loop, versus the covariance-model-column count used by the locked code (stem_pairs). No call or status is changed.
Writes results/gate3/stem_direct.tsv."""
import json,csv,collections,sys
PAIR={("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
M={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included"}
w=csv.writer(open("results/gate3/stem_direct.tsv","w"),delimiter="\t"); w.writerow("acc species phylum class role uaa n33 base cm_stem direct_stem loop arm5 arm3 loop_match".split()); n=0
for l in open("results/gate3/run/out.jsonl"):
    d=json.loads(l); m=M.get(d["acc"]); t=d.get("trny") or {}
    if not m or t.get("status")!="called": continue
    a,b=sorted((int(t["start"]),int(t["end"]))); g=int(t["g34_pos"]); s="".join(x.strip() for x in open("data/sweep/fa/"+str(sum(map(ord,d["acc"]))%200).zfill(3)+f"/{d['acc']}.fna") if x[0]!=">").upper()
    x=s[max(0,a-11):b+10]; off=a-1-max(0,a-11)     # work on a window around the gene
    W=s[a-1:b]; 
    if t["strand"]=="-": W=rc(W); j=b-g
    else: j=g-a
    W=W.replace("T","U")
    pad5=s[max(0,a-6):a-1]; pad3=s[b:b+5]
    full=(rc(s[max(0,a-11):b+10]) if t["strand"]=="-" else s[max(0,a-11):b+10]).replace("T","U")
    lp=10 if a-11>=0 else a-1
    j2=j+(min(10,len(s)-b) if t["strand"]=="-" else lp)
    loop=full[j2-2:j2+5]; s5=full[j2-7:j2-2]; s3=full[j2+5:j2+10]
    ok=loop==t["loop"]
    direct=sum((s5[i],s3[4-i]) in PAIR for i in range(5)) if len(s5)==5 and len(s3)==5 and ok else ""
    w.writerow([d["acc"],m["species"],m["phylum"],m["class"],m["role"],m["uaa"],m["n33"],t["N33"],t["stem_pairs"],direct,t["loop"],s5,s3,ok]); n+=1
print("rows",n)
