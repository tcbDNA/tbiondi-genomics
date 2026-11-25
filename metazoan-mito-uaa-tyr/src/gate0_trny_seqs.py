"""Gate 0: print tRNA-Tyr gene sequence + anticodon loop for UNANNOTATED CONTROL accessions using tRNAscan-SE -O coordinates (controls only)."""
import subprocess, re
CTRL=["PP708078.1","PP708091.1","OZ181722.1","OZ576062.1","MF374802.1","FN908482.1"]
rc=lambda s:s.translate(str.maketrans("ACGT","TGCA"))[::-1]
for a in CTRL:
    seq="".join(l.strip() for l in open(f"data/ncbi/controls/{a}.fna") if l[0]!=">")
    r=subprocess.run(["tRNAscan-SE","-O","-q","--brief",f"data/ncbi/controls/{a}.fna"],capture_output=True,text=True)
    for l in r.stdout.splitlines():
        f=l.split()
        if len(f)<9 or f[4]!="Tyr": continue
        b,e=int(f[2]),int(f[3]); t=seq[b-1:e] if b<e else rc(seq[e-1:b])
        for m in re.finditer("(?=GTA)",t):
            i=m.start(); st5,lp,st3=t[i-7:i-2],t[i-2:i+5],t[i+5:i+10]
            pair=sum(x+y in("AT","TA","GC","CG","GT","TG") for x,y in zip(st5,rc(st3)[::-1].translate(str.maketrans("ACGT","TGCA"))[::-1])) if len(st5)==5 and len(st3)==5 else -1
            p=sum((x+y) in("AT","TA","GC","CG","GT","TG") for x,y in zip(st5,st3[::-1])) if len(st5)==5 and len(st3)==5 else -1
            print(a,f"{b}-{e}",f"score={f[8]}",len(t),"local_idx_N33=",i,"N33=",t[i-1],"stem5",st5,"loop",lp,"stem3",st3,"pairs",p)
        print("   ",t)
