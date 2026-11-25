"""Cold review 2, check 2 part 2 (POST-HOC): the C. infrequens copy on OZ318123.1 (whole mitogenome, 99.8-99.9%) carries 20 of the N33 reads although nuclear depth is ~1.8x. Depth along it and beyond its edges, reads across the edges, and combined depth per mitogenome position. appends to results/posthoc2/craniella_depth.txt"""
import subprocess,re,collections,statistics as st
M="env/tools/mapenv/bin/samtools"; bam="results/work/r2/cran/infr_hifi.bam"; B="env/prefix/bin/"; o=open("results/posthoc2/craniella_depth2.txt","w"); P=lambda *x:(print(*x),print(*x,file=o)); run=lambda *a: subprocess.run([M]+list(a),capture_output=True,text=True).stdout
def depth(c,a,b): return {int(x.split("\t")[1]):int(x.split("\t")[2]) for x in run("depth","-aa","-G","2304","-r",f"{c}:{a}-{b}",bam).splitlines()}
P("POST-HOC check 2, part 2")
r=subprocess.run([B+"blastn","-task","blastn","-query","results/work/r2/cran/OZ318138.1.q","-db","results/work/numt/GCA_974939545.1","-evalue","1e-20","-dust","no","-max_target_seqs","2000","-num_threads","4","-outfmt","6 sseqid pident length qstart qend sstart send"],capture_output=True,text=True).stdout
H=[(f[0],float(f[1]),int(f[2]),*map(int,f[3:7])) for f in (l.split("\t") for l in r.splitlines())]
big=[h for h in H if h[0]!="OZ318138.1" and h[2]>=3000]
for h in sorted(big,key=lambda h:-h[2])[:14]:
    a,b=sorted(h[5:7]); d=depth(h[0],max(1,a-3000),b+3000); ins=[d[x] for x in range(a,b+1) if x in d]; l=[d[x] for x in range(max(1,a-3000),a-200) if x in d]; q=[d[x] for x in range(b+200,b+3000) if x in d]
    span=[0,0]
    for i,e in enumerate((a,b)):
        for f in (x.split("\t") for x in run("view","-F","2304",bam,f"{h[0]}:{e}-{e}").splitlines()):
            s=int(f[3]); z=s+sum(int(n) for n,k in re.findall(r"(\d+)([MDN=X])",f[5]))-1; span[i]+=(s<=e-500 and z>=e+500)
    P("  ",h[0],f"{a}-{b}",f"id {h[1]} len {h[2]} mt {h[3]}-{h[4]}","median depth inside",st.median(ins),"| outside left",st.median(l) if l else None,"right",st.median(q) if q else None,"| reads crossing the edges by >=500 nt each way: left",span[0],"right",span[1])
# the OZ318123 block as a whole
a,b=6113449,6131541; d=depth("OZ318123.1",a-30000,b+30000)
P("  OZ318123.1 block 6,113,449-6,131,541: median inside",st.median([d[x] for x in range(a,b+1)]),"| 5-30 kb left",st.median([d[x] for x in range(a-30000,a-5000)]),"| 5-30 kb right",st.median([d[x] for x in range(b+5000,b+30000)]))
# combined depth by mitogenome position: record + OZ318123 copy
dm=depth("OZ318138.1",1,18285); pos={}
for c,idn,ln,qs,qe,ss,se in H:
    if c!="OZ318123.1" or ln<8000: continue
    step=1 if se>ss else -1
    for k in range(qs,qe+1): pos[k]=ss+step*(k-qs)
for nm,u,v in (("COX1",11592,13115),("ND1",13615,14544),("trnY",14818,14890),("ND2",15261,16178),("ND5",16352,18154),("ATP6/COX region 5000-6000",5000,6000),("whole",1,18285)):
    x=[dm.get(k,0) for k in range(u,v+1)]; y=[d.get(pos[k],0) if k in pos else 0 for k in range(u,v+1)]; P(f"   mitogenome {nm} {u}-{v}: record {st.median(x)}, OZ318123 copy {st.median(y)}, sum {st.median([i+j for i,j in zip(x,y)])}")
# are there fixed differences between record and copy, and do reads split cleanly?
P("  circularity: reads whose alignment ends at the record's end (<=50 nt) and have a supplementary part at its start:")
n=0
for f in (l.split("\t") for l in run("view","-F","2304",bam,"OZ318138.1").splitlines()):
    if any(x.startswith("SA:Z:OZ318138.1") for x in f[11:]): n+=1
P("   reads split within OZ318138.1 (origin-crossing):",n)
