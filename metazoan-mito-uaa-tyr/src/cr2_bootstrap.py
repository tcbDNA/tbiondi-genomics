"""Cold review 2, check 8 / item 18 (POST-HOC): figure 3 left panel. Genome-clustered bootstrap 95% CI (resample the 65 source genomes, keep all their perturbation replicates) beside the plotted Wilson intervals.
usage: cr2_bootstrap.py -> results/posthoc2/fig3_bootstrap.tsv"""
import json,collections,random,math
def wil(k,n,z=1.96): p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; e=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-e,c+e
R=[json.loads(l) for l in open("results/gate2/planted_panel.jsonl")]; CX=[json.loads(l) for l in open("results/gate2/planted_codetta.jsonl")]
T=lambda u:u["n_taa_internal"]>=5 and u["n_taa_genes"]>=3 and u["n_taa_tyrcol"]>=3 and u["llr_Y"]>=3
print("codetta row keys",list(CX[0])[:14],len(CX))
P=collections.defaultdict(dict)
for r in R: P[r["frac"]].setdefault(r["base"],{})[r["rep"]]=bool(T(r))
cod={}
for c in CX:
    k=(c.get("base") or c.get("panel_id","").replace("PLANTED_",""),c.get("frac"),c.get("rep",0)); cod[k]=c.get("TAA")
B=10000; rng=random.Random(20261125); out=open("results/posthoc2/fig3_bootstrap.tsv","w"); out.write("series\tplanted_frac\tgenomes\tobservations\tcalled\trate\twilson_lo\twilson_hi\tboot_lo\tboot_hi\twidth_ratio\n")
def boot(d):
    g=sorted(d); v=[(sum(x.values()) if isinstance(x,dict) else x[0],len(x) if isinstance(x,dict) else x[1]) for x in (d[i] for i in g)]; n=len(v); s=[]
    for _ in range(B):
        a=b=0
        for _ in range(n): q=v[rng.randrange(n)]; a+=q[0]; b+=q[1]
        s.append(a/b)
    s.sort(); return s[int(.025*B)],s[int(.975*B)-1]
def line(nm,f,d):
    k=sum(sum(x.values()) for x in d.values()); n=sum(len(x) for x in d.values()); lo,hi=wil(k,n); bl,bh=boot(d); out.write("\t".join(map(str,[nm,f,len(d),n,k,round(k/n,4),round(lo,4),round(hi,4),round(bl,4),round(bh,4),round((bh-bl)/(hi-lo),2)]))+"\n"); print(nm,f,f"{k}/{n}={k/n:.3f} Wilson {100*lo:.1f}-{100*hi:.1f} cluster bootstrap {100*bl:.1f}-{100*bh:.1f}")
for f in (0.05,0.1,0.25,0.5,1.0): line("counter",f,P[f])
for f in (0.05,0.1,0.25,0.5,1.0):
    d={}
    for b,x in P[f].items():
        if 0 in x: d[b]={0:bool(x[0] and cod.get((b,f,0))=="Y")}
    line("counter+Codetta(rep 0)",f,d)
# within-genome agreement at 10%
import statistics
x=[sum(v.values()) for v in P[0.1].values()]; print("10% level: replicates called per genome",sorted(collections.Counter(x).items()))
