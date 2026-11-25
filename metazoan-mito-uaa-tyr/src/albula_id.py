"""Cold review follow-up: which bonefish is it? Compare a query mitochondrial sequence (CYTB or COI region) with all GenBank Albula mitochondrial records (species by GenBank label), nearest records by blastn identity over >=200 bp.
usage: albula_id.py QUERY.fa LABEL"""
import sys,json,subprocess,collections,os
q,label=sys.argv[1:3]; B="env/prefix/bin/"; W="results/work/albula/"; os.makedirs(W,exist_ok=True)
R=json.load(open("results/work/albula_mt.json")); S={r[0]:r for r in R}
open(W+"db.fa","w").write("".join(f">{r[0]}\n{r[5]}\n" for r in R)); subprocess.run([B+"makeblastdb","-in",W+"db.fa","-dbtype","nucl","-out",W+"db"],capture_output=True)
r=subprocess.run([B+"blastn","-task","blastn","-query",q,"-db",W+"db","-evalue","1e-20","-max_target_seqs","1000","-num_threads","4","-outfmt","6 sseqid pident length qstart qend"],capture_output=True,text=True).stdout
best={}
for l in r.splitlines():
    s,p,n,a,b=l.split("\t"); n=int(n)
    if n<200: continue
    if s not in best or n>best[s][1]: best[s]=(float(p),n)
by=collections.defaultdict(list)
for s,(p,n) in best.items(): by[S[s][1]].append((p,n,s,S[s][4],S[s][3][:50]))
print("==",label)
for sp,v in sorted(by.items(),key=lambda x:-max(y[0] for y in x[1]))[:10]: v.sort(reverse=True); print(f"   {sp}: records {len(v)}, best {v[0][0]:.1f}% over {v[0][1]} bp ({v[0][2]}, {v[0][3]}), median {sorted(y[0] for y in v)[len(v)//2]:.1f}%")
