"""Cold review follow-up 1: reads of SRR13036352 (results/work/albula/hi2016.reads.fa, from albula_reads_id.py) carrying CYTB 25-mers diagnostic for A. glossodonta (CM165726) or A. virgata (KJ910045), over the 502 bp window the A. virgata record covers."""
import json,subprocess
W="results/work/albula/"; K=25; rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
R={r[0]:r for r in json.load(open("results/work/albula_mt.json"))}
km=lambda s:{s[i:i+K] for i in range(len(s)-K+1) if set(s[i:i+K])<=set("ACGT")}
gl="".join(l.strip() for l in open(W+"CM165726_CYTB.fa") if l[0]!=">"); vg=R["KJ910045.1"][5]; open(W+"v.fa","w").write(">v\n"+vg+"\n")
al=subprocess.run(["env/prefix/bin/blastn","-task","blastn","-query",W+"v.fa","-subject",W+"CM165726_CYTB.fa","-outfmt","6 sstart send"],capture_output=True,text=True).stdout.split(); a,b=sorted(map(int,al[:2]))
G,V=km(gl),km(vg); GW=km(gl[a-1:b])-V; VW=V-G; ng=nv=nb=0
for l in open(W+"hi2016.reads.fa"):
    if l[0]==">": continue
    s=l.strip().upper(); k=km(s)|km(rc(s)); g,v=bool(k&GW),bool(k&VW); ng+=g; nv+=v; nb+=g and v
print("diagnostic 25-mers: glossodonta",len(GW),"virgata",len(VW)); print("reads with glossodonta-diagnostic",ng,"virgata-diagnostic",nv,"both",nb)
