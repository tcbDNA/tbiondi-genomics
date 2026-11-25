"""Collect MTRF1L-like candidates: (a) proteomes by blastp against the human refs, (b) unannotated assemblies by miniprot with bivalve/nematode refs.
Writes data/mtrf/candidates.faa and results/work/mtrf_candidates.tsv. Control-panel species only (names come from the panel-linked files)."""
import sys,os,subprocess,glob,csv,re,zipfile,io
import sys as _s; _s.path.insert(0,"src"); import mtrf_guard
B="env/prefix/bin/"; REF="data/mtrf/human_refs.faa"; OUT=[]; FA=open("data/mtrf/candidates.faa","w"); seen=set()
def emit(src,sp,pid,s,ev): 
    if (src,pid) in seen: return
    seen.add((src,pid)); FA.write(f">{src}|{pid}\n{s}\n"); OUT.append([src,sp,pid,len(s),ev])
def proteome(src,sp,faa):
    mtrf_guard.check([src])
    db="results/work/db_x"; subprocess.run([B+"makeblastdb","-in",faa,"-dbtype","prot","-out",db],capture_output=True)
    r=subprocess.run([B+"blastp","-query",REF,"-db",db,"-evalue","1e-8","-max_target_seqs","20","-outfmt","6 qseqid sseqid evalue bitscore sseq","-num_threads","4"],capture_output=True,text=True)
    # full sequences
    seqs={}; cur=None
    for l in open(faa):
        if l[0]==">": cur=l[1:].split()[0]; seqs[cur]=""
        else: seqs[cur]+=l.strip()
    best={}
    for l in r.stdout.splitlines():
        q,s,e,b,_=l.split("\t")
        if s not in best or float(b)>best[s][1]: best[s]=(q,float(b),float(e))
    for s,(q,b,e) in sorted(best.items(),key=lambda x:-x[1][1])[:6]:
        seq=seqs.get(s)
        if seq: emit(src,sp,s,seq.rstrip("*"),e)
for f,sp in (("Ach","Anomia chinensis"),("Hsq","Heteranomia squamula"),("Pmu","Plicatula muricata"),("Pvi","Placuna vitream")): proteome("FIG_"+f,sp,f"data/figshare/{f}.chr.pep")
for l in open("results/audit/mtrf_bg_fetch.tsv"):
    pid,sp,acc,st=l.rstrip("\n").split("\t")
    if st=="ok": proteome("BG_"+acc,sp,f"data/mtrf/bg/{acc}.faa")
FA.flush(); FA.close()
csv.writer(open("results/work/mtrf_candidates.tsv","w"),delimiter="\t").writerows([["source","species","protein","length","evalue"]]+OUT); print("candidates",len(OUT))
