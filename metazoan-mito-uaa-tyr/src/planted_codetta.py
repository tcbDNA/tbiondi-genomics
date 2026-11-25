"""Codetta (OXPHOS profiles, -m, Lin thresholds e=1e-5, r=0.7) on planted sequences, replicate 0 only."""
import sys,json,os,glob,re
sys.path.insert(0,"src"); import panel_guard
from codetta_run import one
from concurrent.futures import ProcessPoolExecutor,as_completed
fs=sorted(glob.glob("results/work/planted/*__0.fna")); ok=panel_guard.allowed()
jobs=[]
for f in fs:
    b=os.path.basename(f)[:-4]; pid,frac,rep=b.rsplit("__",2); assert pid in ok, pid
    if float(frac)>0: jobs.append((f,"metazoan_mito_proteins.hmm","1e-5","0.7","0.01",0,pid,float(frac)))
def run(j):
    f,db,e,r,m,mode,pid,frac=j; d=one((f,db,e,r,m,mode)); d.update(base=pid,frac=frac); return d
if __name__=="__main__":
    with ProcessPoolExecutor(3) as ex, open("results/work/planted_codetta.jsonl","w") as o:
        for fu in as_completed([ex.submit(run,j) for j in jobs]): o.write(json.dumps(fu.result())+"\n")
    print(len(jobs))
