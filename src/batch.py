"""Resumable parallel runner. usage: batch.py {trny|uaa} filelist out.jsonl workers
Every input gets exactly one JSON line; exceptions are recorded with a reason, never dropped."""
import sys,os,json,csv,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
sys.path.insert(0,"src"); import panel_guard
def one(a):
    kind,f,circ=a; pid=os.path.basename(f).replace(".fna","")
    try:
        if kind=="trny":
            import trny_extract as T; r=T.call(f)
        else:
            import uaa_caller as U; r=U.summarise(f,circular=circ); r["status"]="ok"
    except Exception as e: r=dict(status="exception",reason=f"{type(e).__name__}: {e}"[:200])
    r["panel_id"]=pid; return r
if __name__=="__main__":
    kind,fl,out,n=sys.argv[1:5]; fs=[l.strip() for l in open(fl) if l.strip()]
    ids=[os.path.basename(f).replace(".fna","") for f in fs]; panel_guard.check(ids)
    topo={r["panel_id"]:r["topology"] for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
    done=set()
    if os.path.exists(out): done={json.loads(l)["panel_id"] for l in open(out)}
    todo=[(kind,f,topo.get(i.split("|")[0],"circular")!="linear") for f,i in zip(fs,ids) if i not in done]
    print(kind,"todo",len(todo),"done",len(done),flush=True)
    with ProcessPoolExecutor(int(n)) as ex, open(out,"a") as o:
        for k,fu in enumerate(as_completed([ex.submit(one,t) for t in todo])):
            o.write(json.dumps(fu.result())+"\n"); o.flush()
            if k%50==0: print(k,flush=True)
    print("finished",flush=True)
