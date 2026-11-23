"""Run Codetta on the control panel. usage: codetta_run.py list out.jsonl workers db evalue prob max_fraction"""
import sys,os,json,subprocess,tempfile,shutil,csv
from concurrent.futures import ProcessPoolExecutor,as_completed
sys.path.insert(0,"src"); import panel_guard
RES=os.path.abspath("data/pfam/mito"); RES_P=os.path.abspath("env/tools/codetta/resources"); PY=os.path.abspath("env/prefix/bin/python"); COD=os.path.abspath("env/tools/codetta/codetta.py")
TCAG="".join(a+b+c+" " for a in "TCAG" for b in "TCAG" for c in "TCAG").split()
def one(a):
    f,db,ev,pr,mf,mode=a; pid=os.path.basename(f).replace(".fna",""); d=tempfile.mkdtemp()
    fa=os.path.join(d,"q.fna"); shutil.copy(f,fa)
    try:
        # raw strings before topology extension: append first 3 kb for circular replicons via duplicate header? Codetta reads one sequence; keep as is
        r=subprocess.run([PY,"-W","ignore",COD,"q.fna","-p",db,"-m","-e",str(ev),"-r",str(pr),"-f",str(mf),"--align_output","q.al","--inference_output","q.inf","--resource_directory",RES_P if db.startswith("Pfam") else RES],cwd=d,capture_output=True,text=True,timeout=1800)
        line=[l for l in r.stdout.splitlines() if l.startswith("Genetic code:")]
        if not line: return dict(panel_id=pid,status="error",reason=(r.stderr or r.stdout)[-200:])
        code=line[0].split(":",1)[1].strip(); t=TCAG
        return dict(panel_id=pid,status="ok",TAA=code[t.index("TAA")],TAG=code[t.index("TAG")],TGA=code[t.index("TGA")],TAT=code[t.index("TAT")],TAC=code[t.index("TAC")],code=code)
    except Exception as e: return dict(panel_id=pid,status="error",reason=str(e)[:150])
    finally: shutil.rmtree(d,ignore_errors=True)
if __name__=="__main__":
    fs=[l.strip() for l in open(sys.argv[1]) if l.strip()]; panel_guard.check([os.path.basename(f).replace(".fna","") for f in fs]); out=sys.argv[2]; n=int(sys.argv[3]); db,ev,pr,mf=sys.argv[4:8]
    done={json.loads(l)["panel_id"] for l in open(out)} if os.path.exists(out) else set(); todo=[(f,db,ev,pr,mf,0) for f in fs if os.path.basename(f).replace(".fna","") not in done]
    with ProcessPoolExecutor(n) as ex,open(out,"a") as o:
        for fu in as_completed([ex.submit(one,t) for t in todo]): o.write(json.dumps(fu.result())+"\n"); o.flush()
