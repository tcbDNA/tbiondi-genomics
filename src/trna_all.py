"""All mitochondrial tRNA calls (ARAGORN -mt and tRNAscan-SE -O) per genome; panel-guarded. usage: trna_all.py filelist out.jsonl workers"""
import sys,os,json,re,subprocess
from concurrent.futures import ProcessPoolExecutor,as_completed
sys.path.insert(0,"src"); import panel_guard
B="env/prefix/bin/"
def one(f):
    pid=os.path.basename(f).replace(".fna",""); d=dict(panel_id=pid,aragorn=[],trnascan=[],err="")
    try:
        r=subprocess.run([B+"aragorn","-mt","-gc5","-w",f],capture_output=True,text=True,timeout=300)
        for l in r.stdout.splitlines():
            m=re.match(r"\s*\d+\s+mtRNA-(\S+)\s+(c?)\[(-?\d+),(\d+)\]\s+(\d+)\s+\((\w+)\)",l)
            if m: d["aragorn"].append([int(m[3]),int(m[4]),"-" if m[2] else "+",m[1],m[6].upper()])
        r=subprocess.run([B+"tRNAscan-SE","-O","-Q","-q","--brief",f],capture_output=True,text=True,timeout=900,stdin=subprocess.DEVNULL)
        for l in r.stdout.splitlines():
            x=l.split()
            if len(x)>=9: a,b=int(x[2]),int(x[3]); d["trnascan"].append([min(a,b),max(a,b),"+" if a<b else "-",x[4],x[5].upper(),float(x[8])])
    except Exception as e: d["err"]=str(e)[:120]
    return d
if __name__=="__main__":
    fs=[l.strip() for l in open(sys.argv[1]) if l.strip()]; panel_guard.check([os.path.basename(f).replace(".fna","") for f in fs]); out=sys.argv[2]
    done={json.loads(l)["panel_id"] for l in open(out)} if os.path.exists(out) else set(); todo=[f for f in fs if os.path.basename(f).replace(".fna","") not in done]
    with ProcessPoolExecutor(int(sys.argv[3])) as ex,open(out,"a") as o:
        for k,fu in enumerate(as_completed([ex.submit(one,f) for f in todo])):
            o.write(json.dumps(fu.result())+"\n"); o.flush()
            if k%40==0: print(k,len(todo),flush=True)
    print("finished",flush=True)
