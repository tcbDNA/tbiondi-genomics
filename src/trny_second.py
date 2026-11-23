"""Independent second callers for trnY: ARAGORN -mt and tRNAscan-SE -O. Emits GUA-anticodon positions (G34, 1-based, strand). One JSON line per genome."""
import sys,os,json,re,subprocess
from concurrent.futures import ProcessPoolExecutor,as_completed
sys.path.insert(0,"src"); import panel_guard
B="env/prefix/bin/"
def aragorn(fa,L):
    out=[]; r=subprocess.run([B+"aragorn","-mt","-gc5","-w",fa],capture_output=True,text=True,timeout=300)
    for l in r.stdout.splitlines():
        m=re.match(r"\s*\d+\s+mtRNA-(\S+)\s+(c?)\[(-?\d+),(\d+)\]\s+(\d+)\s+\((\w+)\)",l)
        if m and m.group(6).lower()=="gta": 
            a,b,k=int(m.group(3)),int(m.group(4)),int(m.group(5)); pos=(a+k-1) if not m.group(2) else (b-k+1); out.append((pos,"-" if m.group(2) else "+",m.group(1)))
    return out
def trnascan(fa):
    out=[]; r=subprocess.run([B+"tRNAscan-SE","-O","-q","--detail","-o","/dev/stdout","-f","/dev/null",fa],capture_output=True,text=True,timeout=900)
    return r
def scan2(fa):
    r=subprocess.run([B+"tRNAscan-SE","-O","-Q","-q","--brief",fa],capture_output=True,text=True,timeout=900,stdin=subprocess.DEVNULL); out=[]
    for l in r.stdout.splitlines():
        f=l.split()
        if len(f)>=9 and f[5].upper()=="GTA": out.append((int(f[2]),int(f[3]),float(f[8])))
    return out
def one(f):
    pid=os.path.basename(f).replace(".fna",""); d=dict(panel_id=pid)
    try: d["aragorn"]=aragorn(f,0)
    except Exception as e: d["aragorn"]=[]; d["aragorn_err"]=str(e)[:100]
    try: d["trnascan"]=scan2(f)
    except Exception as e: d["trnascan"]=[]; d["trnascan_err"]=str(e)[:100]
    return d
if __name__=="__main__":
    fs=[l.strip() for l in open(sys.argv[1]) if l.strip()]; panel_guard.check([os.path.basename(f).replace(".fna","") for f in fs]); out=sys.argv[2]
    done={json.loads(l)["panel_id"] for l in open(out)} if os.path.exists(out) else set(); todo=[f for f in fs if os.path.basename(f).replace(".fna","") not in done]
    with ProcessPoolExecutor(int(sys.argv[3])) as ex,open(out,"a") as o:
        for fu in as_completed([ex.submit(one,f) for f in todo]): o.write(json.dumps(fu.result())+"\n"); o.flush()
