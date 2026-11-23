"""Gate 3 uniform pipeline (no pre-filter; Codetta is run on every genome). Same functions as Gate 2; U.analyse is memoised per genome (pure function).
usage: gate3_run.py manifest.tsv outdir workers   manifest columns: acc path topology
Checkpointed: one JSON line per genome in outdir/out.jsonl; rerun resumes. Progress in outdir/progress.log. Module failures are recorded per module, never dropped."""
import sys,os,json,csv,re,subprocess,time,shutil
sys.path.insert(0,"src"); import panel_guard
B="env/prefix/bin/"
def trnas(fa):
    ar=[];ts=[];sec=dict(aragorn=[],trnascan=[])
    r=subprocess.run([B+"aragorn","-mt","-gc5","-w",fa],capture_output=True,text=True,timeout=600)
    for l in r.stdout.splitlines():
        m=re.match(r"\s*\d+\s+mtRNA-(\S+)\s+(c?)\[(-?\d+),(\d+)\]\s+(\d+)\s+\((\w+)\)",l)
        if m:
            ar.append([int(m[3]),int(m[4]),"-" if m[2] else "+",m[1],m[6].upper()])
            if m[6].lower()=="gta": a,b,k=int(m[3]),int(m[4]),int(m[5]); sec["aragorn"].append([(a+k-1) if not m[2] else (b-k+1),"-" if m[2] else "+",m[1]])
    r=subprocess.run([B+"tRNAscan-SE","-O","-Q","-q","--brief",fa],capture_output=True,text=True,timeout=1800,stdin=subprocess.DEVNULL)
    for l in r.stdout.splitlines():
        x=l.split()
        if len(x)>=9:
            a,b=int(x[2]),int(x[3]); ts.append([min(a,b),max(a,b),"+" if a<b else "-",x[4],x[5].upper(),float(x[8])])
            if x[5].upper()=="GTA": sec["trnascan"].append([a,b,float(x[8])])
    return dict(aragorn=ar,trnascan=ts),sec
def process(j):
    acc,path,topo=j; import trny_extract as T, uaa_caller as U, aux_checks as AX, codetta_run as C
    if not hasattr(U,"_memo"):
        U._memo={}; U._raw=U.analyse
        def an(fa,circular=True,E=1e-5,tyr_thr=0.5):
            k=(fa,circular,E)
            if k not in U._memo: U._memo.clear(); U._memo[k]=U._raw(fa,circular,E)
            return U._memo[k]
        U.analyse=an
    d=dict(acc=acc,err={}); circ=topo!="linear"; t0=time.time()
    def mod(n,f):
        try: d[n]=f()
        except Exception as e: d["err"][n]=f"{type(e).__name__}: {e}"[:160]; d[n]=None
    mod("trny",lambda:T.call(path)); mod("t",lambda:trnas(path)); tr,sec=d.pop("t") or ({},{}); d["trna"]=tr; d["second"]=sec
    mod("uaa",lambda:U.summarise(path,circular=circ)); mod("first",lambda:AX.terminal(path,circ,"first",tr)); mod("tagfirst",lambda:AX.terminal(path,circ,"tagfirst",tr)); mod("oe",lambda:AX.codon_oe(path,circ))
    mod("pcg",lambda:[[k["ns"]+1,k["ne"],k["strand"],k["gene"],k["frame"]] for k in U.analyse(path,circ)[2]])
    mod("codetta",lambda:C.one((path,"metazoan_mito_proteins.hmm","1e-5","0.7","0.01",0)))
    d["sec"]=round(time.time()-t0,1); return d
if __name__=="__main__":
    from concurrent.futures import ProcessPoolExecutor,as_completed
    man,od,n=sys.argv[1],sys.argv[2],int(sys.argv[3]); os.makedirs(od,exist_ok=True); J=[(r["acc"],r["path"],r["topology"]) for r in csv.DictReader(open(man),delimiter="\t")]; panel_guard.check([j[0] for j in J])
    out=od+"/out.jsonl"; done=set()
    if os.path.exists(out):
        good=[]
        for l in open(out):
            try: done.add(json.loads(l)["acc"]); good.append(l)
            except Exception: pass          # a torn last line after a crash is redone
        open(out+".tmp","w").writelines(good); os.replace(out+".tmp",out)
    todo=[j for j in J if j[0] not in done]; t0=time.time(); P=open(od+"/progress.log","a",buffering=1); P.write(f"{time.strftime('%F %T')} START total {len(J)} done {len(done)} todo {len(todo)} workers {n}\n")
    it=iter(todo); k=0; nerr=0
    with ProcessPoolExecutor(n,max_tasks_per_child=200) as ex, open(out,"a",buffering=1) as o:
        fut=set()
        def fill():
            while len(fut)<n*4:
                try: fut.add(ex.submit(process,next(it)))
                except StopIteration: return
        fill()
        while fut:
            f=next(as_completed(fut)); fut.discard(f)
            d=f.result(); o.write(json.dumps(d)+"\n"); k+=1; nerr+=bool(d["err"])
            if k%100==0 or k==len(todo):
                r=k/(time.time()-t0); fr=shutil.disk_usage(".").free/1e9; P.write(f"{time.strftime('%F %T')} {len(done)+k}/{len(J)} rate {r*60:.1f}/min ETA {(len(todo)-k)/r/3600:.1f} h module_errors {nerr} free {fr:.0f} GB\n")
                if fr<31: P.write("STOP: free disk near 30 GB\n"); ex.shutdown(cancel_futures=True); sys.exit(2)
            fill()
    P.write(f"{time.strftime('%F %T')} FINISHED {len(done)+k}/{len(J)} module_errors {nerr}\n")
