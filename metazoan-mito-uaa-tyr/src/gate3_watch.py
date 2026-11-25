"""Interim watcher at the user's request (continuous looks, logged in results/gate3/interim/LOG.md). Uses the committed classify functions on newly finished genomes only.
Emits one line per hit (Tyr, weak, or non-U33 resolved/stem_audit_review in a Stop genome) and when the sweep finishes. State in watch_state.json."""
import sys,json,csv,time,os
sys.path.insert(0,"src"); import gate3_classify as K
OUT="results/gate3/run/out.jsonl"; ST="results/gate3/interim/watch_state.json"
M={r["acc"]:r for r in csv.DictReader(open("results/gate3/sweep_audit.tsv"),delimiter="\t")}
FE={}
for l in open("data/sweep/features_fixed.jsonl"): x=json.loads(l); FE[x["acc"]]=x["f"]
pos=json.load(open(ST))["pos"] if os.path.exists(ST) else 0
if pos==0:
    # start after the look 4 snapshot (same file, earlier lines already reported)
    n=sum(1 for _ in open("results/gate3/interim/look4/out.jsonl"))
    with open(OUT,"rb") as f:
        for _ in range(n): f.readline()
        pos=f.tell()
nt=0
while True:
    with open(OUT,"rb") as f:
        f.seek(pos); lines=f.readlines()
    if lines and not lines[-1].endswith(b"\n"): lines=lines[:-1]
    hits=[]
    for l in lines:
        try: d=json.loads(l)
        except Exception: continue
        m=M.get(d["acc"])
        if not m: continue
        c,note=K.uaa_class(d); z=K.n33(d,FE.get(d["acc"])); u=d.get("uaa") or {}; nt+=1
        tag=m["role"] if m["role"]!="confirmatory" else ""
        if c in("Tyr","weak"): hits.append(f"HIT {c} {d['acc']} {m['species']} {m['phylum']}/{m['class'] or m['order']} {m['family']} {m['provenance']} {tag} TAA {u.get('n_taa_internal')} genes {u.get('n_taa_genes')} tyrcol {u.get('n_taa_tyrcol')} LLR {u.get('llr_Y')} internalTAG {u.get('n_tag_internal')} N33 {z['status']} {z['N33']}")
        elif c=="Stop" and z["status"] in("resolved","stem_audit_review") and z["N33"]!="U": hits.append(f"HIT non-U33 Stop {d['acc']} {m['species']} {m['phylum']}/{m['class'] or m['order']} {m['family']} {m['provenance']} {tag} {z['N33']} {z['loop']} {z['status']} score {z['score']} agree {z['n_agree']}")
        if d["err"]: hits.append(f"MODULE_ERROR {d['acc']} {m['species']} {d['err']}")
    pos+=sum(len(l) for l in lines); json.dump({"pos":pos},open(ST,"w"))
    for h in hits: print(h,flush=True)
    last=open("results/gate3/run/progress.log").read().strip().splitlines()[-1]
    if "FINISHED" in last or "GAVE UP" in last or "stopped: disk floor" in last: print("SWEEP ENDED:",last,"| watcher processed",nt,"new genomes",flush=True); break
    time.sleep(300)
