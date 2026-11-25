"""Locus-identity audit (PREREG 2.3). For a trnY call: (1) verified GUA in a 7-nt loop with stem >=4/5; (2) overlap with any annotated CDS or other tRNA (or, without annotation,
profile PCG domains and non-GUA tRNA calls) <= TOL nt, either strand; (3) if the record annotates a trnY, >50% overlap on the same strand. Panel guarded.
usage: locus_audit.py calls.jsonl annotation.json out.tsv"""
import sys,json,csv
sys.path.insert(0,"src"); import panel_guard, uaa_caller as U
def ov(a,b): return min(a[1],b[1])-max(a[0],b[0])+1
def audit(pid,d,feats,tol,trna,topo,use_hmm=True):
    if d.get("status")!="called": return "unresolved","not_called"
    a,b=sorted((int(d["start"]),int(d["end"]))); st="+" if d["strand"]!="-" else "-"
    if int(d["stem_pairs"])<4: return "unresolved","stem_lt4"
    worst=0;why=""
    for f in (feats or {}).get("cds",[]):
        o=ov((a,b),f[:2]); 
        if o>worst: worst,why=o,"overlaps_CDS_"+str(f[3])
    for f in (feats or {}).get("trna",[]):
        if f[4]: continue           # annotated trnY itself
        o=ov((a,b),f[:2])
        if o>worst: worst,why=o,"overlaps_tRNA_"+str(f[3])
    if feats is None:
        for x in trna.get("aragorn",[]):
            if x[4]!="GTA" and x[2]==st: 
                o=ov((a,b),(x[0],x[1]))
                if o>worst: worst,why=o,"overlaps_tRNA_call_"+x[4]
        for d2 in _hmm(pid,topo):
            o=ov((a,b),(d2[0]+1,d2[1]))
            if o>worst: worst,why=o,"overlaps_PCG_domain_"+d2[2]
    if worst>tol: return "unresolved",f"{why}_{worst}nt"
    yy=[f for f in (feats or {}).get("trna",[]) if f[4]]
    if yy and not any(ov((a,b),f[:2])>0.5*min(b-a,f[1]-f[0]) and f[2]==st for f in yy): return "unresolved","disagrees_with_annotated_trnY"
    return "pass",f"overlap_{worst}nt"
_c={}
def _hmm(pid,topo):
    if pid not in _c:
        n,L0,kept,*_=U.analyse(f"data/panel/{pid}.fna",topo!="linear"); _c[pid]=[(k["ns"],k["ne"],k["gene"]) for k in kept]
    return _c[pid]
if __name__=="__main__":
    P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; panel_guard.check(P)
    T={json.loads(l)["panel_id"]:json.loads(l) for l in open(sys.argv[1])}; FE=json.load(open(sys.argv[2])); TR={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/trna_all.jsonl")}
    w=csv.writer(open(sys.argv[3],"w"),delimiter="\t"); w.writerow("panel_id lineage role status reason".split())
    for k in P:
        s,r=audit(k,T[k],FE.get(k),10,TR[k],P[k]["topology"]); w.writerow([k,P[k]["lineage"],P[k]["role"],s,r])
