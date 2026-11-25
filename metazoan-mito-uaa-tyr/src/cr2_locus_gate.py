"""Cold review 2, check 1 (POST-HOC, 2026-11-25; locked outputs are read, never written). Uses the locked functions gate3_classify.n33 / trny_decide.agree.
 part A: which second caller corroborates, and on which strand, for the 17 screen positives and the 20 panel positives.
 part B: sensitivity (a) tRNAscan-SE corroboration must be on the same strand; (b) if the top CM arm is not 'resolved', retry the next plausible arm (GUA, 7-nt loop, stem>=3) of the same cmsearch pass, best score first; (ab) both;
         (c) extra, for review item 6: unannotated records also use tRNAscan-SE non-GUA calls (both strands) and ARAGORN non-GUA calls on both strands in the overlap check; annotated records get the same calls added.
usage: cr2_locus_gate.py  -> results/posthoc2/locus_gate_*.tsv, locus_gate_summary.txt"""
import sys,os,json,csv,collections,copy
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"
import trny_decide as TD, gate3_classify as GC, trny_extract as TE, locus_audit as LA
from concurrent.futures import ProcessPoolExecutor
O="results/posthoc2/"; os.makedirs(O,exist_ok=True)
locked_agree=TD.agree
def strands(d,s):
    """detail of corroboration: ARAGORN GUA within 2 nt of G34 by strand; tRNAscan GUA calls with >50% overlap by strand"""
    a,b=sorted((int(d["start"]),int(d["end"]))); st="+" if d["strand"]!="-" else "-"; out=[]
    for x in s.get("aragorn",[]):
        if abs(x[0]-d["g34_pos"])<=2: out.append(("ARAGORN",x[1],"same" if x[1]==st else "OPPOSITE","G34 %+d"%(x[0]-d["g34_pos"])))
    for p,q,sc in s.get("trnascan",[]):
        if min(b,max(p,q))-max(a,min(p,q))>0.5*(b-a): t="+" if p<q else "-"; out.append(("tRNAscan-SE",t,"same" if t==st else "OPPOSITE","score %.1f"%sc))
    return out
def agree_same(d,s):
    if d["status"]!="called" or not d.get("g34_pos"): return (False,False)
    st="+" if d["strand"]!="-" else "-"; s2=dict(s); s2["trnascan"]=[x for x in s.get("trnascan",[]) if ("+" if x[0]<x[1] else "-")==st]
    return locked_agree(d,s2)
def with_agree(f,fn):
    TD.agree=f
    try: return fn()
    finally: TD.agree=locked_agree
def arms(j):
    """all plausible arms of the pass that produced the locked call; returns list sorted by score"""
    acc,path,reason=j; hits,err=TE.run(path,reason=="sens")
    if err: return acc,None,str(err)
    ok=[h for h in hits if h["anticodon"]=="GUA" and h["loop_ok"] and h["stem_pairs"]>=3]; ok.sort(key=lambda h:-h["score"]); return acc,ok,""
def key(z): return (z["status"],z["N33"] if z["status"]!="unresolved" else "")
if __name__=="__main__":
    out=open(O+"locus_gate_summary.txt","w"); W=lambda *x:(print(*x,flush=True),print(*x,file=out,flush=True))
    W("POST-HOC (cold review 2, check 1). Locked results unchanged.")
    M={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")}
    FE={}
    for l in open("data/sweep/features_fixed.jsonl"): x=json.loads(l); FE[x["acc"]]=x["f"]
    R={}
    for l in open("results/gate3/run/out.jsonl"):
        x=json.loads(l); R[x["acc"]]=x
    PATH={r["acc"]:r["path"] for r in csv.DictReader(open("results/gate3/manifest.tsv"),delimiter="\t")}
    C=[a for a,m in M.items() if m["status"]=="included" and m["role"]=="confirmatory" and m["provenance"]!="named_extra"]
    hits=[a for a in C if M[a]["uaa"]=="Stop" and M[a]["n33"]=="resolved" and M[a]["base"]!="U"]
    base={a:GC.n33(R[a],FE.get(a)) for a in C}
    bad=[a for a in C if (base[a]["status"],base[a]["N33"])!=(M[a]["n33"],M[a]["base"])]
    W("confirmatory genomes",len(C),"| locked statuses reproduced from out.jsonl; mismatches:",len(bad),bad[:5])
    # part A
    A=csv.writer(open(O+"locus_gate_corroboration.tsv","w"),delimiter="\t"); A.writerow("set acc species strand score stem N33 ARAGORN_same ARAGORN_opposite tRNAscan_same tRNAscan_opposite detail locked_status".split())
    def rowA(st,acc,sp,d,s,status):
        x=strands(d,s) if d.get("status")=="called" else []; n=lambda c,k:sum(1 for y in x if y[0]==c and y[2]==k)
        A.writerow([st,acc,sp,d.get("strand",""),d.get("score",""),d.get("stem_pairs",""),d.get("N33",""),n("ARAGORN","same"),n("ARAGORN","OPPOSITE"),n("tRNAscan-SE","same"),n("tRNAscan-SE","OPPOSITE"),"; ".join(" ".join(map(str,y)) for y in x),status])
    for a in sorted(hits,key=lambda a:M[a]["species"]): rowA("screen17",a,M[a]["species"],R[a]["trny"],R[a]["second"],M[a]["n33"])
    P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; PT={}; 
    for l in open("results/gate2/trny_panel_final.jsonl"): x=json.loads(l); PT[x["panel_id"]]=x
    PS={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/trny_second.jsonl")}; PF=json.load(open("results/gate2/panel_features_fixed.json")); TR={json.loads(l)["panel_id"]:json.loads(l) for l in open("results/gate2/trna_all.jsonl")}
    pos=[k for k in P if P[k]["role"]=="pos"]
    def panel_status(k,d):
        g=TD.decide(d,PS.get(k,{})); s,r=LA.audit(k,d,PF.get(k),10,TR[k],P[k]["topology"]); return ("pass" if g["final"]=="called" and s=="pass" else "fail",g["final"],g["reason"],s,r,d.get("N33","") if g["final"]=="called" else "")
    pbase={k:panel_status(k,PT[k]) for k in P}
    W("panel positives:",len(pos),"| gate+audit pass",sum(pbase[k][0]=="pass" for k in pos),"| gate called",sum(pbase[k][1]=="called" for k in pos))
    for k in pos: rowA("panel20",k,P[k]["species"],PT[k],PS.get(k,{}),"%s (gate %s; audit %s %s)"%(pbase[k][0],pbase[k][1],pbase[k][3],pbase[k][4]))
    # part B(a)
    va={a:with_agree(agree_same,lambda a=a:GC.n33(R[a],FE.get(a))) for a in C}
    chg=[a for a in C if key(va[a])!=key(base[a])]
    tsonly=[a for a in C if R[a].get("trny") and R[a]["trny"].get("status")=="called" and any(y[0]=="tRNAscan-SE" and y[2]=="OPPOSITE" for y in strands(R[a]["trny"],R[a]["second"] or {}))]
    W("\n(a) same-strand tRNAscan-SE. genomes whose top arm overlaps an opposite-strand tRNAscan GUA call:",len(tsonly),"| calls that change status:",len(chg),collections.Counter((base[a]["status"],base[a]["N33"],va[a]["status"],M[a]["uaa"]) for a in chg))
    W("   affected among the 17:",[(a,M[a]["species"]) for a in chg if a in hits])
    pa={k:with_agree(agree_same,lambda k=k:panel_status(k,PT[k])) for k in P}; W("   affected panel genomes:",[(k,P[k]["role"],pbase[k][:3],pa[k][:3]) for k in P if pa[k][0]!=pbase[k][0] or pa[k][1]!=pbase[k][1]])
    # part B(b)
    todo=[(a,PATH[a],R[a]["trny"]["reason"]) for a in C if R[a].get("trny") and R[a]["trny"].get("status")=="called" and base[a]["status"]!="resolved" and int(R[a]["trny"].get("n_ok",1))>1]
    ptodo=[(k,f"data/panel/{k}.fna",PT[k]["reason"]) for k in P if PT[k].get("status")=="called" and pbase[k][0]!="pass" and int(PT[k].get("n_ok",1))>1]
    W("\n(b) retry. confirmatory calls not resolved with a second plausible arm:",len(todo),"| panel:",len(ptodo))
    AR={}; fails=[]
    with ProcessPoolExecutor(8) as ex:
        for acc,ok,e in ex.map(arms,todo+ptodo,chunksize=4):
            if ok is None: fails.append((acc,e))
            else: AR[acc]=ok
    W("   cmsearch reruns",len(AR),"failures",fails)
    same=sum(1 for a,_,_ in todo if a in AR and AR[a] and abs(AR[a][0]["score"]-float(R[a]["trny"]["score"]))<0.05 and AR[a][0]["start"]==R[a]["trny"]["start"]); W("   reruns whose top arm equals the locked call:",same,"of",len(todo))
    def retry(a,agf):
        for i,h in enumerate(AR[a]):
            d=dict(R[a]); d["trny"]=dict(R[a]["trny"],**h); z=with_agree(agf,lambda:GC.n33(d,FE.get(a)))
            if z["status"]=="resolved": return i,h,z
        return None
    for tag,agf,bs in (("b",locked_agree,base),("ab",agree_same,va)):
        T=csv.writer(open(O+f"locus_gate_retry_{tag}.tsv","w"),delimiter="\t"); T.writerow("acc species phylum uaa locked_status locked_reason top_score top_N33 new_rank new_score new_N33 new_stem new_start new_end new_strand n_agree".split()); n=collections.Counter(); ch=[]
        for a,_,_ in todo:
            if a not in AR: continue
            z=retry(a,agf)
            if z: i,h,y=z; ch.append(a); n[(M[a]["uaa"],base[a]["status"],base[a]["N33"],"->"+y["N33"])]+=1; T.writerow([a,M[a]["species"],M[a]["phylum"],M[a]["uaa"],base[a]["status"],base[a]["reason"],R[a]["trny"]["score"],R[a]["trny"]["N33"],i+1,h["score"],y["N33"],h["stem_pairs"],h["start"],h["end"],h["strand"],y["n_agree"]])
        W(f"({tag}) newly resolved:",len(ch)); [W("    ",k,v) for k,v in sorted(n.items(),key=lambda x:-x[1])]
        st=[a for a in C if M[a]["uaa"]=="Stop"]; res=[a for a in st if bs[a]["status"]=="resolved" or a in ch]
        W(f"   Stop denominator under ({tag}): {len(res)} (locked 28,583); U33/non-U33:",collections.Counter("U" if (bs[a]["status"]=="resolved" and bs[a]["N33"]=="U") or (bs[a]["status"]!="resolved" and retry(a,agf)[2]["N33"]=="U") else "nonU" for a in res))
        W("   new non-U33 Stop:",[(a,M[a]["species"],retry(a,agf)[2]["N33"],retry(a,agf)[1]["score"]) for a in ch if M[a]["uaa"]=="Stop" and retry(a,agf)[2]["N33"]!="U"])
        W("   Tyr genome CM003310.1 retried:", "CM003310.1" in AR, retry("CM003310.1",agf) if "CM003310.1" in AR else "no second plausible arm")
    for k,_,_ in ptodo:
        if k in AR:
            z=[(i+1,h["score"],h["N33"],panel_status(k,dict(PT[k],**h))) for i,h in enumerate(AR[k])]; W("   panel retry",k,P[k]["role"],P[k]["species"],[x for x in z if x[3][0]=="pass"] or "no arm passes")
    # note also: does any of the 17 have a runner-up? (listed for the record)
    W("   n_ok of the 17:",{a:R[a]["trny"].get("n_ok") for a in hits})
    # (c)
    def locusc(d,feats):
        t=d["trny"]; a,b=sorted((int(t["start"]),int(t["end"]))); w,why,dis,has=locked_locus(d,feats)
        for x in (d.get("trna") or {}).get("trnascan",[])+(d.get("trna") or {}).get("aragorn",[]):
            if x[4].upper()=="GTA": continue
            o=GC.ov((a,b),(x[0],x[1]))
            if o>w: w,why=o,"call_"+x[4]+x[2]
        return w,why,dis,has
    locked_locus=GC.locus; GC.locus=locusc; vc={a:GC.n33(R[a],FE.get(a)) for a in C}; GC.locus=locked_locus
    cc=[a for a in C if key(vc[a])!=key(base[a])]; W("\n(c) extra: all non-GUA ARAGORN and tRNAscan-SE calls on either strand added to the overlap test: status changes",len(cc),"| of the 17:",[(a,M[a]["species"],vc[a]["reason"]) for a in cc if a in hits])
    ccs=[a for a in cc if M[a]["uaa"]=="Stop"]; W("   resolved Stop lost:",sum(base[a]["status"]=="resolved" for a in ccs))
