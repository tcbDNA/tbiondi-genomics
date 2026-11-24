"""POST-HOC, not confirmatory: can the stem_audit_review calls of a class (Aves, Testudines, ...) be treated as real tRNAs? Criteria fixed 2026-11-23 before any final result (after interim look 2 showed the bird counts).
A class is re-included only if, over its stem_audit_review calls (n>=10):
 1 location: >=95% of calls in records with an annotated trnY overlap it (>50%, same strand);
 2 corroboration: >=95% have >=2 of 2 independent callers (ARAGORN, tRNAscan-SE) at the same G34 locus;
 3 shared mismatch: >=90% share the same non-pairing stem position(s) (positions 1-5 from the top of the arm) with the same base pair;
 4 not mis-extraction: the 5 bp 5' arm and 3' arm each have >=80% identity to the class consensus, and the loop has the 7-nt form with GUA in >=99%;
 5 sister contrast (reported, not required): share of resolved calls of the class's sister taxa having 5/5 pairs.
Fail any of 1-4 -> calls stay excluded. usage: gate3_posthoc_stem.py  (reads mitogenome_master.tsv and out.jsonl; panel test: --panel)"""
import sys,csv,json,collections
sys.path.insert(0,"src")
WCGU={("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
def arm(path,d):
    s="".join(l.strip() for l in open(path) if l[0]!=">").upper().replace("T","U"); g=int(d["g34_pos"])
    if d["strand"]!="-": w=s[g-8:g+9]
    else: w=s[g-10:g+7].translate(str.maketrans("ACGUN","UGCAN"))[::-1]
    return w if len(w)==17 else ""
def evaluate(name,items):
    """items: list of dict(arm, callers(0-2), annotated(None|bool))"""
    it=[x for x in items if x["arm"]]; n=len(it); out=[f"{name}: n={n}"]
    if n<10: return False,out+["  fewer than 10 calls: not assessed"]
    a=[x["annotated"] for x in it if x["annotated"] is not None]; c1=sum(a)/len(a) if a else 0; c2=sum(x["callers"]==2 for x in it)/n
    mm=collections.Counter(); co5=collections.Counter(); co3=collections.Counter()
    for x in it:
        w=x["arm"]; p=list(zip(w[0:5],w[12:17][::-1])); mm[tuple((i,q) for i,q in enumerate(p) if q not in WCGU)]+=1; co5[w[0:5]]+=1; co3[w[12:17]]+=1
    top,k=mm.most_common(1)[0]; c3=k/n
    def ident(cn,f):
        cons="".join(collections.Counter(f(x)[i] for x in it).most_common(1)[0][0] for i in range(5)); return cons,sum(sum(a==b for a,b in zip(f(x),cons))/5 for x in it)/n
    s5,i5=ident(co5,lambda x:x["arm"][0:5]); s3,i3=ident(co3,lambda x:x["arm"][12:17]); c4=sum(x["arm"][5:12][2:5]=="GUA" for x in it)/n
    ok=(c1>=0.95 if a else False) and c2>=0.95 and c3>=0.90 and i5>=0.8 and i3>=0.8 and c4>=0.99
    out+=[f"  1 annotated trnY agreement {c1:.3f} (n={len(a)})",f"  2 both callers agree {c2:.3f}",f"  3 shared mismatch {c3:.3f} pattern {[(i+1,''.join(q)) for i,q in top]} top arms {mm.most_common(3)}"[:260],f"  4 consensus 5'arm {s5} identity {i5:.3f}; 3'arm {s3} identity {i3:.3f}; GUA in loop {c4:.3f}",f"  VERDICT: {'re-include (post-hoc)' if ok else 'stay excluded'}"]
    return ok,out
if __name__=="__main__":
    if "--panel" in sys.argv:
        import trny_decide as TD; P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; G={json.loads(l)["acc"]:json.loads(l) for l in open("results/work/gate3_paneltest/out.jsonl")}; FE=json.load(open("results/gate2/panel_features_fixed.json")); import gate3_classify as K
        TAX={}
        for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
            f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; TAX[f[1]]=f
        rows=[]
        for k in G:
            z=K.n33(G[k],FE.get(k)); t=G[k]["trny"]
            if z["status"]=="stem_audit_review":
                cl=(TAX.get(P[k]["species"]) or [""]*9)[6]; f=FE.get(k); yy=[y for y in (f or {}).get("trna",[]) if y[4]]; ag=sum(TD.agree(t,G[k]["second"])); s,e=sorted((t["start"],t["end"]))
                rows.append((cl,dict(arm=arm(f"data/panel/{k}.fna",t),callers=ag,annotated=(any(min(e,y[1])-max(s,y[0])>0.5*min(e-s,y[1]-y[0]) for y in yy) if yy else None))))
        d=collections.defaultdict(list); [d[c].append(x) for c,x in rows]
        for c,x in d.items(): print("\n".join(evaluate(c or "(none)",x)[1]))
        sys.exit()
    import trny_decide as TD
    R={json.loads(l)["acc"]:json.loads(l) for l in open("results/gate3/run/out.jsonl")}; FE={json.loads(l)["acc"]:json.loads(l)["f"] for l in open("data/sweep/features_fixed.jsonl")}; rows=[r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included"]
    d=collections.defaultdict(list); cls=lambda r: r["class"] or r["order"]
    for r in rows:
        if r["n33"]=="stem_audit_review":
            t=R[r["acc"]]["trny"]; f=FE.get(r["acc"]); yy=[y for y in (f or {}).get("trna",[]) if y[4]]; s,e=sorted((t["start"],t["end"])); st="+" if t["strand"]!="-" else "-"
            d[cls(r)].append(dict(arm=arm("data/sweep/fa/"+str(sum(map(ord,r["acc"]))%200).zfill(3)+f"/{r['acc']}.fna",t),callers=sum(TD.agree(t,R[r["acc"]]["second"])),annotated=(any(min(e,y[1])-max(s,y[0])>0.5*min(e-s,y[1]-y[0]) and y[2]==st for y in yy) if yy else None),acc=r["acc"],base=r["base"],uaa=r["uaa"],role=r["role"],species=r["species"]))
    with open("results/gate3/posthoc_stem.txt","w") as o:
        print("POST-HOC sensitivity analysis: stem_audit_review re-inclusion. Not confirmatory; criteria in src/gate3_posthoc_stem.py.",file=o); keep=set()
        for c,x in sorted(d.items(),key=lambda z:-len(z[1])):
            ok,out=evaluate(c,x); print("\n".join(out),file=o)
            if ok: keep.add(c); print("   re-included bases",dict(collections.Counter(y["base"] for y in x)),"non-U33 (Stop):",[(y["acc"],y["species"],y["base"]) for y in x if y["base"]!="U" and y["uaa"]=="Stop"],file=o)
        print("re-included classes:",sorted(keep),file=o)
        st=[r for r in rows if r["role"]=="confirmatory" and r["uaa"]=="Stop" and (r["n33"]=="resolved" or (r["n33"]=="stem_audit_review" and cls(r) in keep))]; nu=[r for r in st if r["base"]!="U"]; print(f"T1 with re-included classes: non-U33 {len(nu)}/{len(st)}; primary (resolved only): {sum(r['n33']=='resolved' and r['base']!='U' for r in st)}/{sum(r['n33']=='resolved' for r in st)}",file=o)
    print(open("results/gate3/posthoc_stem.txt").read()[:3000])
