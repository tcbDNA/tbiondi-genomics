"""Cold review 2 (POST-HOC, found during check 1): Bufotes pewzowi NC_047225.1 tRNA-Tyr. ARAGORN reads Tyr(GUA), tRNAscan-SE reads Thr(GGU) at the same coordinates and strand. Compares the two framings and the congeners. writes results/posthoc2/bufotes_framing.txt"""
import csv,json
RC=str.maketrans("ACGTN","TGCAN"); rc=lambda s:s.translate(RC)[::-1]; o=open("results/posthoc2/bufotes_framing.txt","w"); W=lambda *x:(print(*x),print(*x,file=o)); PAIR={("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
M={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t") if r["status"]=="included" and (r["genus"] in("Bufotes","Bufo","Duttaphrynus","Strauchbufo"))}; R={}
for l in open("results/gate3/run/out.jsonl"):
    x=json.loads(l)
    if x["acc"] in M: R[x["acc"]]=x
W("POST-HOC")
for a,m in sorted(M.items(),key=lambda z:(z[1]["genus"]!="Bufotes",z[1]["species"])):
    x=R[a]; t=x["trny"]
    if t.get("status")!="called": W(a,m["species"],t.get("status")); continue
    s="".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,a))%200).zfill(3)+f"/{a}.fna") if l[0]!=">").upper(); u,v=sorted((t["start"],t["end"])); q=s[u-1:v]; q=rc(q) if t["strand"]=="-" else q; q=q.replace("T","U"); g=(v-t["g34_pos"]) if t["strand"]=="-" else t["g34_pos"]-u
    ts=[c for c in x["trna"]["trnascan"] if c[0]<v and c[1]>u]
    arm=lambda p:(q[p-7:p-2],q[p-2:p+5],q[p+5:p+10],sum((i,j) in PAIR for i,j in zip(q[p-7:p-2],q[p+5:p+10][::-1])))
    W(a,m["species"],m["n33"],m["base"],"| score",t["score"],"| Tyr framing",arm(g),"| tRNAscan-SE at the locus:",[(c[3],c[4],c[5],c[2]) for c in ts],"|",q)
    if a=="NC_047225.1": W("    one-base-shifted framing (anticodon GGU):",arm(g-1))
