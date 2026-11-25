"""[P] Post-hoc (2026-11-25): pairwise identity of tRNA-Tyr among tylenchids, by boundary set. Checks REPORT 7.2.3 and the preprint.
Boundaries: ARAGORN (-mt -gc5, the study's options; calls read from results/gate3/run/out.jsonl and results/gate3_unverified/run/out.jsonl), the covariance-model call (trny start/end in the same outputs; results/gate2/trny_panel_final.jsonl
for Radopholus), and the GenBank annotation (Radopholus FN313571.1 915-970, strand +). ARAGORN makes no Tyr(GUA) call in Radopholus FN313571.1 (results/gate2/trny_second.jsonl; rerun here), so no identity with Radopholus can be on ARAGORN boundaries for both.
The first figures (2026-11-24; 51/57, 50/58, 38/55, 39/53) used ARAGORN for R. reniformis only and the covariance-model boundaries for every other sequence, in one MAFFT G-INS-i alignment of seven sequences: block "original".
Identity = matches / columns where both have a base; gap columns are given beside it. Pairwise global alignment: Biopython (match 1, mismatch -1, open -2, extend -0.5). Writes results/posthoc2/tylench_trny_identity.tsv"""
import sys,os,csv,json,subprocess,itertools,re
sys.path.insert(0,"src"); import uaa_caller as U
from Bio.Align import PairwiseAligner
def sp(a): return f"data/sweep/fa/{str(sum(map(ord,a))%200).zfill(3)}/{a}.fna"
P={"R_reniformis":("CM003310.1",sp("CM003310.1")),"P_goodeyi":("ON525593.1","data/sweep_unverified/fa/169/ON525593.1.fna"),"Radopholus":("FN313571.1","data/panel/FN313571.1.fna"),"P_vulnus":("GQ332425.1",sp("GQ332425.1")),"M_incognita":("NC_024097.1",sp("NC_024097.1")),"M_graminicola":("NC_056772.1",sp("NC_056772.1")),"Globodera_I":("KU726971.1",sp("KU726971.1"))}
AC={v[0]:k for k,v in P.items()}; AR={};CM={}
for f in ("results/gate3/run/out.jsonl","results/gate3_unverified/run/out.jsonl"):
    for l in open(f):
        if any(a in l[:60] for a in AC):
            d=json.loads(l)
            if d["acc"] in AC:
                k=AC[d["acc"]]; AR[k]=[x[:3] for x in d["trna"]["aragorn"] if x[4]=="GTA"]; t=d["trny"]
                if t and t.get("start"): CM[k]=[*sorted((int(t["start"]),int(t["end"]))),t["strand"]]
for l in open("results/gate2/trny_panel_final.jsonl"):
    d=json.loads(l)
    if d["panel_id"]=="FN313571.1": CM["Radopholus"]=[*sorted((d["start"],d["end"])),d["strand"]]
r=subprocess.run(["env/prefix/bin/aragorn","-mt","-gc5","-w",P["Radopholus"][1]],capture_output=True,text=True).stdout; print("Radopholus, ARAGORN -mt -gc5, calls with anticodon gta:",[l.strip() for l in r.splitlines() if re.match(r"\s*\d+\s+mtRNA",l) and "(gta)" in l] or "none")
def grab(k,c): s=U.load(P[k][1])[1]; x=s[c[0]-1:c[1]]; return U.rc(x) if c[2]=="-" else x
S={}
for k in P:
    if AR.get(k): S[k+"|ARAGORN"]=grab(k,AR[k][0])
    if k in CM: S[k+"|CM"]=grab(k,CM[k])
S["Radopholus|GenBank"]=grab("Radopholus",[915,970,"+"])
for k,v in S.items(): print(k.ljust(24),len(v),(AR if "ARAGORN" in k else CM).get(k.split("|")[0]) if "GenBank" not in k else "915-970",v)
A=PairwiseAligner(); A.mode="global"; A.match_score=1; A.mismatch_score=-1; A.open_gap_score=-2; A.extend_gap_score=-0.5
def ratio(x,y):
    c=[(a,b) for a,b in zip(x,y) if not(a=="-" and b=="-")]; u=[q for q in c if "-" not in q]; return sum(a==b for a,b in u),len(u),len(c)
def mafft(d):
    t=subprocess.run(["env/prefix/bin/mafft","--quiet","--globalpair","--maxiterate","1000","-"],input="".join(f">{k}\n{v}\n" for k,v in d.items()),capture_output=True,text=True,check=True).stdout; o={}
    for b in t.split(">")[1:]: h,q=b.split("\n",1); o[h.strip()]=q.replace("\n","").upper()
    return o
out=[]
orig={"R_reniformis|ARAGORN":S["R_reniformis|ARAGORN"],**{k+"|CM":S[k+"|CM"] for k in ("P_goodeyi","Radopholus","P_vulnus","M_graminicola","M_incognita","Globodera_I")}}; print("P_goodeyi CM call",CM["P_goodeyi"],"(the 2026-11-24 command used 8660-8716)")
orig["P_goodeyi|CM"]=grab("P_goodeyi",[8660,8716,"+"])
sets={"original (R. ren ARAGORN, others CM; MAFFT of 7)":orig,
 "ARAGORN where called, Radopholus on CM (MAFFT)":{k:S[k] for k in ("R_reniformis|ARAGORN","P_goodeyi|ARAGORN","Radopholus|CM","P_vulnus|ARAGORN","M_incognita|ARAGORN")},
 "ARAGORN where called, Radopholus on GenBank (MAFFT)":{k:S[k] for k in ("R_reniformis|ARAGORN","P_goodeyi|ARAGORN","Radopholus|GenBank","P_vulnus|ARAGORN","M_incognita|ARAGORN")},
 "covariance model for all (MAFFT)":{k+"|CM":S[k+"|CM"] for k in ("R_reniformis","P_goodeyi","Radopholus","P_vulnus","M_incognita")}}
for nm,d in sets.items():
    m=mafft(d); print("==",nm)
    for a,b in itertools.combinations(d,2):
        if not({a.split("|")[0],b.split("|")[0]}&{"R_reniformis","P_goodeyi","Radopholus"}): continue
        x=ratio(m[a],m[b]); p=A.align(d[a],d[b])[0]; y=ratio(p[0],p[1])
        out.append(dict(set=nm,a=a,b=b,len_a=len(d[a]),len_b=len(d[b]),mafft=f"{x[0]}/{x[1]}",mafft_pct=round(100*x[0]/x[1],1),mafft_gap_columns=x[2]-x[1],pairwise=f"{y[0]}/{y[1]}",pairwise_pct=round(100*y[0]/y[1],1),pairwise_gap_columns=y[2]-y[1]))
        if {a.split("|")[0],b.split("|")[0]}<={"R_reniformis","P_goodeyi","Radopholus","P_vulnus","M_incognita"}: o=out[-1]; print("  ",a,"vs",b,"| MAFFT",o["mafft"],o["mafft_pct"],"gaps",o["mafft_gap_columns"],"| pairwise",o["pairwise"],o["pairwise_pct"],"gaps",o["pairwise_gap_columns"])
with open("results/posthoc2/tylench_trny_identity.tsv","w",newline="") as f: w=csv.DictWriter(f,list(out[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(out)
