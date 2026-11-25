"""Cold review 2, check 5 (POST-HOC): which Tier 1-3 assignments would change if the exact 95% bound on the U33 share (not the point share) had to be <= a threshold. Imports src/tiers.py with its writer redirected; nothing is changed. writes results/posthoc2/read_power_tiers.txt"""
import sys,io,contextlib,csv,collections
sys.path.insert(0,"src"); from scipy.stats import beta
src=open("src/tiers.py").read().replace('open("results/cold_review/tiers.tsv","w")','open("/dev/null","w")'); ns={}
with contextlib.redirect_stdout(io.StringIO()): exec(compile(src,"tiers","exec"),ns)
E,BASE,LIN=ns["E"],ns["BASE"],ns["LIN"]; R={r["acc"]:r for r in csv.DictReader(open("results/mitogenome_master.tsv"),delimiter="\t")}
def run(th):
    def ok(r):
        if not r: return False
        n=r[0]+r[1]; return r[0]>=10 and (beta.ppf(.95,r[1]+1,n-r[1])<=th if th else r[1]<=.01*n)
    t={a:"dropped" if e.get("drop") else "Tier 1" if e["ind_rec"] or ok(e["ind_reads"]) else "Tier 2" if ok(e["own"]) else "" for a,e in E.items()}
    for a in E:
        if t[a]=="":
            L=LIN[a]
            if L["n"]>=1 and L["uniform"] and next(iter(L["bases"]))==BASE[a] and any(t[x] in("Tier 1","Tier 2") for x in E if x!=a and R[x]["species"] in L["species"]): t[a]="Tier 3"
    return {a:(v or "unconfirmed") for a,v in t.items()}
b=run(None); out=open("results/posthoc2/read_power_tiers.txt","w"); print("POST-HOC. Tier 1 genomes that rest on reads alone (no second record):",[R[a]["species"] for a,e in E.items() if b[a]=="Tier 1" and not e["ind_rec"]],file=out)
for th in (.01,.02,.05,.10):
    c=run(th); print(f"bound<={th:.0%}:",[(R[a]["species"],b[a],"->",c[a]) for a in E if b[a]!=c[a]],"totals",dict(collections.Counter((BASE[a],c[a]) for a in E)),file=out)
out.close(); print(open("results/posthoc2/read_power_tiers.txt").read())
