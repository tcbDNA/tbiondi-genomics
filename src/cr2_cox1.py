"""Cold review 2, check 8 / item 15 (POST-HOC): reconcile 10,294 (tsv/log) with 10,291 (REPORT 7.2.4). Reads results/gate3/tylench_cox1.tsv; writes results/posthoc2/cox1_reconcile.txt"""
import csv,collections
R=list(csv.DictReader(open("results/gate3/tylench_cox1.tsv"),delimiter="\t")); o=open("results/posthoc2/cox1_reconcile.txt","w"); W=lambda *x:(print(*x),print(*x,file=o))
W("POST-HOC. rows",len(R),"distinct acc",len({r["acc"] for r in R}))
any_=lambda r:int(r["taa_any"])
for lab,sel in (("all",lambda r:True),("len<=2500",lambda r:int(r["len"])<=2500),("cols>=1",lambda r:int(r["cols"])>=1)): W(lab,sum(sel(r) for r in R))
def counts(S,tag):
    G=lambda g:[r for r in S if r["genus"]==g]
    W(f"[{tag}] n={len(S)}"); 
    for g in ("Rotylenchulus","Hoplolaimus","Rotylenchus","Meloidogyne","Heterodera","Helicotylenchus","Scutellonema","Radopholus"): x=G(g); W("   ",g,"records with >=1 internal UAA",sum(any_(r)>=1 for r in x),"of",len(x),"; >=2:",sum(any_(r)>=2 for r in x))
    oth=[r for r in S if r["family"] not in("Hoplolaimidae","Pratylenchidae")]; W("    outside Hoplolaimidae+Pratylenchidae: >=2 UAA",sum(any_(r)>=2 for r in oth),"of",len(oth),"; >=1:",sum(any_(r)>=1 for r in oth))
counts(R,"as shipped")
dup=[a for a,n in collections.Counter(r["acc"] for r in R).items() if n>1]; W("duplicate accessions:",dup)
# the report's denominators are each 1 (Meloidogyne), 1 (Scutellonema) and 2 (outside) lower: list candidates for the 3 removed records
W("species field with odd names:",[ (r["acc"],r["species"],r["family"],r["len"]) for r in R if not r["family"] or not r["species"]])
W("\nFilter that reproduces the report (10,291; Meloidogyne 1,101; Scutellonema 194; outside 7,584): the query of 2026-11-24 kept cols>=1 AND domain score>=60 (one Scutellonema record fails both). The cols>=1 condition alone removes the same three; its printed label said only domain score >=60, and REPORT 7.2.4 states neither.")
F=[r for r in R if int(r["cols"])>=1]; counts(F,"cols>=1"); counts([r for r in R if float(r["domscore"])>=60],"domain score>=60 (does NOT reproduce the report)")
W("records with 0 conserved-Tyr columns:",[(r["acc"],r["species"],r["domscore"]) for r in R if int(r["cols"])<1])
W("the two outside-family records with exactly 1 UAA:",[(r["acc"],r["species"],r["family"],r["taa_any"],r["domscore"]) for r in R if r["family"] not in("Hoplolaimidae","Pratylenchidae") and any_(r)>=1])
W("the retrieved record that did not align (10,295 retrieved, 10,294 aligned):")
import glob,re
ids=set()
for f in glob.glob("data/tylench_cox1/b*.gb"):
    for b in open(f).read().split("\n//"):
        m=re.search(r"^ACCESSION\s+(\S+)",b,re.M)
        if m: ids.add((m[1],re.search(r"^DEFINITION\s+(.*)",b,re.M)[1],re.search(r"^LOCUS\s+\S+\s+(\d+)",b,re.M)[1]))
have={r["acc"] for r in R}; W("  ",[x for x in ids if x[0] not in have],"records on disk",len(ids))
