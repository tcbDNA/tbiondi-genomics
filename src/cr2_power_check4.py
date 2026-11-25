"""Cold review 2, checks 4 and 5 (POST-HOC): exact one-sided 95% upper bounds for the competitive counts in results/posthoc2/secondcopy/summary.tsv. writes results/posthoc2/read_power_check4.tsv"""
import csv
from scipy.stats import beta
w=csv.writer(open("results/posthoc2/read_power_check4.tsv","w"),delimiter="\t"); w.writerow("acc run record_base record_reads U33_reads n upper95_U33_share supports_1pct".split())
for r in csv.DictReader(open("results/posthoc2/secondcopy/summary.tsv"),delimiter="\t"):
    k=int(r[r["record"]]); u=int(r["U"]); n=k+u; U=beta.ppf(.95,u+1,n-u); w.writerow([r["acc"],r["run"],r["record"]+"33",k,u,n,f"{U:.4f}","yes" if U<=.01 else "no"])
