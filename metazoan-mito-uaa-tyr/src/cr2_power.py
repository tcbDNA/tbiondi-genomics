"""Cold review 2, check 5 / item 11 (POST-HOC): exact 95% upper bound (one-sided Clopper-Pearson; 1-0.05^(1/n) when the competing count is 0) on the share of the competing allele (U33 against the record base) for every read check of the 17 hits. writes results/posthoc2/read_power.tsv"""
import csv
from scipy.stats import beta
R=[("NC_037360.1","M. kubotai","C","SRR39401262","independent WGS, exact 31-mer",572,1,0),("NC_037360.1","M. kubotai","C","SRR4556446","amplicon, exact 31-mer",10,0,0),
("CM165726.1","A. glossodonta","C","SRR13036352","independent RNA-seq, aligned",1797,1,0),("CM165726.1","A. glossodonta","C","SRR38145622","own specimen HiFi, aligned",32,0,0),
("NC_013627.1","S. beanii","C","SRR36227421","independent, record 31-mer",0,0,0),("NC_013627.1","S. beanii","C","SRR36227421","independent, own-flank 31-mer",272,0,0),
("NC_013436.1","S. sector","C","SRR22396609","independent, exact 31-mer",363,0,0),("NC_013628.2","St. hypomelas","C","SRR38094758","independent, exact 31-mer",85,0,1),
("CM163399.1","L. verrillii","C","SRR37037826","independent, exact 31-mer",622,1,0),("PZ495457.1","S. lanceolatoides","C","SRR38368693","own specimen, exact 31-mer",419,0,0),
("PV357903.1","C. meadi","C","SRR25465297","own specimen, exact 31-mer",88,0,0),("PX804232.1","O. asakusae","C","SRR36182985","own specimen, exact 31-mer",200,0,0),
("KU726971.1","G. ellingtonae","C","SRR15202050","population 1, exact 31-mer",176,0,1),("KU726971.1","G. ellingtonae","C","SRR15202047","population 2, exact 31-mer",271,0,0),
("OZ318138.1","C. infrequens","G","SRR33975184","independent RNA-seq, aligned >=98%/60 bp",49,0,0),("OZ318138.1","C. infrequens","G","SRR33975184","independent RNA-seq, aligned >=90%/45 bp",55,0,0),("OZ318138.1","C. infrequens","G","SRR33975184","independent RNA-seq, exact 31-mer",66,0,0),
("OZ318138.1","C. infrequens","G","ERR14835550","own specimen HiFi, exact 31-mer",40,0,0),("OZ318138.1","C. infrequens","G","ERR14792883","RNA-seq, exact 31-mer",2,0,0),
("OZ407450.1","C. zetlandica","G","ERR14842678","own specimen HiFi, exact 31-mer",55,0,0),("OZ407450.1","C. zetlandica","G","ERR14857142","own specimen RNA-seq, exact 31-mer",35,0,0),
("NC_047225.1","B. pewzowi","G","SRR9586127","RAD, exact 31-mer",0,0,0),("NC_047225.1","B. pewzowi","G","SRR9586188","RAD, exact 31-mer",0,0,0)]
w=csv.writer(open("results/posthoc2/read_power.tsv","w"),delimiter="\t"); w.writerow("acc species record_base run kind record_reads U33_reads other_base_reads n observed_U33_share upper95_U33_share supports_1pct_claim".split())
for a,sp,b,run,kind,k,u,x in R:
    n=k+u
    if n==0: up="1 (no reads)"; ok="no (uninformative)"; ob="NA"
    else: U=float(beta.ppf(0.95,u+1,n-u)); up=f"{U:.4f}"; ok="yes" if U<=0.01 else "no"; ob=f"{u/n:.4f}"
    w.writerow([a,sp,b+"33",run,kind,k,u,x,n,ob,up,ok]); print(f"{sp:18s} {run:12s} {kind:42s} {k}/{u}  upper {up}  <=1%: {ok}")
print("reads needed for an upper bound of 1% with no U33 read:",next(n for n in range(1,1000) if 1-0.05**(1/n)<=0.01),"; with one U33 read:",next(n for n in range(2,2000) if beta.ppf(.95,2,n-1)<=0.01))
print("no read check at all: M. nanus NC_015546.1, Pterothrissus NC_005796.1, Xyrias KY021065.1, Toxocara NC_082198.1 (no SRA data)")
