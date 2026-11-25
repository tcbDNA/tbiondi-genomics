"""Build armless nematode trnY CMs from PANEL nematodes' annotated trnY, leave-one-order-out. T-arm pairs removed from the MITOS structure."""
import json,csv,subprocess,sys,collections,re
sys.path.insert(0,"src"); import panel_guard
P={r["panel_id"]:r for r in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}
A=json.load(open("results/gate2/trny_diag/annotated_trny_panel.json")); B="env/prefix/bin/"; D="results/gate2/trny_diag/"
nem={k:v[0]["seq"] for k,v in A.items() if P[k]["lineage"]=="Nematoda" and v and v[0]["arms"] and v[0]["arms"][0][0]>=4}; panel_guard.check(nem)
MT={"NC_023208.1":"Tylenchina","NC_021120.1":"Tylenchina","NC_056772.1":"Tylenchina","NC_026556.1":"Tylenchina","NC_024097.1":"Tylenchina","NC_025291.1":"Tylenchina","NC_061210.1":"Tylenchina","NC_001328.1":"Rhabditina"}
TAX={}
for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
    f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; TAX.setdefault(f[1],f[5])
order={k:(TAX.get(P[k]["species"]) or P[k]["note"].split("/")[-1] or "?") for k in nem}
print(len(nem),collections.Counter(order.values()))
folds=collections.defaultdict(list)
for k,o in order.items(): folds[o].append(k)
json.dump({"order":order,"folds":folds},open(D+"nem_folds.json","w"))
def build(name,keys):
    open(D+name+".fa","w").write("".join(f">{k}\n{nem[k]}\n" for k in keys))
    subprocess.run([B+"cmalign","--notrunc","-o",D+name+".sto","data/mitos_ref/refseq89m/ncRNA/trnY.cm",D+name+".fa"],check=True,capture_output=True)
    L=open(D+name+".sto").read().split("\n"); ss="".join(l.split()[2] for l in L if l.startswith("#=GC SS_cons"))
    # T-arm = last '<...>' hairpin: blank it
    hp=[m.span() for m in re.finditer(r"<[<.\-]*[_,.~]+[>.\-]*>",ss)]
    a,b=hp[-1]; new=ss[:a]+re.sub(r"[<>]",".",ss[a:b])+ss[b:]; i=0; out=[]
    for l in L:
        if l.startswith("#=GC SS_cons"): w=len(l.split()[2]); out.append(l.replace(l.split()[2],new[i:i+w])); i+=w
        elif l.startswith("#=GC RF") or l.startswith("#=GR"): continue
        else: out.append(l)
    open(D+name+".2.sto","w").write("\n".join(out)); subprocess.run([B+"cmbuild","-F",D+name+".cm",D+name+".2.sto"],check=True,capture_output=True)
for o,ks in folds.items(): build("nem_minus_"+re.sub(r"\W","",o),[k for k in nem if k not in ks])
build("nem_all",list(nem))
r=subprocess.run([B+"cmcalibrate","--forecast","--cpu","8",D+"nem_all.cm"],capture_output=True,text=True); print(r.stdout[-700:])
print(subprocess.run([B+"cmemit","-c",D+"nem_all.cm"],capture_output=True,text=True).stdout)
