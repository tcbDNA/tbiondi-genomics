"""EXPLORATORY / POST-HOC. Builds the UNVERIFIED-flagged species set in its own directories (data/sweep_unverified, results/gate3_unverified). Uses the same filters and the fixed annotation parser. Does not touch Gate 3 files."""
import csv,glob,json,os,re,sys,warnings,collections; warnings.filterwarnings("ignore")
sys.path.insert(0,"src"); from Bio import SeqIO; import gate3_features_fix as FF
TAX={}
for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
    f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]; TAX[f[0]]=f
REP={r["accession"]:r for r in csv.DictReader(open("results/insdc_only_set.tsv"),delimiter="\t") if r["eligible"]=="no" and "unverified" in r["title"].lower()}
P=list(csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")); bi=lambda s:" ".join(re.sub(r"[^A-Za-z. ]"," ",s).split()[:2]).lower(); PSP={bi(r["species"]) for r in P}|{"plicatula australis","anomia chinensis","placuna placenta"}; PID={r["panel_id"] for r in P}
KNOWN={"anomia","enigmonia","placuna","plicatula","cephalodiscus","radopholus","pseudotanais"}; BAN={"PP708078","PP708091","OZ576062"}
done=set(); man=[]; aud=[]; F=open("data/sweep_unverified/features.jsonl","w"); os.makedirs("data/sweep_unverified/fa",exist_ok=True)
for fn in sorted(glob.glob("data/sweep_unverified/gb/*.gb")):
    for r in SeqIO.parse(fn,"genbank"):
        a=r.id
        if a not in REP or a in done: continue
        done.add(a); R=REP[a]
        try: s=re.sub(r"[^ACGTN]","N",str(r.seq).upper())
        except Exception: s=""
        L=len(s); nN=s.count("N"); ex=""
        if not s: ex="no_sequence"
        elif a.split(".")[0] in BAN: ex="excluded_by_prereg"
        elif not 8000<=L<=200000: ex="length_outside_8_200kb"
        elif nN>0.01*L: ex="ambiguous_gt1pct"
        t=TAX.get(R["taxid"]) or [""]*9; sp=t[2] or t[1] or R["species"]
        role="panel" if a in PID else "panel_species" if bi(sp) in PSP or bi(R["species"]) in PSP else "known_reassigned_genus" if sp.split(" ")[0].lower() in KNOWN else "unverified_exploratory"
        if not ex:
            sh="data/sweep_unverified/fa/"+str(sum(map(ord,a))%200).zfill(3); os.makedirs(sh,exist_ok=True); p=f"{sh}/{a}.fna"; open(p,"w").write(f">{a}\n{s}\n"); man.append([a,p,r.annotations.get("topology","")]); F.write(json.dumps(dict(acc=a,f=FF.feats(r)))+"\n")
        aud.append([a,R["species"],sp,R["taxid"],t[7],t[6],t[5],t[4],t[3],"UNVERIFIED",r.annotations.get("topology",""),L,nN,role,"included" if not ex else ex,r.annotations.get("comment","")[:0]])
F.close(); [aud.append([a,R["species"],R["species"],R["taxid"],"","","","","","UNVERIFIED","","","","","excluded:not_fetched",""]) for a,R in REP.items() if a not in done]
w=csv.writer(open("results/gate3_unverified/manifest.tsv","w"),delimiter="\t"); w.writerow(["acc","path","topology"]); w.writerows(man)
w=csv.writer(open("results/gate3_unverified/sweep_audit.tsv","w"),delimiter="\t"); w.writerow("acc organism species taxid phylum class order family genus provenance topology length nN role status note".split()); w.writerows(aud)
print("UNVERIFIED species",len(REP),"included",len(man),collections.Counter(x[14] for x in aud),"roles",collections.Counter(x[13] for x in aud if x[14]=="included"))
