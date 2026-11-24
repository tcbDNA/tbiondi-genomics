"""Gate 4: collect mtRF1a (MTRF1L) candidate proteins. Proteomes by blastp (as Gate 2 mtrf_candidates.py); assemblies and TSA sets by miniprot against RF1 references.
Output: data/mtrf/g4/cand.faa (>SOURCE|id). Needs GATE3_UNLOCK=1 (mtrf_guard). usage: gate4_models.py [SOURCE ...] (default: all)"""
import sys,os,subprocess,re,gzip,collections
sys.path.insert(0,"src"); import mtrf_guard
B="env/prefix/bin/"; G="data/mtrf/g4"; W="results/work/g4"; os.makedirs(W,exist_ok=True); os.makedirs("results/gate4",exist_ok=True)
REF=G+"/refs.faa"; N="results/work/numt"
def fa(p):
    S={};c=None
    for l in (gzip.open(p,"rt") if p.endswith(".gz") else open(p)):
        if l[0]==">": c=l[1:].split()[0]; S[c]=[]
        else: S[c].append(l.strip())
    return {k:"".join(v) for k,v in S.items()}
if not os.path.exists(REF):
    R=fa("data/mtrf/miniprot_refs.faa"); H=fa("data/mtrf/human_refs.faa"); L={}
    import csv
    for r in csv.DictReader(open("results/work/mtrf_labels.tsv"),delimiter="\t"): L[r["protein"]]=r["best_human"]
    C=fa("data/mtrf/candidates.faa"); n=0
    for k,s in C.items():
        if L.get(k)=="HUMAN_MTRF1L" and 250<=len(s)<=500 and k not in R: R[k]=s; n+=1
    R["HUMAN_MTRF1L"]=H["HUMAN_MTRF1L"]; open(REF,"w").write("".join(f">{k}\n{s}\n" for k,s in R.items())); print("refs",len(R))
PROT={"PROT_Albula_glossodonta_BYU":"Albula_glossodonta_BYU.faa","PROT_Megalops":"Megalops_cyprinoides.faa","PROT_Anguilla":"Anguilla_anguilla.faa","PROT_Lycodopsis":"Lycodopsis_pacificus.faa","PROT_Zoarces":"Zoarces_viviparus.faa","PROT_M_hapla":"Meloidogyne_hapla.faa"}
GEN={"GEN_R_reniformis":G+"/Rotylenchulus_reniformis.fna","GEN_P_vulnus":G+"/Pratylenchus_vulnus.fna","GEN_H_columbus":G+"/Hoplolaimus_columbus.fna","GEN_R_robustus":G+"/Rotylenchus_robustus.fna","GEN_C_apion":G+"/Cinachyrella_apion.fna",
 "GEN_A_glossodonta_hap1":N+"/GCA_056783665.1.fna","GEN_Lycenchelys":N+"/GCA_056747025.1.fna","GEN_Cr_infrequens":N+"/GCA_974939545.1.fna","GEN_Cr_zetlandica":N+"/GCA_978397455.1.fna"}
TSA=[m for m in("GGVO","GGVQ","GGVP","GGVR","GGVV")]
want=sys.argv[1:]; out={}
def sel(k): return not want or k in want
for k,f in PROT.items():
    if not sel(k): continue
    mtrf_guard.check([k]); p=f"{G}/{f}"; db=f"{W}/{k}"; S=fa(p)
    subprocess.run([B+"makeblastdb","-in",p,"-dbtype","prot","-out",db],capture_output=True)
    r=subprocess.run([B+"blastp","-query",REF,"-db",db,"-evalue","1e-8","-max_target_seqs","20","-num_threads","4","-outfmt","6 sseqid bitscore"],capture_output=True,text=True).stdout; best=collections.Counter()
    for l in r.splitlines(): s,b=l.split("\t"); best[s]=max(best[s],float(b))
    for s,b in best.most_common(8): out[f"{k}|{s}"]=S[s].rstrip("*")
    print(k,"proteins",len(S),"hits",len(best),flush=True)
def mp(src,fna,label):
    mtrf_guard.check([src]); o=f"{W}/{src}.gff"
    if not os.path.exists(o):
        r=subprocess.run([B+"miniprot","-t","6","--gff","--trans","-G","60k",fna,REF],capture_output=True,text=True); open(o,"w").write(r.stdout)
    H=[];cur=None
    for l in open(o):
        if l.startswith("##PAF"): cur=l.rstrip("\n").split("\t")
        elif l.startswith("##STA") and cur:
            f=cur; H.append((float(re.search(r"AS:i:(-?\d+)",l0)[1]) if False else 0,f[1],f[6],f[5],int(f[7+0]) if False else int(f[7]) if False else 0,l.split("\t")[1].strip(),f))
    res=[]
    for _,q,ctg,st,_,prot,f in H:
        try: score=float(f[12] if False else re.search(r"\tAS:i:(-?\d+)",chr(9).join(f))[1])
        except Exception: score=0
        res.append((score,q,f[6],f[5],int(f[8]),int(f[9]),prot.rstrip("*")))
    res.sort(key=lambda x:-x[0]); used=[]; n=0
    for sc,q,c,st,a,b,pr in res:
        if len(pr)<100: continue
        k=[u for u in used if c==u[0] and min(b,u[2])-max(a,u[1])>0]
        if len(k)>=4 or (k and any(u[3]==q for u in k)): continue
        used.append((c,a,b,q)); out[f"{src}|{c}:{a}-{b}({st})~{q[-14:]}"]=pr; n+=1
        if n>=40: break
    print(src,"alignments",len(res),"loci kept",n,flush=True)
for k,f in GEN.items():
    if sel(k): mp(k,f,k)
if sel("TSA_R_reniformis"):
    f=f"{W}/TSA_R_reniformis.fna"
    if not os.path.exists(f):
        with open(f,"w") as o:
            for m in TSA:
                for k,s in fa(f"{G}/tsa.{m}.fsa_nt.gz").items(): o.write(f">{k}\n{s}\n")
    mp("TSA_R_reniformis",f,"tsa")
mode="a" if want else "w"; S=fa(G+"/cand.faa") if want and os.path.exists(G+"/cand.faa") else {}
S.update(out); open(G+"/cand.faa","w").write("".join(f">{k}\n{s}\n" for k,s in S.items())); print("candidates total",len(S),"new",len(out),flush=True)
