"""Figures for REPORT.md (SVG, PDF, PNG). fig1 tree of non-U33 Stop and UAA=Tyr lineages; fig2 mtRF1a loop alignment; fig3 detection sensitivity. All numbers come from committed tables; sources are printed."""
import csv,json,os,subprocess,sys,math,collections
sys.path.insert(0,"src")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"svg.fonttype":"none","pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
O="reports/figures/"; BL,OR,AQ,INK,MUT,SUR="#2a78d6","#eb6834","#1baf7a","#0b0b0b","#52514e","#fcfcfb"
def save(fig,n):
    for e in("svg","pdf","png"): fig.savefig(f"{O}{n}.{e}",dpi=300,bbox_inches="tight",facecolor="white")
    plt.close(fig)
# ---------- fig1 ----------
def lineages(names):
    out={}; want={n.lower() for n in names}; keep={}
    for l in open("data/ncbi/taxonomy/rankedlineage.dmp"):
        f=[x.strip() for x in l.rstrip("\n").rstrip("|").split("\t|\t")]
        if f[1].lower() in want: keep[f[1].lower()]=f
    for n in names:
        f=keep.get(n.lower()); out[n]=dict(species=n,genus=f[3],family=f[4],order=f[5],cls=f[6],phylum=f[7]) if f else None
    return out
TI={r["acc"]:(" ".join(r["tier"].split()[:2]) if r["tier"].startswith("Tier") else r["tier"].split()[0]) for r in csv.DictReader(open("results/cold_review/tiers.tsv"),delimiter="\t")}
SP={"NC_037360.1":"Microdevario kubotai","NC_015546.1":"Microdevario nanus","CM165726.1":"Albula glossodonta","NC_005796.1":"Pterothrissus gissu","NC_013627.1":"Serrivomer beanii","NC_013436.1":"Serrivomer sector","PZ495457.1":"Serrivomer lanceolatoides","NC_013628.2":"Stemonidium hypomelas","CM163399.1":"Lycenchelys verrillii","PV357903.1":"Coloconger meadi","PX804232.1":"Ophichthus asakusae","KY021065.1":"Xyrias revulsus","NC_082198.1":"Toxocara apodemi","KU726971.1":"Globodera ellingtonae","OZ407450.1":"Craniella zetlandica","OZ318138.1":"Craniella infrequens","NC_047225.1":"Bufotes pewzowi"}
BA={r["acc"]:r["base"][0] for r in csv.DictReader(open("results/cold_review/tiers.tsv"),delimiter="\t")}
T=[(SP[a],BA[a],{"unconfirmed":"unc","dropped":"drop"}.get(TI[a],TI[a])) for a in SP]+[("Rotylenchulus reniformis","T","cand"),("Anomia chinensis","T","panelC"),("Enigmonia aenigmatica","T","panelC"),("Placuna vitream","T","panelC"),("Plicatula muricata","T","panelC"),("Cephalodiscus hodgsoni","T","panelC"),("Radopholus similis","T","panelU"),("Pseudotanais sp.","T","panelU"),("Pratylenchus goodeyi","T","weak")]
# calibration positives: how many of the species' panel genomes pass the molecular gates (tRNA-Tyr gate and locus audit); computed from the post-hoc corroboration table
PG=collections.defaultdict(list)
for r_ in csv.DictReader(open("results/posthoc2/locus_gate_corroboration.tsv"),delimiter="\t"):
    if r_["set"]=="panel20":
        k=r_["species"].replace("Plicatuala","Plicatula"); k="Pseudotanais sp." if "Pseudotanais" in k else "Anomia chinensis" if "Anomia_chinensis" in k else "Placuna vitream" if "vitream" in k else "Enigmonia aenigmatica" if "Enigmonia" in k else k; PG[k].append(r_["locked_status"].startswith("pass"))
PG={k:(sum(v),len(v)) for k,v in PG.items()}; print("panel gate passes",PG)
names=[t[0] for t in T]; L=lineages(names); miss=[n for n in names if L[n] is None]; print("no taxonomy line:",miss)
if "Pseudotanais sp." in miss: L["Pseudotanais sp."]=dict(species="Pseudotanais sp.",genus="Pseudotanais",family="Pseudotanaidae",order="Tanaidacea",cls="Malacostraca",phylum="Arthropoda")
if "Placuna vitream" in miss: L["Placuna vitream"]=dict(species="Placuna vitream",genus="Placuna",family="Placunidae",order="Pectinida",cls="Bivalvia",phylum="Mollusca")
RK=["phylum","cls","order","family","genus"]
rows=sorted(T,key=lambda t:[L[t[0]][r] or "~" for r in RK]+[t[0]]); Y={t[0]:len(rows)-1-i for i,t in enumerate(rows)}
fig,ax=plt.subplots(figsize=(9.2,8.2)); X={r:i for i,r in enumerate(["root"]+RK+["tip"])}; ax.set_xlim(-0.2,11.6); ax.set_ylim(-8.2,len(rows)-0.2)
def node(r,key): return tuple((L[n][k]) for k in RK[:RK.index(r)+1]) if False else None
groups=collections.OrderedDict(); ys={}
for depth,r in enumerate(RK):
    g=collections.OrderedDict()
    for t in rows: g.setdefault(tuple(L[t[0]][k] for k in RK[:depth+1]),[]).append(Y[t[0]])
    ys[r]={k:sum(v)/len(v) for k,v in g.items()}; groups[r]=g
for depth,r in enumerate(RK):                           # vertical connectors at each rank and horizontals to children
    prev=ys[RK[depth-1]] if depth else None
    for k,y in ys[r].items():
        x=X[r]
        if depth==0: ax.plot([0,x],[y,y],color=MUT,lw=0.8)
        else:
            py=prev[k[:-1]]; ax.plot([X[RK[depth-1]],x],[y,y],color=MUT,lw=0.8) if False else None
            ax.plot([X[RK[depth-1]],x],[y,y],color=MUT,lw=0.8); ax.plot([X[RK[depth-1]]]*2,[py,y],color=MUT,lw=0.8)
for p,v in groups["phylum"].items(): ax.plot([0,0],[min(v),max(v)],color=MUT,lw=0)  # spacing only
ax.plot([0,0],[min(ys["phylum"].values()),max(ys["phylum"].values())],color=MUT,lw=0.8)
for t in rows:
    n,b,s=t; y=Y[n]; gk=tuple(L[n][k] for k in RK); ax.plot([X["genus"],X["tip"]-0.1],[y,y],color=MUT,lw=0.8) if abs(ys["genus"][gk]-y)<1e-9 else (ax.plot([X["genus"],X["genus"]],[ys["genus"][gk],y],color=MUT,lw=0.8),ax.plot([X["genus"],X["tip"]-0.1],[y,y],color=MUT,lw=0.8))
    col={"C":BL,"G":OR,"T":AQ}[b]; mk=dict(color=col,ms=7,mew=1.4); 
    if s=="Tier 1": ax.plot(6.0,y,"o",**mk)
    elif s in("panelC","panelU"):
        a_,b_=PG[n]; ax.plot(6.0,y,"^",ms=8,color=col,mfc=col if a_==b_ else "white" if a_==0 else col,mew=1.4,fillstyle="full" if a_ in(0,b_) else "left",markerfacecoloralt="white",mec=col)
    elif s=="Tier 2": ax.plot(6.0,y,"D",ms=6,color=col)
    elif s=="Tier 3": ax.plot(6.0,y,"s",ms=6,color=col)
    elif s=="unc": ax.plot(6.0,y,"o",ms=7,mfc="white",mec=col,mew=1.4)
    elif s=="drop": ax.plot(6.0,y,"x",ms=7,mec=col,mew=1.8)
    elif s=="cand": ax.plot(6.0,y,"*",ms=11,mfc="white",mec=col,mew=1.4)
    elif s=="weak": ax.plot(6.0,y,"o",ms=7,mfc="white",mec=col,mew=1.4,ls="none")
    tag={("C","Tier 1"):"C33 Stop, Tier 1 same-species independent",("C","Tier 2"):"C33 Stop, Tier 2 own-specimen reads",("C","Tier 3"):"C33 Stop, Tier 3 lineage-backed only",("C","unc"):"C33 Stop, unconfirmed",("C","drop"):"C33 Stop, dropped",("G","Tier 1"):"G33 Stop, Tier 1 same-species independent",("G","Tier 2"):"G33 Stop, Tier 2 own-specimen reads",("G","unc"):"G33 Stop, unconfirmed",("T","panelC"):"UAA=Tyr, panel, C33",("T","panelU"):"UAA=Tyr, panel, U33",("T","cand"):"UAA=Tyr, U33; unconfirmed candidate, stem-blocked",("T","weak"):"UAA=Tyr?, U33; supported weak (UNVERIFIED)"}[(b,s)]
    if s in("panelC","panelU"): a_,b_=PG[n]; tag+=f"; {a_} of {b_} pass the gates"+{"Enigmonia aenigmatica":" (audit: overlap)","Plicatula muricata":" (stem 3/5)","Radopholus similis":" (no 2nd caller)"}.get(n,"")
    ax.text(6.25,y,n,va="center",style="italic",color=INK,fontsize=7.5); ax.text(9.0,y,tag,va="center",color=MUT,fontsize=6.6)
for r in RK: ax.text(X[r],len(rows)-0.35,{"cls":"class"}.get(r,r),ha="center",va="bottom",fontsize=6.5,color=MUT)
for p,y in ys["phylum"].items(): ax.text(0.04,y+0.16,p[0],fontsize=6.2,color=MUT)
txt=["Screen positives (blue C33, orange G33; post-hoc tiers, REPORT 7.2.1): circle Tier 1, diamond Tier 2, square Tier 3, hollow circle unconfirmed, cross dropped. Hollow star: unconfirmed candidate.",
"Calibration positives (green triangles, not new): filled, all panel genomes shown pass the tRNA-Tyr gate and audit; half-filled, some pass; hollow, none pass (N33 is then a raw call).",
"C33 + UAA=Stop: Tier 1, 6 genomes (4 fish lineages, 4 origins); Tier 2, 3; Tier 3, 2; unconfirmed, 2; dropped, 1 (14). Tiers 1 and 2: 9 genomes, 6 origins.",
"G33 + UAA=Stop: Tier 1, 1 genome; Tier 2, 1 (2 species of 1 genus, Craniella; 1 origin); unconfirmed, 1 (Bufotes) (3).",
"11 gains on the NCBI hierarchy, all 17 screen positives (one dropped and three unconfirmed included).",
"UAA=Tyr, panel, about 5 origins. Supported weak: P. goodeyi (hollow green circle, UNVERIFIED record). Non-panel Tyr: Rotylenchulus, unconfirmed single-genome candidate.",
"Background not drawn: 28,566 U33 Stop genomes. Origin counts come from the full taxonomy table, not from this pruned display.",
"The tree is NCBI taxonomy (phylum to genus), not a molecular phylogeny, so origin counts are approximate; branch lengths have no meaning."]
for i,s in enumerate(txt): ax.text(-0.1,-1.3-i*0.8,s,fontsize=6.6,color=INK if i in(2,3,4) else MUT)
ax.axis("off"); ax.set_title("Non-U33 tRNA-Tyr with UAA=Stop, and UAA=Tyr lineages",loc="left",fontsize=10,color=INK,pad=16); save(fig,"fig1_tree")
# ---------- fig2 ----------
H="EGGVHRVQRVPKTEKQGRVHTSTMTVA"; G={}
for r in csv.DictReader(open("results/gate4/loop.tsv"),delimiter="\t"):
    if r["covered"]=="True" and r["label"]=="MTRF1L": G.setdefault(r["source"],r["loop"])
A=[("ref","Homo sapiens MTRF1L",H,"reference","ref"),
 ("tyr","Anomia chinensis",None,"EGGVHRFQRETRVKQARNFFSCTITIA","Gate 2"),("tyr","Anomia ephippium",None,"EGGVHRFQRETRVKKARTFYSCTITIA","Gate 2"),("tyr","Enigmonia aenigmatica",None,"EGGVQRFQRETRVKVPRVFFSCTITVA","Gate 2"),("tyr","Placuna vitream",None,"EGGVQRFQRETRVKSSRTFFSCTITIA","Gate 2"),
 ("tyr","Cephalodiscus hodgsoni (reads, partial)",None,"HKVQRIPKTENKNRIHSST","Gate 2"),("tyr","Cephalodiscus sp. (reads, partial)",None,"HKVQRIPKTERNNRIHSST","Gate 2"),
 ("tyr","Plicatula muricata",None,"EGGVHRVQRKPVTDRNNKRQTSTITVA","Gate 2"),("tyr","Plicatula plicata (reads)",None,"EGGVHRVQRKPVTDRNNKRQTSTITVA","Gate 4"),
 ("tyr","Radopholus similis (partial; manual splice)",None,"KRLPVTNSKGKGVLHTST","Gate 2"),
 ("tyrS","Rotylenchulus reniformis (candidate)",None,G.get("GEN_R_reniformis"),"Gate 4"),
 ("stop","Pratylenchus vulnus",None,G.get("GEN_P_vulnus"),"Gate 4"),("stop","Meloidogyne hapla",None,G.get("PROT_M_hapla"),"Gate 4"),("stopq","Rotylenchus robustus",None,G.get("GEN_R_robustus"),"Gate 4"),
 ("stop","Heteranomia squamula",None,"EGGVHRVIRVPSTSKADKVHTSTVHVA","Gate 2"),("stop","Pododesmus macrochisma",None,"EGGVHRVQRVPVTEQKGRVHTSTITLS","Gate 2"),("stop","Rhabdopleura annulata (reads)",None,"EGGVHRVKRVPKTERKGRIHTSTSSVA","Gate 2"),
 ("tb","Albula glossodonta",None,G.get("PROT_Albula_glossodonta_BYU"),"Gate 4"),("tb","Lycenchelys verrillii",None,G.get("GEN_Lycenchelys"),"Gate 4"),("tb","Craniella zetlandica",None,G.get("GEN_Cr_zetlandica"),"Gate 4"),("tb","Craniella infrequens",None,G.get("GEN_Cr_infrequens"),"Gate 4"),
 ("tbs","Megalops cyprinoides (sister, U33)",None,G.get("PROT_Megalops"),"Gate 4"),("tbs","Lycodopsis pacificus (sister, U33)",None,G.get("PROT_Lycodopsis"),"Gate 4"),("tbs","Cinachyrella apion (sister, U33)",None,G.get("GEN_C_apion"),"Gate 4")]
R=[(a[0],a[1],a[2] if a[0]=="ref" else a[3],a[4]) for a in A]; assert all(x[2] for x in R),[x[1] for x in R if not x[2]]
fa="results/work/g4/fig2.faa"; open(fa,"w").write("".join(f">{i}\n{s}\n" for i,(_,n,s,_) in enumerate(R)))
aln=subprocess.run(["env/prefix/bin/mafft","--quiet","--globalpair","--maxiterate","1000",fa],capture_output=True,text=True).stdout; AL={};c=None
for l in aln.splitlines():
    if l[0]==">": c=int(l[1:]); AL[c]=""
    else: AL[c]+=l.strip()
h=AL[0]; pos={};p=195
for i,ch in enumerate(h):
    if ch!="-": p+=1; pos[p]=i
assert h.replace("-","")==H
n=len(R); W=len(h); fig,ax=plt.subplots(figsize=(10.4,6.6)); ax.set_xlim(-15.4,W+3.4); ax.set_ylim(-6.2,n+2.2); lab={"tyr":"UAA=Tyr (Gate 2, calibration)","tyrS":"unconfirmed candidate","stop":"UAA=Stop sisters","stopq":"UAA=Stop sisters","tb":"[E] Tier B: non-U33 Stop taxa","tbs":"[E] Tier B: U33 sisters","ref":""}
for c203 in(203,216):
    x=pos[c203]; ax.add_patch(plt.Rectangle((x-0.5,-0.7),1,n+0.9,color="#f0efec",zorder=0)); ax.text(x,n+0.55,str(c203),ha="center",fontsize=7.5,weight="bold",color=INK)
for k in(196,222): ax.text(pos[k],n+0.55,str(k),ha="center",fontsize=6.5,color=MUT)
# evidence class of the UAA status of each row (REPORT 6): panel call / confirmatory call / exploratory (UNVERIFIED sweep) call / no call in this study (literature or inferred) / unconfirmed candidate
EV={}
for (g_,nm_,*_) in R:
    EV[nm_]="none" if nm_.split(" (")[0] in("Cephalodiscus sp.","Rhabdopleura annulata","Plicatula plicata","Rotylenchus robustus") else "expl" if nm_.startswith("Meloidogyne hapla") else "cand" if nm_.startswith("Rotylenchulus") else "conf" if g_ in("tb","tbs") or nm_.startswith("Pratylenchus vulnus") else "panel"
MK={"panel":("s",INK,5),"conf":("o",INK,5),"expl":("o","white",5),"none":("X","white",6),"cand":("*","white",9)}
prev=None; y=n-1; seps=[]
for i,(g,nm,s,src) in enumerate(R):
    yy=n-1-i; grp=lab[g]
    if prev is not None and grp!=prev and grp: seps.append(yy+0.5)
    prev=grp if grp else prev
    ax.text(-1.5,yy,nm,ha="right",va="center",fontsize=7,style="italic" if i else "normal",color=INK)
    if i: m_=EV[nm]; ax.plot(-0.8,yy,MK[m_][0],ms=MK[m_][2],mfc=MK[m_][1],mec=INK,mew=1.0,clip_on=False)
    for j,ch in enumerate(AL[i]):
        col=INK; bg=None; w="normal"
        if ch!="-" and j==pos[216] and ch=="S": bg=OR; w="bold"
        if ch!="-" and j==pos[216] and ch=="T": bg="#dff3ea"
        if ch!="-" and j==pos[203] and ch!="Q": bg="#cfe0f6"; w="bold"
        if ch!="-" and j==pos[203] and ch=="Q": bg=None
        if bg: ax.add_patch(plt.Rectangle((j-0.46,yy-0.46),0.92,0.92,color=bg,zorder=1))
        ax.text(j,yy,ch,ha="center",va="center",family="DejaVu Sans Mono",fontsize=7.6,color=col,weight=w,zorder=2)
    if i in(0,): pass
for s in seps: ax.plot([-13.8,W],[s,s],color=MUT,lw=0.5)
def block(y0,y1,t): ax.text(-11.3,(y0+y1)/2,t,ha="left",va="center",fontsize=6.6,color=MUT)
yb=lambda g0:[n-1-i for i,(g,*_) in enumerate(R) if lab[g]==g0]
for t in dict.fromkeys(lab.values()):
    v=yb(t)
    if t and v: pass
for k_,(m_,lb_) in enumerate((("panel","calibration panel call"),("conf","confirmatory Stop call"),("expl","exploratory Stop call (UNVERIFIED sweep)"),("none","no mitogenome call here: literature or inferred"),("cand","Tyr call, unconfirmed candidate"))):
    xk=(-15.1,-6.3,1.6,-15.1,0.2)[k_]; yk=-1.35 if k_<3 else -2.15; ax.plot(xk+0.2,yk,MK[m_][0],ms=MK[m_][2],mfc=MK[m_][1],mec=INK,mew=1.0,clip_on=False); ax.text(xk+0.65,yk,lb_,fontsize=6.6,va="center",color=INK)
txt=["Marks beside the names: evidence class of the UAA status. mtRF1a sequences in the [E] block and in the tylenchids are exploratory.",
"Orange cell: S at human 216 (T216S). Green cell: T at 216. Blue cell: residue other than Q at 203. Human numbering NP_061914.3. Radopholus +2 aa (placement approximate).",
"Loop length against tylenchid Stop sisters (PREREG 2.4): Rotylenchulus 27 aa against 26 in R. robustus, P. vulnus and M. hapla (+1 aa); +0 against human.",
"Partial sequences (Cephalodiscus, Radopholus): gaps are missing data. The two P. vulnus models with +14 aa read through a 47-nt intron and are not used (REPORT 6). H. columbus is omitted (bacterial).",
"[E] exploratory block: Albula, Lycenchelys, Craniella (non-U33 Stop) and three sisters. Rotylenchulus is an unconfirmed single-genome candidate (stem-blocked). T191S is not isolable in Anomia/Placuna (GATE_0)."]
for i,s in enumerate(txt): ax.text(-15.1,-2.95-i*0.75,s,fontsize=6.4,color=MUT)
for t,ys_ in [(t,[n-1-i for i,(g,*_) in enumerate(R) if lab[g]==t]) for t in dict.fromkeys(lab.values()) if t]:
    if ys_: ax.text(-11.3,max(ys_)+0.52,t,fontsize=7,weight="bold",color=INK,va="bottom") if False else None
LEN={"Rotylenchulus reniformis (candidate)":"27 aa: +1","Pratylenchus vulnus":"26 aa","Meloidogyne hapla":"26 aa","Rotylenchus robustus":"26 aa","Homo sapiens MTRF1L":"27 aa"}
for i,(g,nm,sq,src) in enumerate(R):
    if nm in LEN: ax.text(W+0.1,n-1-i,LEN[nm],va="center",fontsize=6.8,color=INK if "27 aa: +1" in LEN[nm] else MUT,weight="bold" if "+1" in LEN[nm] else "normal")
ax.axis("off"); ax.set_title("mtRF1a loop (human 196–222): UAA=Tyr lineages, Stop sisters, and Tier B non-U33 Stop taxa",loc="left",fontsize=10,color=INK)
# block labels at the far left as coloured bars
bars=[("ref",""),]
for (gset,tx,cl) in ((("tyr","tyrS"),"UAA = Tyr",AQ),(("stop","stopq"),"UAA = Stop",BL),(("tb","tbs"),"[E] exploratory (Tier B)",OR)):
    v=[n-1-i for i,(g,*_) in enumerate(R) if g in gset]; ax.add_patch(plt.Rectangle((-15.1,min(v)-0.45),0.16,max(v)-min(v)+0.9,color=cl)); ax.text(-14.45,(max(v)+min(v))/2,tx,rotation=90,va="center",ha="center",fontsize=7,color=INK)
save(fig,"fig2_mtrf1a_loop"); print("fig2 alignment columns",W,"pos203",pos[203],"pos216",pos[216])
# ---------- fig3 ----------
def wil(k,n,z=1.96): p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; e=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-e,c+e
B3=list(csv.DictReader(open("results/posthoc2/fig3_bootstrap.tsv"),delimiter="\t")); fr=[5,10,25,50,100]; ser={}
for r_ in B3: ser.setdefault(r_["series"],[]).append((int(r_["called"]),int(r_["observations"]),float(r_["boot_lo"]),float(r_["boot_hi"]),int(r_["genomes"])))
cnt=[x[:2] for x in ser["counter"]]; cod=[x[:2] for x in ser["counter+Codetta(rep 0)"]]; UP0=1-0.025**(1/65); print("exact two-sided 95% upper limit on 0 of 65 genomes",round(UP0,4))
fig,(a1,a2)=plt.subplots(1,2,figsize=(10,4.2),gridspec_kw=dict(width_ratios=[1.15,1])); X=range(5)
for nm,lb,cl,mk,off in (("counter","Counter rule alone (5 replicates per genome; 1 at 100%)",BL,"o",-0.11),("counter+Codetta(rep 0)","With Codetta agreeing (replicate 0)",OR,"s",0.11)):
    d=ser[nm]; y=[v[0]/v[1] for v in d]; lo=[y[i]-v[2] for i,v in enumerate(d)]; hi=[v[3]-y[i] for i,v in enumerate(d)]; lo[0]=0; hi[0]=UP0
    a1.errorbar([x+off for x in X],y,yerr=[lo,hi],color=cl,marker=mk,ms=6.5,lw=0,elinewidth=1.6,capsize=3,label=lb,ls="none")
for x,(k,n) in enumerate(cnt): a1.annotate(f"{k}/{n}"+(" (0 of 65 genomes)" if k==0 else ""),(x-0.11,k/n),textcoords="offset points",xytext=[(16,-13),(-21,-3),(-22,-3),(-24,-3),(-22,-9)][x],ha="center",fontsize=6.8,color=MUT)
a1.set_xticks(list(X)); a1.set_xticklabels([f"{f}%" for f in fr]); a1.set_ylim(-0.12,1.08); a1.set_ylabel("Fraction of genomes called Tyr"); a1.set_xlim(-0.6,4.5); a1.set_xlabel("Conserved Tyr codons replaced by TAA (categories, not to scale)"); a1.legend(frameon=False,loc="lower right",fontsize=7); a1.set_title("Planted-null detection (95% intervals over 65 source genomes)",loc="left",fontsize=9,color=INK)
LO=list(csv.DictReader(open("results/gate2/lolo.tsv"),delimiter="\t")); S=collections.OrderedDict()
for r in LO:
    if r["unit"]=="lineage" and r["proc"]!="N": S.setdefault(r["held_out"],{})[r["proc"]]=(int(r["recalled"]),int(r["genomes"]))
LIN={r_["panel_id"]:r_["lineage"] for r_ in csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")}; END=collections.defaultdict(lambda:[0,0])
for r_ in csv.DictReader(open("results/posthoc2/locus_gate_corroboration.tsv"),delimiter="\t"):
    if r_["set"]=="panel20": l_=LIN[r_["acc"]]; END[l_][1]+=1; END[l_][0]+=r_["locked_status"].startswith("pass")
print("panel pass rate by lineage",{k:tuple(v) for k,v in END.items()},"total",sum(v[0] for v in END.values()),"of",sum(v[1] for v in END.values()))
ln=sorted(S); procs=["Fixed","M","R","END"]; heads={"Fixed":"UAA rule\n(in-sample)","M":"held-out:\nmidpoint","R":"held-out:\n0.8 x min","END":"panel pass rate\nfor the molecular\ngates"}
for i,l in enumerate(ln):
    for j,q in enumerate(procs):
        k,n=END[l] if q=="END" else S[l][q]; f=k/n; a2.add_patch(plt.Rectangle((j-0.45,len(ln)-1-i-0.4),0.9,0.8,facecolor="#2a78d6" if f==1 else "#9cc3ef" if f>=0.5 else "#dbe9fa" if f>0 else "#f0efec",edgecolor="white")); a2.text(j,len(ln)-1-i,f"{k}/{n}",ha="center",va="center",fontsize=7.5,color="white" if f==1 else INK)
a2.text(len(procs)-1,-0.75,f"{sum(v[0] for v in END.values())} / {sum(v[1] for v in END.values())}",ha="center",fontsize=7.5,weight="bold",color=INK); a2.text(0,-0.75,f"{sum(S[l]['Fixed'][0] for l in ln)} / {sum(S[l]['Fixed'][1] for l in ln)}",ha="center",fontsize=7.5,weight="bold",color=INK); a2.text(-0.65,-0.75,"all",ha="right",fontsize=7.2,color=MUT)
a2.set_xlim(-0.6,len(procs)-0.4); a2.set_ylim(-1.1,len(ln)-0.3); a2.set_xticks(range(len(procs))); a2.set_xticklabels([heads[q] for q in procs],fontsize=6.8); a2.set_yticks(range(len(ln))); a2.set_yticklabels(ln[::-1],fontsize=7,style="italic"); a2.tick_params(length=0); a2.set_title("Recall of the 20 known positives (genomes)",loc="left",fontsize=9,color=INK)
for sp in a2.spines.values(): sp.set_visible(False)
fig.text(0.01,-0.03,"Left: 10-100%, percentile bootstrap that resamples the 65 source genomes and keeps their replicates together (10,000 resamples). At 5%, 0 of 65 genomes: upper limit of the two-sided exact (Clopper-Pearson) 95% interval, 5.5%.",fontsize=6.6,color=MUT)
fig.text(0.01,-0.065,"Intervals do not account for relatedness among source genomes. The two series use different replicate sets. Right: the last column is the UAA rule plus the tRNA-Tyr gate and locus audit (PREREG 2.3) on known positives; it does not test",fontsize=6.6,color=MUT)
fig.text(0.01,-0.095,"the novelty or replication parts of PREREG 3. The UAA rule was set with all lineages in view. False Tyr with thresholds from the other phyla (29 held-out phyla): 0 of 340, upper bound 1.1%. Sources: results/gate2, results/posthoc2 [P].",fontsize=6.6,color=MUT)
fig.tight_layout(); save(fig,"fig3_sensitivity"); print("figures:",sorted(os.listdir(O)))
