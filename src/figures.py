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
T=[("Microdevario kubotai","C","sup"),("Microdevario nanus","C","sup"),("Albula glossodonta","C","sup"),("Pterothrissus gissu","C","sup"),("Serrivomer beanii","C","sup"),("Serrivomer sector","C","sup"),("Serrivomer lanceolatoides","C","sup"),("Stemonidium hypomelas","C","sup"),("Lycenchelys verrillii","C","sup"),
   ("Coloconger meadi","C","unc"),("Ophichthus asakusae","C","unc"),("Xyrias revulsus","C","unc"),("Toxocara apodemi","C","unc"),("Globodera ellingtonae","C","drop"),
   ("Craniella zetlandica","G","sup"),("Craniella infrequens","G","sup"),("Bufotes pewzowi","G","drop"),
   ("Rotylenchulus reniformis","T","new"),("Anomia chinensis","T","panel"),("Enigmonia aenigmatica","T","panel"),("Placuna vitream","T","panel"),("Plicatula muricata","T","panel"),("Cephalodiscus hodgsoni","T","panel"),("Radopholus similis","T","panel"),("Pseudotanais sp.","T","panel"),("Pratylenchus goodeyi","T","weak")]
names=[t[0] for t in T]; L=lineages(names); miss=[n for n in names if L[n] is None]; print("no taxonomy line:",miss)
if "Pseudotanais sp." in miss: L["Pseudotanais sp."]=dict(species="Pseudotanais sp.",genus="Pseudotanais",family="Pseudotanaidae",order="Tanaidacea",cls="Malacostraca",phylum="Arthropoda")
if "Placuna vitream" in miss: L["Placuna vitream"]=dict(species="Placuna vitream",genus="Placuna",family="Placunidae",order="Pectinida",cls="Bivalvia",phylum="Mollusca")
RK=["phylum","cls","order","family","genus"]
rows=sorted(T,key=lambda t:[L[t[0]][r] or "~" for r in RK]+[t[0]]); Y={t[0]:len(rows)-1-i for i,t in enumerate(rows)}
fig,ax=plt.subplots(figsize=(9.2,8.2)); X={r:i for i,r in enumerate(["root"]+RK+["tip"])}; ax.set_xlim(-0.2,11.6); ax.set_ylim(-6.0,len(rows)-0.2)
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
    if s in("sup","panel"): ax.plot(6.0,y,"o",**mk)
    elif s=="unc": ax.plot(6.0,y,"o",ms=7,mfc="white",mec=col,mew=1.4)
    elif s=="drop": ax.plot(6.0,y,"x",ms=7,mec=col,mew=1.8)
    elif s=="new": ax.plot(6.0,y,"*",ms=11,color=col)
    elif s=="weak": ax.plot(6.0,y,"o",ms=7,mfc="white",mec=col,mew=1.4,ls="none")
    tag={("C","sup"):"C33, Stop, supported",("C","unc"):"C33, Stop, unconfirmed",("C","drop"):"C33, Stop, dropped",("G","sup"):"G33, Stop, supported",("G","drop"):"G33, Stop, dropped",("T","panel"):"UAA=Tyr, calibration panel",("T","new"):"UAA=Tyr, new; sensitivity (blocked by the stem audit)",("T","weak"):"UAA=Tyr?, supported weak (UNVERIFIED)"}[(b,s)]
    ax.text(6.25,y,n,va="center",style="italic",color=INK,fontsize=7.5); ax.text(9.0,y,tag,va="center",color=MUT,fontsize=6.6)
for r in RK: ax.text(X[r],len(rows)-0.35,{"cls":"class"}.get(r,r),ha="center",va="bottom",fontsize=6.5,color=MUT)
for p,y in ys["phylum"].items(): ax.text(0.04,y+0.16,p[0],fontsize=6.2,color=MUT)
txt=["Genomes and origins are counted separately (origins: NCBI taxonomy, minimum changes; genomes: confirmatory species).",
"C33 + UAA=Stop:   supported 9 genomes, 4 fish lineages (4 origins)   |   unconfirmed 4 genomes, 4 origins   |   dropped 1 genome (Globodera), 1 origin.",
"G33 + UAA=Stop:   supported 2 species in 1 genus (Craniella), 1 origin   |   dropped 1 genome (Bufotes), 1 origin.     All 17 non-U33: 17 genomes, 11 origins.",
"UAA=Tyr: calibration panel, about 5 origins (not new).   Non-panel: 1 genome (Rotylenchulus), blocked by the locked stem count (REPORT 3.6); P. goodeyi is supported weak only.",
"The tree is NCBI taxonomy (phylum to genus), not a molecular phylogeny, so origin counts are approximate.","Branch lengths have no meaning. Status words come from post-hoc audits."]
for i,s in enumerate(txt): ax.text(-0.1,-1.3-i*0.8,s,fontsize=6.6,color=INK if i in(1,2,3) else MUT)
ax.axis("off"); ax.set_title("Non-U33 tRNA-Tyr with UAA=Stop, and UAA=Tyr lineages",loc="left",fontsize=10,color=INK); save(fig,"fig1_tree")
# ---------- fig2 ----------
H="EGGVHRVQRVPKTEKQGRVHTSTMTVA"; G={}
for r in csv.DictReader(open("results/gate4/loop.tsv"),delimiter="\t"):
    if r["covered"]=="True" and r["label"]=="MTRF1L": G.setdefault(r["source"],r["loop"])
A=[("ref","Homo sapiens MTRF1L",H,"reference","ref"),
 ("tyr","Anomia chinensis",None,"EGGVHRFQRETRVKQARNFFSCTITIA","Gate 2"),("tyr","Anomia ephippium",None,"EGGVHRFQRETRVKKARTFYSCTITIA","Gate 2"),("tyr","Enigmonia aenigmatica",None,"EGGVQRFQRETRVKVPRVFFSCTITVA","Gate 2"),("tyr","Placuna vitream",None,"EGGVQRFQRETRVKSSRTFFSCTITIA","Gate 2"),
 ("tyr","Cephalodiscus hodgsoni (reads, partial)",None,"HKVQRIPKTENKNRIHSST","Gate 2"),("tyr","Cephalodiscus sp. (reads, partial)",None,"HKVQRIPKTERNNRIHSST","Gate 2"),
 ("tyr","Plicatula muricata",None,"EGGVHRVQRKPVTDRNNKRQTSTITVA","Gate 2"),("tyr","Plicatula plicata (reads)",None,"EGGVHRVQRKPVTDRNNKRQTSTITVA","Gate 4"),
 ("tyr","Radopholus similis (partial; manual splice)",None,"KRLPVTNSKGKGVLHTST","Gate 2"),
 ("tyrS","Rotylenchulus reniformis [S]",None,G.get("GEN_R_reniformis"),"Gate 4"),
 ("stop","Pratylenchus vulnus",None,G.get("GEN_P_vulnus"),"Gate 4"),("stop","Meloidogyne hapla",None,G.get("PROT_M_hapla"),"Gate 4"),("stopq","Rotylenchus robustus (no mitogenome call)",None,G.get("GEN_R_robustus"),"Gate 4"),
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
n=len(R); W=len(h); fig,ax=plt.subplots(figsize=(10.4,6.6)); ax.set_xlim(-13.6,W+0.6); ax.set_ylim(-3.9,n+2.2); lab={"tyr":"UAA=Tyr (Gate 2, calibration)","tyrS":"UAA=Tyr, sensitivity label [S]","stop":"UAA=Stop sisters","stopq":"UAA=Stop sisters","tb":"[E] Tier B: non-U33 Stop taxa","tbs":"[E] Tier B: U33 sisters","ref":""}
for c203 in(203,216):
    x=pos[c203]; ax.add_patch(plt.Rectangle((x-0.5,-0.7),1,n+0.9,color="#f0efec",zorder=0)); ax.text(x,n+0.55,str(c203),ha="center",fontsize=7.5,weight="bold",color=INK)
for k in(196,222): ax.text(pos[k],n+0.55,str(k),ha="center",fontsize=6.5,color=MUT)
prev=None; y=n-1; seps=[]
for i,(g,nm,s,src) in enumerate(R):
    yy=n-1-i; grp=lab[g]
    if prev is not None and grp!=prev and grp: seps.append(yy+0.5)
    prev=grp if grp else prev
    ax.text(-0.6,yy,nm,ha="right",va="center",fontsize=7,style="italic" if i else "normal",color=INK)
    for j,ch in enumerate(AL[i]):
        col=INK; bg=None; w="normal"
        if ch!="-" and j==pos[216] and ch=="S": bg=OR; w="bold"
        if ch!="-" and j==pos[216] and ch=="T": bg="#dff3ea"
        if ch!="-" and j==pos[203] and ch!="Q": bg="#cfe0f6"; w="bold"
        if ch!="-" and j==pos[203] and ch=="Q": bg=None
        if bg: ax.add_patch(plt.Rectangle((j-0.46,yy-0.46),0.92,0.92,color=bg,zorder=1))
        ax.text(j,yy,ch,ha="center",va="center",family="DejaVu Sans Mono",fontsize=7.6,color=col,weight=w,zorder=2)
    if i in(0,): pass
for s in seps: ax.plot([-12.0,W],[s,s],color=MUT,lw=0.5)
def block(y0,y1,t): ax.text(-11.3,(y0+y1)/2,t,ha="left",va="center",fontsize=6.6,color=MUT)
yb=lambda g0:[n-1-i for i,(g,*_) in enumerate(R) if lab[g]==g0]
for t in dict.fromkeys(lab.values()):
    v=yb(t)
    if t and v: pass
txt=["Orange cell: S at human 216 (T216S). Green cell: T at 216. Blue cell: residue other than Q at 203. Human numbering NP_061914.3; insertions relative to human take extra columns (Radopholus +2 aa; its placement by MAFFT is approximate).",
"Partial sequences (Cephalodiscus, Radopholus) start or end inside 196–222; gaps there are missing data. [S] sensitivity: Rotylenchulus is blocked by the trnY stem audit. [E] exploratory block: Albula, Lycenchelys, Craniella (non-U33 Stop) and three sisters.",
"H. columbus is omitted (bacterial contaminant). Rotylenchus robustus has no mitogenome and no UAA call; its Stop status is inferred."]
for i,s in enumerate(txt): ax.text(-13.4,-1.5-i*0.75,s,fontsize=6.4,color=MUT)
for t,ys_ in [(t,[n-1-i for i,(g,*_) in enumerate(R) if lab[g]==t]) for t in dict.fromkeys(lab.values()) if t]:
    if ys_: ax.text(-11.3,max(ys_)+0.52,t,fontsize=7,weight="bold",color=INK,va="bottom") if False else None
ax.axis("off"); ax.set_title("mtRF1a loop (human 196–222): UAA=Tyr lineages, Stop sisters, and Tier B non-U33 Stop taxa",loc="left",fontsize=10,color=INK)
# block labels at the far left as coloured bars
bars=[("ref",""),]
for (gset,tx,cl) in ((("tyr","tyrS"),"UAA = Tyr",AQ),(("stop","stopq"),"UAA = Stop",BL),(("tb","tbs"),"[E] exploratory (Tier B)",OR)):
    v=[n-1-i for i,(g,*_) in enumerate(R) if g in gset]; ax.add_patch(plt.Rectangle((-13.3,min(v)-0.45),0.16,max(v)-min(v)+0.9,color=cl)); ax.text(-12.55,(max(v)+min(v))/2,tx,rotation=90,va="center",ha="center",fontsize=7,color=INK)
save(fig,"fig2_mtrf1a_loop"); print("fig2 alignment columns",W,"pos203",pos[203],"pos216",pos[216])
# ---------- fig3 ----------
def wil(k,n,z=1.96): p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; e=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-e,c+e
fr=[5,10,25,50,100]; cnt=[(0,325),(79,325),(291,325),(310,325),(63,65)]; cod=[(0,65),(16,65),(59,65),(63,65),(63,65)]
fig,(a1,a2)=plt.subplots(1,2,figsize=(10,4),gridspec_kw=dict(width_ratios=[1.15,1])); X=range(5)
for d,lb,cl,mk,off in((cnt,"Counter rule alone",BL,"o",-0.05),(cod,"With Codetta agreeing (replicate 0)",OR,"s",0.05)):
    y=[k/n for k,n in d]; lo=[y[i]-wil(*d[i])[0] for i in X]; hi=[wil(*d[i])[1]-y[i] for i in X]; a1.errorbar([x+off for x in X],y,yerr=[lo,hi],color=cl,marker=mk,ms=6,lw=1.6,capsize=2,label=lb)
for x,(k,n) in enumerate(cnt): a1.annotate(f"{k}/{n}",(x-0.05,k/n),textcoords="offset points",xytext=[(0,-13),(-19,11),(-21,2),(-14,10),(-3,9)][x],ha="center",fontsize=6.8,color=MUT)
a1.set_xticks(list(X)); a1.set_xticklabels([f"{f}%" for f in fr]); a1.set_ylim(-0.12,1.08); a1.set_ylabel("Fraction of genomes called Tyr"); a1.set_xlabel("Conserved Tyr codons replaced by TAA"); a1.legend(frameon=False,loc="lower right",fontsize=7); a1.set_title("Planted-null detection (Wilson 95% intervals)",loc="left",fontsize=9,color=INK)
LO=list(csv.DictReader(open("results/gate2/lolo.tsv"),delimiter="\t")); S=collections.OrderedDict(); pr=["Fixed"]
for r in LO:
    if r["unit"]=="lineage": S.setdefault(r["held_out"],{})[r["proc"]]=(int(r["recalled"]),int(r["genomes"]))
procs=list(dict.fromkeys(p for v in S.values() for p in v)); ln=list(S); print("LOLO procedures",procs)
for i,l in enumerate(ln):
    for j,p in enumerate(procs):
        k,n=S[l].get(p,(None,None))
        if k is None: continue
        f=k/n; a2.add_patch(plt.Rectangle((j-0.45,len(ln)-1-i-0.4),0.9,0.8,facecolor="#2a78d6" if f==1 else "#9cc3ef" if f>=0.5 else "#dbe9fa",edgecolor="white")); a2.text(j,len(ln)-1-i,f"{k}/{n}",ha="center",va="center",fontsize=7.5,color="white" if f==1 else INK)
a2.set_xlim(-0.6,len(procs)-0.4); a2.set_ylim(-0.7,len(ln)-0.3); a2.set_xticks(range(len(procs))); a2.set_xticklabels([{"Fixed":"Fixed rule\n(in-sample)","N":"Negatives-only\nthresholds","M":"Midpoint\nthresholds","R":"0.8 × minimum\nthresholds"}.get(p,p) for p in procs],fontsize=6.8); a2.set_yticks(range(len(ln))); a2.set_yticklabels(ln[::-1],fontsize=7,style="italic"); a2.tick_params(length=0); a2.set_title("Leave-one-lineage-out recall (genomes)",loc="left",fontsize=9,color=INK)
for sp in a2.spines.values(): sp.set_visible(False)
fig.text(0.01,-0.03,"The fixed rule was set by hand with all lineages in view, so its column is in-sample; the other three set thresholds without the held-out lineage. False Tyr: 0 of 340 panel negatives. Sources: GATE_2 §1 and A6, results/gate2/lolo.tsv.",fontsize=6.8,color=MUT)
fig.tight_layout(); save(fig,"fig3_sensitivity"); print("figures:",sorted(os.listdir(O)))
