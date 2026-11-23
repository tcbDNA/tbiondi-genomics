"""Firewall for mtRF1a scripts. Every data source must be registered here against a panel species, or be a logged deviation.
After PREREG.md is committed, GATE3_UNLOCK=1 lifts the restriction (same rule as panel_guard)."""
import csv,os,subprocess,sys
REF={"HUMAN","ECOLI_RF1","ECOLI_RF2"}
# source key -> panel species it stands for
SRC={"FIG_Ach":"Anomia chinensis","FIG_Hsq":"Heteranomia squamula","FIG_Pmu":"Plicatula muricata","FIG_Pvi":"Placuna vitream",
 "GENOME_AnoEph":"Anomia ephippium","GENOME_EniAen":"Enigmonia aenigmatica","GENOME_PodMac":"Pododesmus macrochisma","GENOME_RadSim":"Radopholus similis",
 "SRR1695471":"Cephalodiscus hodgsoni","ERR15194044":"Pseudotanais sp."}
# used in Gate 2 although the species is NOT in the frozen panel: logged in GATE_2_CALIBRATION.md section 7, closed to reuse
DEVIATIONS={"SRR11101524":"Cephalodiscus sp. (SRR11101524)","SRR11101525":"Rhabdopleura annulata (SRR11101525)"}
def _unlocked():
    if os.environ.get("GATE3_UNLOCK")!="1": return False
    return bool(subprocess.run(["git","log","-1","--format=%H","--","PREREG.md"],capture_output=True,text=True).stdout.strip())
def panel_species():
    P=list(csv.DictReader(open("results/gate2_panel.tsv"),delimiter="\t")); return {" ".join(r["species"].split()[:2]).lower() for r in P}, {r["panel_id"] for r in P}
def check(sources):
    if _unlocked(): return
    sp,ids=panel_species(); bg={}
    if os.path.exists("results/audit/mtrf_bg_fetch.tsv"):
        for l in open("results/audit/mtrf_bg_fetch.tsv"):
            f=l.rstrip("\n").split("\t")
            if len(f)==4 and f[0] in ids: bg["BG_"+f[2]]=f[0]
    bad=[]
    for s in set(sources):
        if s in REF or s.startswith("HUMAN") or s in bg: continue
        if s in SRC and " ".join(SRC[s].split()[:2]).lower() in sp: continue
        bad.append(s+(" (logged deviation; not reusable)" if s in DEVIATIONS else ""))
    if bad: sys.exit("FIREWALL mtRF1a: not a panel source: "+", ".join(sorted(bad)))
