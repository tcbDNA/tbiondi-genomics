"""Firewall: callers refuse genomes that are not in the frozen control panel (until PREREG.md is committed and --unlock is given)."""
import csv,os,subprocess,hashlib,sys
PANEL="results/gate2_panel.tsv"
def allowed():
    return {r["panel_id"] for r in csv.DictReader(open(PANEL),delimiter="\t")}
def check(ids):
    if os.environ.get("GATE3_UNLOCK")=="1":
        r=subprocess.run(["git","log","-1","--format=%H","--","PREREG.md"],capture_output=True,text=True)
        if r.stdout.strip(): return
        sys.exit("GATE3_UNLOCK set but PREREG.md is not committed; refusing")
    bad=[i for i in ids if i.split("|")[0] not in allowed() and not i.startswith("PLANTED_")]
    if bad: sys.exit("FIREWALL: not in control panel: "+", ".join(bad[:5]))
