"""Firewall: callers refuse genomes that are not in the frozen control panel (until PREREG.md is committed and --unlock is given)."""
import csv,os,subprocess,hashlib,sys
PANEL="results/gate2_panel.tsv"
def allowed():
    return {r["panel_id"] for r in csv.DictReader(open(PANEL),delimiter="\t")}
def check(ids):
    if os.environ.get("GATE3_UNLOCK")=="1":
        r=subprocess.run(["git","log","-1","--format=%H","--","PREREG.md"],capture_output=True,text=True)
        if not r.stdout.strip(): sys.exit("GATE3_UNLOCK set but PREREG.md is not committed; refusing")
        if subprocess.run(["git","diff","--quiet","HEAD","--","PREREG.md"]).returncode: sys.exit("GATE3_UNLOCK: PREREG.md differs from the committed version; refusing")
        h=hashlib.sha256(open("PREREG.md","rb").read()).hexdigest()
        if h not in open("reports/GATE_2_CALIBRATION.md").read(): sys.exit("GATE3_UNLOCK: PREREG.md SHA256 is not the one recorded in GATE_2_CALIBRATION.md; refusing")
        return
    bad=[i for i in ids if i.split("|")[0] not in allowed() and not i.startswith("PLANTED_")]
    if bad: sys.exit("FIREWALL: not in control panel: "+", ".join(bad[:5]))
