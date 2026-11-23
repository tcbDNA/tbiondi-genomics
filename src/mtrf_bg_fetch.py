"""Download annotated protein sets for panel species (one per genus) via NCBI datasets. Logs failures."""
import json,subprocess,os,zipfile,glob
import sys; sys.path.insert(0,"src"); import panel_guard
C=json.load(open("results/work/mtrf_bg_candidates.json")); panel_guard.check([c[0] for c in C]); seen=set(); sel=[]
for pid,sp,lin,acc,lv in C:
    g=sp.split()[0]
    if g in seen: continue
    seen.add(g); sel.append((pid,sp,lin,acc))
os.makedirs("data/mtrf/bg",exist_ok=True); log=open("results/audit/mtrf_bg_fetch.tsv","w"); log.write("panel_id\tspecies\taccession\tstatus\n")
for pid,sp,lin,acc in sel:
    d=f"data/mtrf/bg/{acc}.faa"
    if not os.path.exists(d):
        z=f"data/mtrf/bg/{acc}.zip"; subprocess.run(["env/prefix/bin/datasets","download","genome","accession",acc,"--include","protein","--filename",z,"--no-progressbar"],capture_output=True)
        try:
            with zipfile.ZipFile(z) as zf:
                n=[x for x in zf.namelist() if x.endswith("protein.faa")]
                if n: open(d,"wb").write(zf.read(n[0]))
        except Exception as e: log.write(f"{pid}\t{sp}\t{acc}\tFAIL {e}\n"); continue
        finally:
            if os.path.exists(z): os.remove(z)
    log.write(f"{pid}\t{sp}\t{acc}\t{'ok' if os.path.exists(d) else 'no_protein_file'}\n"); log.flush()
print(len(sel))
