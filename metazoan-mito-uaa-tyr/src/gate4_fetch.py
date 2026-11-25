"""Gate 4 downloads: protein sets and assemblies via NCBI datasets (zip hashed, then unpacked); TSA FASTA from NCBI FTP. 30 GB floor; everything ledgered."""
import subprocess,os,hashlib,shutil,sys,json,gzip
B="env/prefix/bin/"; G="data/mtrf/g4"; os.makedirs(G,exist_ok=True)
def free(): return shutil.disk_usage(".").free
def led(desc,f,n,sha): open("data/manifests/download_ledger.tsv","a").write(f"4\t{desc}\t{f}\t{n}\t{sha}\ton_disk\n")
def sh256(f):
    h=hashlib.sha256()
    for c in iter(lambda x=open(f,"rb"):x.read(1<<20),b""): h.update(c)
    return h.hexdigest()
PROT={"GCA_019788955.1":"Albula_glossodonta_BYU","GCF_013368585.1":"Megalops_cyprinoides","GCF_013347855.1":"Anguilla_anguilla","GCA_028022725.1":"Lycodopsis_pacificus","GCA_040110945.1":"Zoarces_viviparus","GCA_051171035.1":"Meloidogyne_hapla"}
GEN={"GCA_001026735.1":"Rotylenchulus_reniformis","GCA_054131405.1":"Pratylenchus_vulnus","GCA_016433075.1":"Hoplolaimus_columbus","GCA_034700465.1":"Rotylenchus_robustus","GCA_982936665.1":"Cinachyrella_apion"}
for d,kind in((PROT,"protein"),(GEN,"genome")):
    for a,n in d.items():
        f=f"{G}/{n}.{'faa' if kind=='protein' else 'fna'}"
        if os.path.exists(f): print("have",f,flush=True); continue
        if free()<32e9: print("STOP: disk"); sys.exit(2)
        z=f"{G}/{n}.zip"; subprocess.run([B+"datasets","download","genome","accession",a,"--include",kind if kind=="protein" else "genome","--filename",z,"--no-progressbar"],capture_output=True)
        if not os.path.exists(z): print("FAIL",a,n,flush=True); continue
        out=subprocess.run(["unzip","-p",z,"*.faa" if kind=="protein" else "*genomic.fna"],capture_output=True).stdout
        if len(out)<1000: print("EMPTY",a,n,kind,flush=True); os.remove(z); continue
        open(f,"wb").write(out); led(f"NCBI datasets {a} {kind}",f,os.path.getsize(z),sh256(z)); print("ok",a,n,kind,round(os.path.getsize(f)/1e6,1),"MB",flush=True); os.remove(z)
for m in("GGVO","GGVQ","GGVP","GGVR","GGVV"):
    f=f"{G}/tsa.{m}.fsa_nt.gz"
    if not os.path.exists(f):
        u=f"https://ftp.ncbi.nlm.nih.gov/genbank/tsa/G/tsa.{m}.1.fsa_nt.gz"; subprocess.run(["curl","-sSL","--fail","--retry","5","-o",f,u]); 
        if os.path.exists(f): led(u,f,os.path.getsize(f),sh256(f)); print("ok TSA",m,os.path.getsize(f),flush=True)
        else: print("FAIL TSA",m,flush=True)
print("FETCH DONE",flush=True)
