"""Exact-k-mer read support with the correct allele mapping, streaming a run through ripgrep (single-line FASTA). usage as a module: support(run, patterns{label:31mer list}, cap_spots)
Both strands are searched. Counts are read-level hits (a read counted once per label)."""
import sys,os,subprocess,collections,re
B="env/prefix/bin/"
RC=str.maketrans("ACGTN","TGCAN")
def rc(s): return s.translate(RC)[::-1]
def fetch(run,maxgb=20):
    f=f"results/work/{run}.sra"
    if os.path.exists(f): return f
    sz=subprocess.run(["curl","-sSI",f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"],capture_output=True,text=True).stdout.lower(); m=re.search(r"content-length:\s*(\d+)",sz); gb=int(m[1])/1e9 if m else 0
    if gb>maxgb: return None
    subprocess.run(["sh","src/dl.sh",f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}",f,"3"],stdout=subprocess.DEVNULL); return f if os.path.exists(f) else None
def support(run,labels,cap=0,maxgb=20,keep=False):
    f=fetch(run,maxgb)
    if not f: return None,"not downloaded (size cap or error)"
    pats=collections.defaultdict(set)
    for lab,ks in labels.items():
        for k in ks: pats[k].add(lab); pats[rc(k)].add(lab)
    pf=f"results/work/{run}.pat"; open(pf,"w").write("\n".join(pats)+"\n")
    cmd=f"{B}fastq-dump -Z --fasta 0 {'-X %d'%cap if cap else ''} {f} 2>/dev/null | env/tools/kmatch/kmatch {pf}"
    r=subprocess.run(["bash","-c",cmd],capture_output=True,text=True); cnt=collections.Counter()
    for l in r.stdout.splitlines():
        n,k=l.split()
        for lab in pats.get(k,()): cnt[lab]+=int(n)
    if not keep: os.remove(f); os.remove(pf)
    return cnt,""
def variants(seq,site,strand,length=31):
    """site: 0-based + strand index of N33; returns (record kmer, kmer with tRNA U33)."""
    h=length//2; a=seq[site-h:site+h+1]
    if len(a)!=length: return None,None
    base_plus_for_tU="A" if strand=="-" else "T"
    b=seq[site-h:site]+base_plus_for_tU+seq[site+1:site+h+1]
    return a,b
