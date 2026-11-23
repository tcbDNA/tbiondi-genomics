"""Fetch SRA files from the NCBI open-data bucket, stream FASTA through diamond, keep hits only. usage: mtrf_sra.py RUN..."""
import sys,subprocess,os,shlex,time
B="env/prefix/bin/"; R="data/mtrf/reads/"; DB="results/work/rf1_diamond"
for run in sys.argv[1:]:
    out=R+run+".hits.tsv"
    if os.path.exists(out): print(run,"exists",flush=True); continue
    f=R+run+".sra"; t0=time.time()
    r=subprocess.run(["curl","-sSL","--fail","--retry","3","-C","-","--max-time","2400","-o",f,f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"],capture_output=True,text=True)
    if r.returncode: print("DOWNLOAD_FAIL",run,r.stderr[-150:],flush=True); continue
    h=subprocess.run(["curl","-sSI",f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"],capture_output=True,text=True).stdout.lower(); want=int(h.split("content-length:")[1].split()[0])
    if os.path.getsize(f)!=want: print("SIZE_MISMATCH",run,os.path.getsize(f),want,flush=True); continue
    v=subprocess.run([B+"vdb-validate",f],capture_output=True,text=True); print(run,"vdb-validate rc",v.returncode,flush=True)
    if v.returncode: continue
    print(run,"downloaded",os.path.getsize(f),round(time.time()-t0),"s",flush=True)
    cmd=f"{B}fastq-dump -Z --fasta 0 {f} 2>/dev/null | {B}diamond blastx -d {DB} -q - --sensitive -e 1e-3 -k 1 --threads 6 -b 0.4 -c 4 -f 6 qseqid sseqid pident length evalue qstart qend sstart send bitscore full_qseq --quiet > {out}.tmp"
    r=subprocess.run(["bash","-o","pipefail","-c",cmd],capture_output=True,text=True,timeout=3000)
    if r.returncode: print("DIAMOND_FAIL",run,r.stderr[-200:],flush=True); continue
    os.rename(out+".tmp",out); os.remove(f); print(run,"hits",sum(1 for _ in open(out)),round(time.time()-t0),"s",flush=True)
