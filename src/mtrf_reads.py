"""Resumable read retrieval + diamond blastx against RF1-type proteins. usage: mtrf_reads.py RUN...
Downloads to data/mtrf/reads (verified by byte count), keeps only diamond hits, then deletes FASTQ if --purge."""
import sys,subprocess,os,urllib.request,json
import sys as _s; _s.path.insert(0,"src"); import mtrf_guard
mtrf_guard.check(sys.argv[1:])
RUNS={"SRR1695471":"Cephalodiscus hodgsoni RNA-Seq","SRR11101524":"Cephalodiscus sp. RNA-Seq","SRR11101525":"Rhabdopleura annulata RNA-Seq","ERR15194044":"Pseudotanais sp. a WGS"}
B="env/prefix/bin/"; R="data/mtrf/reads/"; os.makedirs(R,exist_ok=True); DB="results/work/rf1_diamond"
if not os.path.exists(DB+".dmnd"): subprocess.run([B+"diamond","makedb","--in","data/mtrf/loop_in.faa","-d",DB],check=True,capture_output=True)
for run in sys.argv[1:]:
    out=R+run+".hits.tsv"
    if os.path.exists(out): print(run,"done already",flush=True); continue
    rec=json.loads(urllib.request.urlopen(f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run}&result=read_run&fields=fastq_ftp,fastq_bytes,fastq_md5&format=json",timeout=60).read())[0]
    urls=rec["fastq_ftp"].split(";"); sizes=[int(x) for x in rec["fastq_bytes"].split(";")]; print(run,RUNS[run],sizes,flush=True); tmp=out+".tmp"; open(tmp,"w").close(); ok=True
    for i,(u,sz) in enumerate(zip(urls,sizes)):
        f=R+os.path.basename(u)
        for a in range(12):
            if os.path.exists(f) and os.path.getsize(f)==sz: break
            subprocess.run(["curl","-sSL","--fail","--retry","3","-C","-","--max-time","1500","-o",f,"https://"+u],capture_output=True)
        if not(os.path.exists(f) and os.path.getsize(f)==sz): print("DOWNLOAD_INCOMPLETE",run,u,os.path.getsize(f) if os.path.exists(f) else 0,flush=True); ok=False; break
        cmd=f"gunzip -c {f} | awk 'NR%4==1{{split(substr($1,2),a,\" \"); h=a[1]}} NR%4==2{{print \">\"h\"_{i}\"; print}}' | {B}diamond blastx -d {DB} -q - --sensitive -e 1e-3 -k 1 --threads 8 -f 6 qseqid sseqid pident length evalue qstart qend sstart send bitscore full_qseq --quiet >> {tmp}"
        r=subprocess.run(["bash","-o","pipefail","-c",cmd],capture_output=True,text=True)
        if r.returncode: print("DIAMOND_FAIL",run,r.stderr[-200:],flush=True); ok=False; break
        os.remove(f)
    if ok: os.rename(tmp,out); print(run,"hits",sum(1 for _ in open(out)),flush=True)
