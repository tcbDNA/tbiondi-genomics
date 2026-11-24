"""Gate 4 reads: RF1-type reads from a run (fasta from SRA file) by diamond blastx against the Gate 2 RF1 protein set (E<=1e-8), then translate the loop region with the Gate 2 alignment (human 196-222, positions 203/216), as mtrf_reads + mtrf_reads_loop.
Reads whose loop equals human are excluded as contamination; >=3 distinct reads are required. usage: gate4_reads.py RUN SRAFILE LABEL"""
import sys,subprocess,os,re,collections
sys.path.insert(0,"src"); import mtrf_guard
run,sra,label=sys.argv[1:4]; mtrf_guard.check([run]); B="env/prefix/bin/"; W="results/work/g4/"; os.makedirs(W,exist_ok=True)
DB="results/work/rf1_diamond"
if not os.path.exists(DB+".dmnd"): subprocess.run([B+"diamond","makedb","--in","data/mtrf/loop_in.faa","-d",DB],check=True,capture_output=True)
hits=W+run+".hits.tsv"
if not os.path.exists(hits):
    cmd=f"{B}fastq-dump -Z --fasta 0 {sra} 2>/dev/null | {B}diamond blastx -d {DB} --query - -e 1e-8 -k 1 -p 6 -f 6 qseqid sseqid pident length evalue qstart qend sstart send bitscore full_qseq --ultra-sensitive -o {hits}"
    r=subprocess.run(["bash","-c",cmd],capture_output=True,text=True); print(r.stderr[-300:],flush=True)
A={};cur=None
for l in open("data/mtrf/loop.aln"):
    if l[0]==">": cur=int(l[1:]); A[cur]=""
    else: A[cur]+=l.strip()
h=A[0]; hp={};p=0
for i,ch in enumerate(h):
    if ch!="-": p+=1; hp[i]=p
HUM="EGGVHRVQRVPKTEKQGRVHTSTMTVA"
T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; Bs="TCAG"; CODE={a+b+c:T[16*Bs.index(a)+4*Bs.index(b)+Bs.index(c)] for a in Bs for b in Bs for c in Bs}
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
rows=[]; n=0
for l in open(hits):
    f=l.rstrip("\n").split("\t"); sid,ev,qs,qe,ss,se,seq=int(f[1]),float(f[4]),int(f[5]),int(f[6]),int(f[7]),int(f[8]),f[10]; n+=1
    a=A.get(sid); 
    if a is None: continue
    m={};q=0
    for i,ch in enumerate(a):
        if ch!="-": q+=1; m[q]=hp.get(i)
    hs=[m.get(x) for x in range(ss,se+1) if m.get(x)]
    if not hs or min(hs)>205 or max(hs)<220: continue
    s=seq if qs<=qe else rc(seq); st=(min(qs,qe)-1) if qs<=qe else len(seq)-max(qs,qe); pep="".join(CODE.get(s[j:j+3],"X") for j in range(st,len(s)-2,3))
    # residues at human 203 / 216 by position in the aligned HSP (no gaps assumed inside 200-220)
    first=m.get(ss); rows.append((f[0],sid,f[2],ev,first,pep[:70]))
print(label,run,"RF1-type reads",n,"spanning human 205-220:",len(rows))
c=collections.Counter()
for r in rows:
    print("   ",r[0][-18:],"human start",r[4],"id",r[2],"e",r[3],r[5]); c[r[5][:24]]+=1
