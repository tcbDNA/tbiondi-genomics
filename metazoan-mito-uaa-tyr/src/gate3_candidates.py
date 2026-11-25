"""Candidate checks for sufficiency-counterexample candidates (T6 and the user's list: Albula vulpes, Microdevario kubotai, Bufotes pewzowi). Written before the sweep finished.
usage: gate3_candidates.py ACC [ACC...]    -> results/gate3/candidate_<ACC>.txt ; each check is reported separately and left to the reader:
 1 record: circularity, 13 intact PCG families, internal TAA/TAG, second GUA trnY loci (from the sweep line)
 2 conspecific/congeneric mitochondrial records in nuccore (complete and partial >=500 bp): trnY N33 by the same extractor, both callers where the call is non-U
 3 read support: SRA runs of the same species; exact 31-mers centred on N33 for the record allele and the U allele, both strands, streamed (cap MAXSPOTS spots/run; only runs <= MAXGB GB)
 4 NUMT: nuclear assembly of the species (if any) searched with trnY +/-150 bp (blastn), and mitochondrial-depth note
Guards: needs GATE3_UNLOCK (PREREG committed). Network: NCBI only."""
import sys,os,re,json,csv,subprocess,tempfile,collections
sys.path.insert(0,"src"); import ncbi, panel_guard, trny_extract as T, gate3_run as G
MAXSPOTS=15000000; MAXGB=3.0; B="env/prefix/bin/"
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def load(p): return "".join(l.strip() for l in open(p) if l[0]!=">").upper()
def conspecifics(genus,species,W,limit=400):
    q=f'"{genus}"[Organism] AND mitochondrion[filter] AND 500:30000[SLEN]'; n,ids,_=ncbi.search("nuccore",q,limit); W(f"2 nuccore mitochondrial records for genus {genus}: {n} (first {len(ids)} examined)"); rows=[]
    os.makedirs("results/work/cons",exist_ok=True)
    for i in range(0,len(ids),50):
        t=ncbi.call("efetch",db="nuccore",id=",".join(ids[i:i+50]),rettype="fasta",retmode="text").decode(errors="replace")
        for blk in t.split(">")[1:]:
            h,_,s=blk.partition("\n"); a=h.split()[0]; s=re.sub(r"\s","",s).upper()
            p=f"results/work/cons/{a}.fna"; open(p,"w").write(f">{a}\n{s}\n"); d=T.call(p)
            if d.get("status")=="called": rows.append((a,h.split(" ",1)[-1][:70],len(s),d["N33"],d["loop"],d["score"],d["stem_pairs"],d["strand"],d["g34_pos"]))
    import collections; c=collections.Counter((r[3]) for r in rows); W("  with a verified GUA trnY call:",len(rows),"of",len(ids),"N33:",dict(c))
    for r in rows:
        if r[3]!="U" or species in r[1]: W("   ",r[0],r[1],"len",r[2],"N33",r[3],"loop",r[4],"score",r[5],"stem",r[6])
    return rows
def runs(species,W):
    n,ids,_=ncbi.search("sra",f'"{species}"[Organism]',200)
    if not ids: W("3 no SRA runs for",species); return []
    x=ncbi.call("efetch",db="sra",id=",".join(ids),rettype="runinfo",retmode="text").decode().split("\n"); h=x[0].split(","); out=[]
    for l in x[1:]:
        f=l.split(",")
        if len(f)>=len(h)-2 and f[0]: d=dict(zip(h,f)); out.append(d)
    W("3 SRA runs for",species,":",len(out),[(d["Run"],d["LibraryStrategy"],d["LibrarySource"],round(float(d["size_MB"] or 0)/1000,2)) for d in out][:12]); return out
def support(acc,seqp,d,species,W):
    s=load(seqp); g=int(d["g34_pos"]); strand=d["strand"]; site=g-2 if strand!="-" else g   # 0-based index of N33 in the + strand
    n33=s[site]; al=n33; alt={"C":"T","T":"C","A":"T","G":"T"}.get(n33,"T"); wt=s[:site]+alt+s[site+1:]
    k=lambda x: x[site-15:site+16]; alls={"record":k(s),"U/alt":k(wt)}; W("  31-mer at N33: record",alls["record"],"alt",alls["U/alt"],"| record base on + strand",n33)
    if len(alls["record"])<31: W("  N33 too close to the sequence end; read check skipped"); return
    for r in runs(species,W):
        gb=float(r["size_MB"] or 0)/1000
        if gb>MAXGB or r["LibrarySource"] not in("GENOMIC","METAGENOMIC","TRANSCRIPTOMIC","GENOMIC SINGLE CELL","TRANSCRIPTOMIC SINGLE CELL"): W("  run",r["Run"],"skipped (size",round(gb,1),"GB or source)"); continue
        f=f"results/work/{r['Run']}.sra"; subprocess.run(["sh","src/dl.sh",f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{r['Run']}/{r['Run']}",f,"3"])
        if not os.path.exists(f): W("  run",r["Run"],"not downloaded"); continue
        c=collections.Counter(); a,b=alls["record"],alls["U/alt"]; pat={a:"record",rc(a):"record",b:"alt",rc(b):"alt"}
        p=subprocess.Popen([B+"fastq-dump","-Z","--fasta","0","-X",str(MAXSPOTS),f],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
        for l in p.stdout:
            if l[0]!=">":
                for q,nm in pat.items():
                    if q in l.upper(): c[nm]+=1
        p.wait(); os.remove(f); W("  run",r["Run"],r["LibraryStrategy"],"reads carrying the exact 31-mer:",dict(c))
def main(acc):
    M={r["acc"]:r for r in csv.DictReader(open("results/gate3/sweep_audit.tsv"),delimiter="\t")}[acc]; path="data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna"; sp=M["species"]; gen=sp.split()[0]
    with open(f"results/gate3/candidate_{acc}.txt","w") as o:
        W=lambda *x: print(*x,file=o,flush=True); d=None
        for l in open("results/gate3/run/out.jsonl"):
            x=json.loads(l)
            if x["acc"]==acc: d=x
        W(acc,sp,M["phylum"],M["order"],"length",M["length"],"topology",M["topology"],"provenance",M["provenance"]); g={x[3] for x in d["pcg"] or []}
        W("1 PCG families hit",len(g),"internal TAA",d["uaa"]["n_taa_internal"],"internal TAG",d["uaa"]["n_tag_internal"],"GUA loci ARAGORN",len(d["second"]["aragorn"]),"tRNAscan",len(d["second"]["trnascan"]),"verified GUA arms",d["trny"].get("n_ok"),"| trnY",{k:d["trny"].get(k) for k in("score","N33","loop","stem_pairs","start","end","strand")})
        R=[r for r in csv.DictReader(open("results/refseq_mito_inventory.tsv"),delimiter="\t") if r["accession"]==acc]; W("   RefSeq derived_from:",R[0]["derived_from"] if R else "n/a")
        conspecifics(gen,sp,W); support(acc,path,d["trny"],sp,W)
        A=[a for a in csv.DictReader(open("results/assemblies_metazoa.tsv"),delimiter="\t") if a["organism"].startswith(sp)]; W("4 nuclear assemblies for",sp,":",[(a["accession"],a["level"]) for a in A])
        if A:
            z=f"results/work/{acc}_nuc.zip"; subprocess.run([B+"datasets","download","genome","accession",A[0]["accession"],"--include","genome","--filename",z,"--no-progressbar"],capture_output=True)
            if os.path.exists(z):
                subprocess.run(["sh","-c",f"unzip -p {z} '*genomic.fna' > results/work/{acc}_nuc.fna && rm {z}"]); s=load(path); t=d["trny"]; a,b=sorted((t["start"],t["end"])); open(f"results/work/{acc}_q.fa","w").write(f">q\n{s[max(0,a-151):b+150]}\n")
                subprocess.run([B+"makeblastdb","-in",f"results/work/{acc}_nuc.fna","-dbtype","nucl","-out",f"results/work/{acc}_db"],capture_output=True)
                r=subprocess.run([B+"blastn","-query",f"results/work/{acc}_q.fa","-db",f"results/work/{acc}_db","-evalue","1e-5","-outfmt","6 sseqid pident length qstart qend sstart send evalue","-num_threads","4"],capture_output=True,text=True); W("   trnY+/-150 bp vs nuclear assembly (NUMT candidates):",[l for l in r.stdout.splitlines()][:6]); os.remove(f"results/work/{acc}_nuc.fna")
        W("DONE")
if __name__=="__main__":
    panel_guard.check(["__gate3_candidates__"]) if os.environ.get("GATE3_UNLOCK")!="1" else panel_guard.check(["x"])
    for a in sys.argv[1:]: main(a)
