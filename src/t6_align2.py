"""Second pass on reads already extracted (t6/LABEL.reads.fa): relaxed identity/length, counting by identity class. usage: t6_align2.py LABEL ACC G34 STRAND"""
import sys,os,subprocess,collections
label,acc,g34,strand=sys.argv[1:5]; g34=int(g34); B="env/prefix/bin/"; W="results/work/t6/"
mt="".join(l.strip() for l in open("data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna") if l[0]!=">").upper(); pos=g34+1 if strand=="-" else g34-1
w0=max(1,pos-2500); w1=min(len(mt),pos+2500); open(W+label+".win.fa","w").write(f">w\n{mt[w0-1:w1]}\n"); site=pos-w0+1
subprocess.run([B+"makeblastdb","-in",W+label+".reads.fa","-dbtype","nucl","-out",W+label],capture_output=True)
r=subprocess.run([B+"blastn","-task","blastn","-query",W+label+".win.fa","-db",W+label,"-evalue","1e-10","-max_target_seqs","500000","-num_threads","6","-outfmt","6 sseqid pident length qstart qend sstart send qseq sseq"],capture_output=True,text=True).stdout
for mid,ml in((98,1000),(95,500),(90,300),(85,200)):
    c=collections.Counter(); seen=set()
    for l in sorted((x.split("\t") for x in r.splitlines()),key=lambda x:-int(x[2])):
        s,pi,ln,qs,qe=l[0],float(l[1]),int(l[2]),int(l[3]),int(l[4])
        if ln<ml or pi<mid or s in seen or not(qs<=site-10 and qe>=site+10): continue
        seen.add(s); p=qs; b="-"
        for x,y in zip(l[7],l[8]):
            if x!="-":
                if p==site: b=y; break
                p+=1
        c[{"A":"U","C":"G","G":"C","T":"A"}.get(b,b) if strand=="-" else b.replace("T","U")]+=1
    print(f"  identity>={mid}% length>={ml}: reads {sum(c.values())}, tRNA base {dict(c)}")
