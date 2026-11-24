"""EXPLORATORY: TAA in GenBank COX1/COI fragments of Tylenchomorpha (all family), as done for Rotylenchulus. Every fragment is searched in 6 frames against the COX1 profile only (same profile and table-5 translation with TAA/TAG=X as uaa_caller);
records the codon at each conserved-Tyr column (p(Tyr)>=0.5), plus TAA/TAG anywhere in the aligned frame, and /transl_table. Terminal 10 codons of the fragment masked (partial fragment ends). usage: tylench_cox1.py"""
import sys,os,re,collections,json,csv,math,time; sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import ncbi, uaa_caller as U, pyhmmer, numpy as np
CAP=int(sys.argv[1]) if len(sys.argv)>1 else 12
ABC=U.ABC; H=[h for h in U.hmms() if str(h.name)=="COX1"][0]; ME=np.array(H.match_emissions); Yi=ABC.symbols.index("Y")
Q='"Tylenchomorpha"[Organism] AND (COX1[Gene] OR COI[Gene] OR CO1[Gene] OR "cytochrome c oxidase subunit I"[Title] OR "cytochrome oxidase subunit I"[Title] OR "cytochrome c oxidase subunit 1"[Title]) AND 250:2500[SLEN] AND mitochondrion[filter]'
os.makedirs("data/tylench_cox1",exist_ok=True); n,ids,_=ncbi.search("nuccore",Q,100000); print("records",n,"ids",len(ids),flush=True)
rows=[]; OUT=open("results/gate3/tylench_cox1.tsv","w"); OUT.write("acc\tspecies\tfamily\tgenus\tlen\ttransl_table\tcols\ttaa_tyr\ttat_tyr\ttac_tyr\ttag_tyr\tttt_tyr\tother_tyr\ttaa_any\ttag_any\tdomscore\tauthor\n")
per=collections.Counter()
def lineage(b):
    m=re.search(r"ORGANISM\s+(.*?)\nREFERENCE",b,re.S); t=m[1] if m else ""; sp=t.split("\n")[0].strip(); L=re.sub(r"\s+"," "," ".join(t.split("\n")[1:])).rstrip(".").split("; "); return sp,L
for i in range(0,len(ids),200):
    f=f"data/tylench_cox1/b{i//200:04d}.gb"
    if not os.path.exists(f):
        for k in range(5):
            try: t=ncbi.call("efetch",db="nuccore",id=",".join(ids[i:i+200]),rettype="gb",retmode="text").decode(errors="replace"); break
            except Exception as e: time.sleep(5)
        open(f+".tmp","w").write(t); os.replace(f+".tmp",f)
    t=open(f).read(); blk=[b for b in t.split("\n//") if b.strip()]; names=[]; seqs=[]; meta={}
    for b in blk:
        m=re.search(r"^ACCESSION\s+(\S+)",b,re.M); s=re.sub(r"[^a-z]","",b.split("ORIGIN",1)[1].lower()).upper() if "ORIGIN" in b else ""
        if not m or not s: continue
        sp,L=lineage(b); fam=next((x for x in L if x.endswith("idae")),""); tt=re.findall(r"/transl_table=(\d+)",b); au=re.search(r"AUTHORS\s+(\S+)",b)
        meta[m[1]]=(sp,fam,sp.split(" ")[0],len(s),tt[0] if tt else "?",au[1] if au else ""); seqs.append((m[1],re.sub(r"[^ACGT]","N",s)))
    fr=[]; own=[]
    for a,s in seqs:
        for st,x in (("+",s),("-",U.rc(s))):
            for p in range(3):
                cod=[x[j:j+3] for j in range(p,len(x)-2,3)]
                aa="".join(U.TAB.get(c,"X") if set(c)<=set("ACGT") else "X" for c in cod)
                fr.append(pyhmmer.easel.TextSequence(name=f"{len(own)}".encode() if False else f"{len(own)}",sequence=aa).digitize(ABC)); own.append((a,cod))
    if not fr: continue
    best={}
    for res in pyhmmer.hmmsearch(H,fr,cpus=2,E=1e-5,domE=1e-5,incE=1e-5,incdomE=1e-5,bias_filter=False,F1=1.0,F2=1.0,F3=1.0,Z=1,domZ=1):
        for hit in res:
            a,cod=own[int(str(hit.name))]
            for d in hit.domains:
                if d.i_evalue>1e-5: continue
                if a not in best or d.score>best[a][0]: best[a]=(d.score,d.alignment,cod)
    for a,(sc,al,cod) in best.items():
        hp,tp=al.hmm_from-1,al.target_from-1; c=collections.Counter(); anyT=anyG=0; cols=[]
        for hc,tc in zip(al.hmm_sequence,al.target_sequence):
            if hc!=".": hp+=1
            if tc!="-": tp+=1
            if hc!="." and tc!="-" and hc.isalpha():
                cd=cod[tp-1]; 
                if tp<=10 or tp>len(cod)-10 or hp<=10 or hp>H.M-10: continue
                if cd=="TAA": anyT+=1
                elif cd=="TAG": anyG+=1
                if ME[hp,Yi]>=0.5: cols.append(hp); c[cd if cd in("TAA","TAT","TAC","TAG","TTT") else "other"]+=1
        sp,fam,gen,ln,tt,au=meta[a]; rows.append((a,sp,fam,gen,ln,tt,len(cols),c["TAA"],c["TAT"],c["TAC"],c["TAG"],c["TTT"],c["other"],anyT,anyG,round(sc,1),au)); OUT.write("\t".join(map(str,rows[-1]))+"\n")
    if (i//200)%10==0: print(i+200,"of",len(ids),"aligned",len(rows),flush=True)
OUT.close(); print("aligned records",len(rows),flush=True); print("DONE",flush=True)
