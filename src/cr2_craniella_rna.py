"""Cold review 2, check 2 part 4 (POST-HOC): do the RNA-seq G33 reads come from the mitochondrion or a nuclear copy? Reads that span N33 as tRNA-Tyr(GUA) in check 4 (windows with up to 110 nt of flank each side) are aligned by blastn to the species' own assembly; the best bit score on the mitogenome is compared with the best on any nuclear contig.
(Whole-run mapping of the RNA libraries to the nuclear assembly was started and stopped as too slow; the reads that matter are all here.) writes results/posthoc2/craniella_rna.txt"""
import sys,subprocess,collections
sys.path.insert(0,"src"); import cr2_secondcopy as S
from concurrent.futures import ProcessPoolExecutor
o=open("results/posthoc2/craniella_rna.txt","w"); P=lambda *x:(print(*x),print(*x,file=o)); P("POST-HOC check 2, part 4")
if __name__=="__main__":
    with ProcessPoolExecutor(8) as ex:
        for acc,run,asm in (("OZ318138.1","SRR33975184","GCA_974939545.1"),("OZ318138.1","ERR14792883","GCA_974939545.1"),("OZ407450.1","ERR14857142","GCA_978397455.1"),("OZ318138.1","ERR14835550","GCA_974939545.1"),("OZ407450.1","ERR14842678","GCA_978397455.1")):
            J,u,D=S.classify(acc,run,ex=ex); ok=[(h,s) for h,s in J if D[s] and D[s]["cover"] and D[s]["n"]>=40 and D[s]["ident"]>=.8 and D[s]["win"] and D[s]["contig"] and D[s]["ac"]=="GUA"]; fa=f"results/work/r2/cran/{acc}.{run}.n33.fa"; open(fa,"w").write("".join(f">r{i}\n{s}\n" for i,(h,s) in enumerate(ok)))
            r=subprocess.run(["env/prefix/bin/blastn","-task","megablast","-word_size","16","-query",fa,"-db",f"results/work/numt/{asm}","-evalue","1e-5","-dust","no","-max_target_seqs","500","-num_threads","6","-outfmt","6 qseqid sseqid bitscore pident length"],capture_output=True,text=True).stdout; m=collections.defaultdict(float); n=collections.defaultdict(float)
            for l in r.splitlines():
                q,c,b,*_=l.split("\t"); d=m if c==acc else n; d[q]=max(d[q],float(b))
            c=collections.Counter()
            for i,(h,s) in enumerate(ok): q=f"r{i}"; k="mitogenome better" if m[q]>n[q]+1 else "nuclear better" if n[q]>m[q]+1 else "tie"; c[(k,D[s]["b33"])]+=1
            P(acc,run,"reads spanning N33",len(ok),{f"{k[0]}, base {k[1]}":v for k,v in c.most_common()})
