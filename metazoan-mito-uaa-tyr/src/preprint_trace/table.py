"""Table 1: look for each number only in lines about that row's species, accession or runs."""
import re,sys,csv,collections
import assign as A
doc=A.tok.doc; P=A.P; R='reports/REPORT.md'
rep=open(P+R,encoding='utf8').read().split('\n')
T=lambda p:list(csv.DictReader(open(P+p),delimiter='\t'))
i0=next(i for i,x in enumerate(doc) if x[2].startswith('Albula glossodonta (Albulidae)'))
rows=[]; cur=None
for i in range(i0,len(doc)):
    a,b,t=doc[i]
    if a!='Main Tables' and a!='Main Tables': break
    if b=='table cell 1': cur=[i]; rows.append(cur)
    elif b.startswith('table cell') and cur is not None: cur.append(i)
hc={}
for l in open(P+'results/posthoc2/secondcopy/highcopy_summary.txt'):
    m=re.match(r'(\w+): mitochondrial depth (\d+)x',l)
    if m: hc[m[1]]=l
def src_lines(row):
    """lines allowed for a Table 1 row"""
    head=doc[row[0]][2]; acc=re.search(r'([A-Z]{2}_?\d+\.\d)',head)[1]; sp=head.split(' (')[0]; g,e=sp.split()[:2]; root=acc.split('.')[0]
    runs=set(re.findall(r'[SED]RR\d+',' '.join(doc[k][2] for k in row))); out=[]
    for n,l in enumerate(rep,1):
        if root in l or re.search(r'\b'+e+r'\b',l) or any(x in l for x in runs): out.append((R+f':{n}',l))
    for f in ('read_power.tsv','read_power_check4.tsv','conspecific_specimen_audit.tsv','bufotes_framing.txt','conspecific_specimen_notes.txt','secondcopy/highcopy_summary.txt','secondcopy/highcopy_craniella.txt','secondcopy/contig_id_St_hypomelas.txt','craniella_depth2.txt','craniella_depth3.txt'):
        for n,l in enumerate(open(P+'results/posthoc2/'+f),1):
            if acc in l or root in l or re.search(r'\b'+e+r'\b',l) or any(x in l for x in runs) or (f.endswith('summary.txt') and re.search(e,l)): out.append((f'results/posthoc2/{f}:{n}',l))
    import glob
    for f in glob.glob(P+f'reports/dossiers/{acc}_*'):
        for n,l in enumerate(open(f),1): out.append((f.replace(P,'')+f':{n}',l))
    return out
def find(v,pct,lines,ctx):
    w=A.W(ctx)-A.STOPW; best=None
    for q in sorted(A.match.variants(v,pct)):
        for s,l in lines:
            nums=set(x.rstrip(',').replace(',','') for x in re.findall(r'(?<![A-Za-z_\d.])\d[\d,]*(?:\.\d+)?',l))
            hit=q in nums
            if not hit and pct and '.' in q:   # rounding of a 4-decimal share
                try: hit=any(abs(round(float(x)*100,len(q.split('.')[1]))-float(q))<1e-9 for x in nums if re.fullmatch(r'0\.\d{3,}',x))
                except ValueError: pass
            if not hit and '.' in q and float(q)>50: hit=any(re.fullmatch(r'\d\d\.\d+',x) and abs(round(float(x),len(q.split('.')[1]))-float(q))<1e-9 for x in nums)
            if hit:
                sc=len(w&(A.W(l)))+(3 if 'posthoc2/read_power' in s else 0)+(1 if s.startswith(R) else 0)
                if best is None or sc>best[0]: best=(sc,s,l)
    return best
