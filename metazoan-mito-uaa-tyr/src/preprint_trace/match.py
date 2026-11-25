"""Index the numbers in REPORT, PREREG, gate reports, dossiers and result files, for matching a sentence to a source line."""
import re,json,glob,os,sys,collections
import tok,paths
P=paths.ROOT
F=['reports/REPORT.md','PREREG.md']+sorted(glob.glob(P+'results/posthoc2/*'))+sorted(glob.glob(P+'results/posthoc2/secondcopy/*'))+sorted(glob.glob(P+'results/cold_review/*'))+sorted(glob.glob(P+'reports/dossiers/*.md'))+sorted(glob.glob(P+'reports/GATE_*.md'))+[P+'results/gate3/'+x for x in ('tests_fixed.txt','post.txt','posthoc_stem.txt','tylench_cox1_summary.txt','conspecific_b.txt','craniella_audit.txt','eel_origins.txt','notes_albula_records.md','numt_lycenchelys.txt','numt_craniella_infr.txt','numt_craniella_zetl.txt','roty_tag_sites.txt','roty_SRR26324223.txt','roty_ERR744555.txt')]+[P+'results/gate3/interim/LOG.md']+[P+'results/gate3_unverified/tests.txt',P+'results/gate4/loop_report.txt',P+'results/gate4/s216_reads.txt',P+'results/gate4/reads_SRR3350473.txt',P+'results/gate2/lolo.tsv',P+'results/gate2/trny_phylum_failure.tsv']
W=lambda t:set(x.lower() for x in re.findall(r'[A-Za-z][A-Za-z\-]{3,}',t))
STOPW=W('with that this from have were which than each when also only other into same both over read reads record records genome genomes these there their they does against across under more less least most between after before where while tRNA-Tyr mitochondrial UAA Stop')
L=[]
for f in F:
    f=f if f.startswith('/') else P+f
    if not os.path.isfile(f) or f.endswith('title_audit.tsv') or os.path.getsize(f)>900000 or f.endswith('.timestamp'): continue
    r=f.replace(P,''); lines=open(f,encoding='utf8',errors='replace').read().split('\n'); hdr=lines[0] if f.endswith('.tsv') else ''
    for i,l in enumerate(lines):
        if not l.strip(): continue
        parts=[l] if (f.endswith('.tsv') or not r.startswith('reports/')) else re.split(r'(?<=[.;:])\s+(?=[A-Z(\[*"≥0-9`])',l)
        # windows of 3 sentences for prose
        if len(parts)>1: parts=[' '.join(parts[max(0,k-1):k+2]) for k in range(len(parts))]
        for p in parts:
            nums=set(x.rstrip(',').replace(',','') for x in re.findall(r'(?<![A-Za-z_])\d[\d,]*(?:\.\d+)?',p))
            if nums: L.append((r,i+1,p,nums,W(p+' '+hdr)-STOPW))
idx=collections.defaultdict(list)
for k,x in enumerate(L):
    for n in x[3]: idx[n].append(k)
def fl(x):
    try: return float(x)
    except Exception: return None
def variants(v,pct):
    v=v.replace(',',''); out={v}
    if '.' in v: out.add(v.rstrip('0').rstrip('.'))
    if pct and fl(v) is not None:
        d=len(v.split('.')[1]) if '.' in v else 0
        out.add(f"{fl(v)/100:.{d+2}f}"); 
    return out
def best(t,cands,ctx,n=3):
    ns=[(v,p) for v,k,s,p in cands]; score=collections.Counter(); cov=collections.defaultdict(set)
    for v,p in ns:
        for q in variants(v,p):
            for k in idx.get(q,()): cov[k].add(v)
    w=W(t)-STOPW
    for k,c in cov.items():
        rare=sum(1 if len(x.replace(',',''))<3 and '.' not in x else 3 for x in c); sh=len(w&L[k][4])
        score[k]=rare*2+min(sh,8)+(2 if L[k][0].startswith('results/') else 0)+(1 if L[k][0]=='reports/REPORT.md' else 0)-(0 if sh else 8)
    return [(k,cov[k]) for k,_ in score.most_common(40)],w
def assign(t,cands):
    b,w=best(t,cands,None); left=[v for v,k,s,p in cands]; res={}
    for k,c in b:
        new=[v for v in c if v not in res]
        if new and (len(w&L[k][4])>=1):
            for v in new: res[v]=k
    return res
