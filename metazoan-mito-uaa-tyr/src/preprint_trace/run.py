"""Step 3: write the number trace. usage (project folder):
  env/prefix/bin/python src/preprint_trace/extract.py MANUSCRIPT.tex
  env/prefix/bin/python src/preprint_trace/refs.py MANUSCRIPT.tex
  env/prefix/bin/python src/preprint_trace/run.py
Matches for text and Table 1 were read by hand; the sources chosen that way are in sources.py."""
import re,sys,csv,json,collections,os
import assign as A, sources as ov, table as TB, paths
doc=A.tok.doc; P=A.P; R='reports/REPORT.md'; PQ='results/posthoc2/'
refs=list(csv.DictReader(open(A.P+'results/posthoc2/preprint_refs_check.tsv'),delimiter='\t'))
rows=[]; untr=[]
def exists(s):
    ok=[]
    for x in s.split(';'):
        f=re.split(r'[:\s(]',x.strip())[0]
        if re.match(r'(results|reports|src|data|PREREG|env|README)',f): ok.append(os.path.exists(P+f.rstrip('/')))
        elif re.match(r'[a-z_0-9]+\.(tsv|txt|md)$',f): ok.append(os.path.exists(P+PQ+f) or os.path.exists(P+PQ+'secondcopy/'+f))
    return all(ok)
def add(i,v,loc,src,match,basis,ctx):
    if not exists(src): print('NOFILE',i,v,src)
    rows.append((v,loc,src,match,basis,ctx))
VER={'release 238':('reports/GATE_1_OVERLAP.md:8','source line'),'pyhmmer':('env/tools_outside_lockfiles.tsv','version list'),'HMMER v3.4':('env/conda-osx-arm64.lock','lockfile'),'HMMER v3.3.2':('env/conda-hmmer332-osx-64.lock','lockfile'),'HMMER 3.3.2':('data/pfam/Pfam-A_enone.hmm (header and DATE field)','file header'),
'cmsearch':('env/conda-osx-arm64.lock','lockfile'),'refseq89m':('PREREG.md:53','source line'),'ARAGORN':('env/conda-osx-arm64.lock','lockfile'),'tRNAscan':('env/conda-osx-arm64.lock','lockfile'),'OrthoDB':('data/manifests/taxonomy_orthodb.sha256 (odb13v0)','file name'),'miniprot':('env/conda-osx-arm64.lock; env/prefix (0.18-r281)','lockfile and version printed by the tool'),
'DIAMOND':('env/tools_outside_lockfiles.tsv','version list'),'FastTree':('env/tools_outside_lockfiles.tsv','version list'),'BLAST':('env/conda-osx-arm64.lock','lockfile'),'MAFFT':('env/conda-osx-arm64.lock','lockfile'),'minimap2':('env/tools_outside_lockfiles.tsv','version list'),'samtools':('env/tools_outside_lockfiles.tsv','version list'),'SciPy':('env/conda-osx-arm64.lock','lockfile'),'SRA Toolkit':('env/tools_outside_lockfiles.tsv','version list'),'SPAdes':('env/tools_outside_lockfiles.tsv','version list'),'GPT':('author statement','not traceable here; AI disclosure left to the author'),'(v3)':('manuscript version label','label'),'v3 archive':('manuscript version label','label')}
DATE={'2026-11-23':('PREREG.md commit 5aa5774 (git log); reports/GATE_1_OVERLAP.md','git log'),'2026-11-24':(R+':172','source line'),'2026-11-25 07:26':(R+':296; results/posthoc2/roty_criteria.md','source line'),'2026-11-25':(R+':71,172 ; results/posthoc2/README.md','source line'),'2021-11-23':('data/pfam/Pfam-A_enone.hmm (DATE field)','file header')}
HASH={'e23bc13f37ca48647b86cdc6b3d70abc5f7e0b69c71968de05e7b32285b4fa0c':('PREREG.md (sha256 recomputed)','recomputed'),'5aa5774':('git commit, "PREREG locked"','git log'),'863359e':('PREREG.md:20','source line')}
def RH(h,t,s):
    # the repository commits are named in the manuscript after this file is committed, so their values are not written here
    if h in HASH: return HASH[h]+(h,)
    return ('git: the project commit that holds this file, or its public copy','checked by hand after the push','(commit hash)')
accs=set()
for f in ('reports/REPORT.md','PREREG.md','data/manifests/download_ledger.tsv','results/mitogenome_master.tsv'):
    accs|=set(re.findall(r'[A-Z]{1,6}_?\d{5,}',open(P+f,encoding='utf8').read()))
more=''
import glob
for f in glob.glob(P+'reports/*.md')+glob.glob(P+'reports/dossiers/*.md')+[x for x in glob.glob(P+PQ+'*') if os.path.isfile(x)]+glob.glob(P+PQ+'secondcopy/*')+[P+'results/gate3_unverified/master.tsv']: more+=open(f,encoding='utf8',errors='replace').read()
def loc(i):
    a,b,t=doc[i]; 
    return f"{a[:70]} | {b}"
T3=[]
J=list(csv.DictReader(open(P+PQ+'sixfish_polya_junctions.tsv'),delimiter='\t'))
for r in csv.DictReader(open(P+PQ+'sixfish_terminal_uaa.tsv'),delimiter='\t'):
    ok={r[k] for k in ('UAA_at_conserved_Tyr','conserved_Tyr_sites','UAU','UAC','other_codons','encoded_UAA','incomplete','transcript_supported_UAA')}|{j['tail_reads_UAA'] for j in J if j['species']==r['species'] and j['verdict'].startswith('transcript-supp')}
    if 'glossodonta' in r['species']: ok.add('132')   # COX2 AGG tails, sixfish_polya_junctions.tsv codon_counts
    T3.append((r['species'],ok))
refidx=-1; refrec=None
for i,(a,b,t) in enumerate(doc):
    if a=='References':
        if t.startswith('\\item'): refidx+=1; refrec=refs[refidx]
        for v,k,s,p in A.tok.tokens(t):
            if refidx<0: continue
            src='results/posthoc2/preprint_refs_check.tsv ('+refrec['source']+', doi:'+refrec['doi']+')'
            add(i,v,f'References | item {refidx+1}',src,'yes','bibliographic record fetched 2026-11-25',t[max(0,s-25):s+30])
        continue
    if a=='Header' and re.match(r'\{|\\',t) and not re.search(r'C33|v3',t): continue
    if a=='Main Tables' and re.match(r'>\{|\\multicolumn',t): continue
    tks=A.tok.tokens(t)
    if re.match(r'^\{?\\(scriptsize|small|footnotesize|color)',t): tks=[]
    fk=[k for p,ks in A.FE if re.search(p,t) for k in ks]; fv={}
    for k in fk:
        for v in A.F[k][0]: fv.setdefault(A.norm(v),A.F[k][1])
    rest=[]; words=[]
    row=next((r for r in TB.rows[:17] if i in r),None)
    trow=next((r for r in TB.rows[17:] if i in r),None); tt=(doc[trow[0]][2]+' | '+t) if trow else t
    if trow: fk=[k for p,ks in A.FE if re.search(p,tt) for k in ks]; fv={}
    if trow:
        for k in fk:
            for v in A.F[k][0]: fv.setdefault(A.norm(v),A.F[k][1])
    ovs=[(vals,s,bs) for rx,vals,s,bs in ov.OV if re.search(rx,tt)]
    t3=next((x for x in T3 if trow and x[0] in doc[trow[0]][2]),None)
    # per-token
    done_hash=set(); done_ver=set(); done_date=set(); done_acc=set()
    L1=TB.src_lines(row) if row else None
    auto=None
    for v,k,s,pc in tks:
        ctx=t[max(0,s-35):s+35]; LO=loc(i)+(f' ({doc[row[0]][2].split(" (")[0]})' if row and i!=row[0] else '')
        if k=='hash':
            m=next((x for x in re.finditer(r'[0-9a-f]{7,64}',t) if x.start()<=s<x.end()),None); h=m[0]
            if (h,m.start()) in done_hash: continue
            done_hash.add((h,m.start())); add(i,RH(h,t,m.start())[2],LO,*[*RH(h,t,m.start())[:1],'yes',RH(h,t,m.start())[1]],t[max(0,m.start()-30):m.end()]); continue
        if k=='version' and t[max(0,s-7):s]=='minimap': add(i,v,LO,'','label','part of a name, not a quantity',ctx); continue
        if k=='version':
            key=next((x for x in sorted(VER,key=len,reverse=True) if (j:=t.rfind(x.split()[0] if x not in('(v3)','v3 archive','release 238','HMMER v3.4','HMMER v3.3.2','HMMER 3.3.2','SRA Toolkit') else x,max(0,s-24),s+14))>=0),None)
            m=re.match(r'v?[\d.]+(?:-r\d+)?',t[s:]) ; pos=s
            if any(abs(pos-q)<9 for q in done_ver): continue
            done_ver.add(pos); sr,bs=VER.get(key,('','')); add(i,t[s:s+re.match(r'[\d.]+(?:-r\d+)?',t[s:]).end()].rstrip('.'),LO,sr,'yes' if 'author' not in sr else 'untraced',bs,ctx); 
            if not key: untr.append((i,v,ctx,'version'))
            continue
        if k=='date':
            m=next(x for x in re.finditer(r'20\d\d-\d\d-\d\d(?: \d\d:\d\d)?',t) if x.start()<=s<x.end())
            if m.start() in done_date: continue
            done_date.add(m.start()); add(i,m[0],LO,DATE[m[0]][0],'yes',DATE[m[0]][1],ctx); continue
        if k=='accession':
            m=next((x for x in re.finditer(r'(?:[A-Z]{1,6}_?\d{5,}(?:\.\d+)?|PV962|Pfam-A_enone)',t) if x.start()<=s<=x.end()),None)
            if m is None or m.start() in done_acc: continue
            done_acc.add(m.start()); root=m[0].split('.')[0]; ok=root in accs or root in more
            src=R if root in open(P+R,encoding='utf8').read() else 'results/mitogenome_master.tsv' if root in accs else 'reports/dossiers/ or results/posthoc2/'
            add(i,m[0],LO,src,'yes' if ok else 'untraced','accession present',ctx); 
            if not ok: untr.append((i,m[0],ctx,'accession'))
            continue
        if k=='path' or k=='ref' or k=='url': 
            continue
        if k=='cite': add(i,v,LO,'References list (bibliographic record)','yes','citation year',ctx); continue
        if k=='label' or re.search(r'(?:minimap|SHA-|MITOS|BLAST|Pfam|S|T|C|A|B|fig)$',t[max(0,s-7):s]) and len(v)<4 or (re.match(ov.LISTNO,t) and s==1):
            add(i,v,LO,'','label','part of a name, not a quantity',ctx); continue
        if row and re.search(r'Bacterial or phage contigs at',t): add(i,v,LO,PQ+'secondcopy/highcopy_summary.txt:24-28 (contig depths 93.3, 90.2 and 29.3 against 35: 2.7, 2.6 and 0.8 times)','yes','recomputed',ctx); continue
        if row and 'HiFi 56/0' in t and v in('31','9.2'): add(i,v,LO,R+':227','yes','source line',ctx); continue
        # number
        if t3 and i!=trow[0]:
            if A.norm(v) in t3[1]: add(i,v,LO,PQ+'sixfish_terminal_uaa.tsv; sixfish_polya_junctions.tsv','yes','recomputed from the table',ctx)
            else: untr.append((i,v,ctx,'t3'))
            continue
        if row:
            col=int(doc[i][1].split()[-1]); 
            if col<=3: add(i,v,LO,'results/cold_review/tiers.tsv; results/mitogenome_master.tsv','yes','tier or code column',ctx); continue
            h=TB.find(v,pc,L1,t+' '+doc[row[0]][2])
            if col==6 and h: 
                hh=[(q,l) for q,l in L1 if re.match(R+r':(28[4-9]|29[0-2])$',q) and re.search(r'(?<![\d,.])'+v+r'(?![\d])',l)]; src=hh[0][0] if hh else PQ+'secondcopy/highcopy_summary.txt'
                add(i,v,LO,src,'yes','source line',ctx); continue
            if v=='68': add(i,v,LO,PQ+'secondcopy/highcopy_summary.txt:3; '+R+':284','yes','source line',ctx); continue
            if h:
                src=h[1]
                if re.search(R+r':113$|:172$|:14$',src): src={'conspec':PQ+'conspecific_specimen_audit.tsv'}.get('conspec') if re.search('conspec|other records|comparison',t) else (PQ+('read_power_check4.tsv' if col==5 else 'read_power.tsv') if col in(4,5) else src)
                if col==5: src=PQ+'read_power_check4.tsv'
                elif col==4 and (re.search(r'/0|/1|%\)',t) and not re.search('COX1',t[max(0,s-3):s+20])) and not re.search('conspec|congener|Ophichth',t[s:s+22]): src=PQ+'read_power.tsv'
                add(i,v,LO,src,'yes','source cell' if 'tsv' in src else 'source line',ctx)
            else: untr.append((i,v,ctx,'t1'))
            continue
        hit=False
        for vals,sr,bs in ovs:
            if vals!='*' and v in vals.split(','): add(i,v,LO,sr,'yes',bs,ctx); hit=True; break
        if hit: continue
        if A.norm(v) in fv: add(i,v,LO,fv[A.norm(v)],'yes','recomputed from the table',ctx); continue
        for vals,sr,bs in ovs:
            if vals=='*': add(i,v,LO,sr,'yes',bs,ctx); hit=True; break
        if hit: continue
        rest.append((v,pc,ctx,LO))
    if rest:
        c2=t+' '+(doc[i-1][2] if i and doc[i-1][0]==a else '')+' '+a
        r=A.auto(t,c2,[(v,p) for v,p,_,_ in rest])
        for v,p,ctx,LO in rest:
            if v in r: add(i,v,LO,r[v][0],'yes','source line',ctx)
            else: untr.append((i,v,ctx,'num'))
# ---- spelled-out counts
WOV=[(r'two bonefishes|two bonefish records|two species or',PQ+'prior_art_c33_fish.md'),(r'^Two recent studies|two bonefish|neither change|two-part|^These two sites|first two follow',''),(r'two cold reviews','reports/COLD_REVIEW.md; reports/COLD_REVIEW_2.md'),
(r'two-sided binomial|two-sided binomial limits','LABEL'),(r'^Metazoan-wide confirmatory sweep: zero|zero resolved non-panel|zero new reassigned|zero-new-clade|zero genomes met','reports/REPORT.md:99-100,142,168'),(r'the other two are R\. ren','results/mitogenome_master.tsv'),(r'two callable genes','reports/REPORT.md:108'),
(r'^Four checks|^Four scope|^Five constraints|two mandatory caveats|first two follow','SELF'),(r'two Albula records','reports/REPORT.md:71; results/gate3/notes_albula_records.md'),(r'five of the six fish|all five species|^Five retain|five Tier 1 fish',PQ+'conspecific_specimen_audit.tsv; reports/REPORT.md:229'),(r'four additional records',PQ+'conspecific_specimen_audit.tsv'),
(r'does not require zero','LABEL'),(r'two populations',PQ+'read_power.tsv'),(r'the two genes share','reports/REPORT.md:280'),(r'zero mitochondrial genes',PQ+'secondcopy/beanii_contig.txt; reports/REPORT.md:286'),(r'altered zero counts','reports/REPORT.md:75,276'),
(r'two unannotated assemblies|RNA-seq exists for two|other four fish|of four of the six fish|four "annotation only"',PQ+'sixfish_terminal_uaa.tsv'),(r'within five residues|six residues longer|fewer than five',PQ+'sixfish_polyA_criteria.md; sixfish_terminal_genes.tsv; sixfish_polya_junctions.tsv'),
(r'found two copies reading T','reports/REPORT.md:276'),(r'^Three internal TAG','reports/REPORT.md:272; results/gate3/roty_tag_sites.txt'),(r'five criteria',PQ+'roty_criteria.md'),(r'^In pairwise global|^A single MAFFT|^The three are',PQ+'tylench_trny_identity.tsv; reports/REPORT.md:266'),
(r'two UAA=Tyr calibration|S216 occurred in two','reports/REPORT.md:16,133'),(r'excluding two MTRF1|^Two additional T216','reports/REPORT.md:164; '+PQ+'mtrf1a_stop_denominator.tsv'),(r'three canonical tylenchid|^Two anomalous|two frameshifted','reports/REPORT.md:161; '+PQ+'pvulnus_models.txt'),
(r'zero starfish','reports/REPORT.md (T4, Li et al. 2019 figures)'),(r'changed zero screen','reports/REPORT.md:278'),(r'three further lineages','results/cold_review/tiers.tsv (4 lineages less Albula)'),(r'^Three additional eels','results/cold_review/tiers.tsv'),(r'two C33 fish with nuclear data','reports/REPORT.md:154-164'),(r'pre-repeat-resolution contigs in two fish','reports/REPORT.md:115'),
(r'four interim','results/gate3/interim/LOG.md; reports/REPORT.md'),(r'three corrupted','PREREG.md:16; results/mitogenome_master.tsv (OZ576062); results/gate3_unverified/master.tsv (PP708078, PP708091)'),(r'^Six-frame','LABEL'),(r'five early read|six assemblies|two pilots|two tool installs','data/manifests/download_ledger.tsv (rows without a SHA-256)'),(r'^Seven tools','env/tools_outside_lockfiles.tsv'),
(r'six fish|the six|Six fish|Six C33|six best|six transcripts|other six|six or seven|six Tier|the six','results/cold_review/tiers.tsv; '+PQ+'sixfish_terminal_uaa.tsv'),(r'three unconfirmed','results/cold_review/tiers.tsv'),(r'only two screen positives','results/mitogenome_master.tsv (code 5)'),(r'both series|two series',PQ+'fig3_bootstrap.tsv'),(r'four lineages|in four fish lineages','results/cold_review/tiers.tsv'),(r'two species',PQ+'sixfish_terminal_uaa.tsv'),(r'seven eel','results/posthoc2/elopomorpha_c33.tsv'),(r'poly\(A\) (tails|ends)|at most two|^Each ends',PQ+'sixfish_terminal_uaa.tsv'),(r'approximately five origins','reports/REPORT.md:35'),(r'five replicates',PQ+'fig3_bootstrap.tsv'),(r'two relatives|lineage support','results/cold_review/tiers.tsv'),(r'Zero of 47',PQ+'mtrf1a_stop_denominator.tsv')]
PATW=r'\b(?:'+'|'.join(x for x in A.WN if x!='one')+r')\b'
for i,(a,b,t) in enumerate(doc):
    if a in('References','Acknowledgements and AI-Assistance Disclosure') or (a=='Header' and 'C33' not in t): continue
    for m in re.finditer(PATW,t,re.I):
        ctx=t[max(0,m.start()-35):m.end()+35]; near=t[max(0,m.start()-30):m.end()+45]
        hit=next(((rx,s) for rx,s in WOV if any(x.start()<=m.end() and m.start()<=x.end()+2 for x in re.finditer(rx,t))),None) or next(((rx,s) for rx,s in WOV if re.search(rx,t) and s),None)
        fk=[k for p,ks in A.FE if re.search(p,t) for k in ks]; fv={}
        for k in fk:
            for v in A.F[k][0]: fv.setdefault(A.norm(v),A.F[k][1])
        d=str(A.WN[m[0].lower()]); LO=loc(i)
        if hit and hit[1]=='SELF': add(i,m[0],LO,'','yes','counts items listed in the manuscript itself',ctx)
        elif hit and hit[1] in('LABEL',''): add(i,m[0],LO,'','label','wording, not a measured count',ctx)
        elif hit: add(i,m[0],LO,hit[1],'yes','source line',ctx)
        elif d in fv: add(i,m[0],LO,fv[d],'yes','recomputed from the table',ctx)
        else: untr.append((i,m[0],ctx,'word'))
# ---- figure annotations
import html
def svg(f):
    s=open(P+f'reports/figures/{f}.svg',encoding='utf8').read()
    return [html.unescape(re.sub(r'\s+',' ',re.sub(r'<[^>]+>','',x)).strip()) for x in re.findall(r'<text[^>]*>(.*?)</text>',s,re.S)]
def fv(*ks):
    o={}
    for k in ks:
        for v in A.F[k][0]: o.setdefault(A.norm(v),A.F[k][1])
    return o
FB={'fig1_tree':(fv('tiers','panelpass','bases','RS'),[(r'^UAA=Tyr, panel','results/posthoc2/locus_gate_corroboration.tsv; results/gate2/locus_audit.tsv'),(r'11 gains','reports/REPORT.md:99'),(r'6 origins|4 origins|1 origin|about 5 origins|genomes?,|Tier|genus','results/cold_review/tiers.tsv; reports/REPORT.md:262 and T3 counter (src/origins.py)'),(r'REPORT 7\.2\.1','reports/REPORT.md'),(r'3/5','results/posthoc2/locus_gate_corroboration.tsv')]),
'fig2_mtrf1a_loop':(fv('mtrf'),[(r'^(203|216|196|222)$|human 196|human 216|Human numbering|Q at 203|at 216','PREREG.md:10,79'),(r'27 aa|26 aa|\+1 aa|\+0|\+2 aa|\+14|47-nt|PREREG 2\.4|27 aa','reports/REPORT.md:161; results/gate4/loop_report.txt'),(r'NP_061914','PREREG.md:10'),(r'C33|G33|U33','LABEL')]),
'fig3_sensitivity':(fv('b_counter','b_counter+C','cp65','panelpass'),[(r'^\d/\d$|0\.8 x min|14 / 20|20 / 20|20 known','results/gate2/lolo.tsv; results/posthoc2/locus_gate_corroboration.tsv'),(r'10,000|resamples|5\.5%|5 per genome','results/posthoc2/fig3_bootstrap.tsv; src/figures.py'),(r'29 held-out|0 of 340|1\.1%','reports/REPORT.md:48'),(r'PREREG','LABEL'),(r'^0\.\d$|^1\.0$|^\d+%$','LABEL')])}
for f,(d,rx) in FB.items():
    for x in svg(f):
        for v,k,s,p in A.tok.tokens(x):
            LO=f'Figure {"123"[list(FB).index(f)]} | drawn label ({f}.svg)'
            if k in('label','accession','ref','version') or re.search(r'(?:Tier|PREREG|REPORT|Gate|GATE_|§) ?[\d.]*$',x[max(0,s-9):s]): add(-1,v,LO,'','label','part of a name, not a quantity',x[:90]); continue
            if re.fullmatch(r'0\.\d|1\.0|\d+%',x) and f=='fig3_sensitivity' and x in('0.0','0.2','0.4','0.6','0.8','1.0'): add(-1,v,LO,'','label','axis tick',x); continue
            h=next((s2 for r2,s2 in rx if re.search(r2,x)),None)
            if A.norm(v) in d and not (h and h!='LABEL'): add(-1,v,LO,d[A.norm(v)],'yes','recomputed from the table',x[:90])
            elif h and h!='LABEL': add(-1,v,LO,h,'yes','source line',x[:90])
            elif A.norm(v) in d: add(-1,v,LO,d[A.norm(v)],'yes','recomputed from the table',x[:90])
            elif h: add(-1,v,LO,'','label','part of a name, not a quantity',x[:90])
            else: untr.append((f,v,x[:90],'fig'))
rows[:]=[r[:3]+(('no',) if r[4]=='MISMATCH' else (r[3],))+(('differs from source',) if r[4]=='MISMATCH' else (r[4],))+r[5:] for r in rows]
NOTE=[(lambda v,l,c:v=='402','wording corrected: was "402 records are identical in sequence to another record"; 402 counts the extra copies, 768 records in 366 sets share a sequence (REPORT line 25 has the same wording)'),
(lambda v,l,c:'12 weak (11 rej' in c,'corrected: cell read "1 weak"; 12 met weak, 11 rejected, 1 supported'),
(lambda v,l,c:'ERR14792883' in c and ('has only' in c or 'bounds of 78' in c),'added: the range 4.44-7.22% left out this run, although the sentence says every read check'),
(lambda v,l,c:v in('16,902','814') and 'CI' not in c,'matches the locked output; NC_085561 has a merged taxid (2479863, now 3905135) that is code 5 but is counted as unassigned: 16,903 and 813 if merged taxids are followed'),
(lambda v,l,c:v in('56','5.21','31','9.2') and 'HiFi' in c and 'glossodonta' in l,'caveat added: 31 of the 56 reads are mitochondrial (0 of 31, bound 9.2%)'),
(lambda v,l,c:v=='90.8' and 'sector' in c,'REPORT: 90.8% over 97% of the genome; conspecific_specimen_audit.tsv has 91.0% over 16,106 bp (another alignment)'),
(lambda v,l,c:v in('681','666','15','182.09') and 'Data and code' in l,'corrected: the sentence said the checksums left out two 0.9-GB pilots; the pilots are 0.46 and 0.45 GB, inside the total, and 15 rows have no SHA-256. Three ledger rows were added for the reference look-ups (678 before)'),
(lambda v,l,c:v=='two' and 'Albula records' in c,'corrected: read "proved to be mislabelled"; the sources say "probably mislabelled" (OZ556087) and "doubtful" (AP002973)'),
(lambda v,l,c:'Bacterial or phage contigs' in c and v in('0.8','2.7'),'corrected: read "High bacterial DNA load (~3x mt depth)"; the three identified contigs are at 2.7, 2.6 and 0.8 times the mitochondrial depth'),]
out=paths.ROOT+'results/posthoc2/preprint_number_trace.tsv' if len(sys.argv)<2 else sys.argv[1]
with open(out,'w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t',lineterminator='\n'); w.writerow(['number','location','source_file','match','basis','context','note'])
    for v,l,s,m,b,c in rows:
        n=next((t for f,t in NOTE if f(v,l,c)),'')
        w.writerow([v,l,s,m,b,'' if v=='(commit hash)' else re.sub(r'\b(?!5aa5774|863359e|e23bc13)(?=[0-9a-f]*\d)(?=[0-9a-f]*[a-f])[0-9a-f]{6,40}\b','(hash)',' '.join(c.split())) if 'Data and code' in l or l.startswith('Supplementary') or 'Header' in l else ' '.join(c.split()),n])
print(len(rows),collections.Counter(r[3] for r in rows),len(untr),untr)
