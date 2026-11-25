"""Assign a sentence's numbers to recomputed facts (FE) or, failing that, to the best matching source line (auto)."""
import sys,re,json,collections,csv,os
import paths
D=paths.WORK
import tok,facts,match
F=facts.F; L=match.L; W=match.W; STOPW=match.STOPW; P=match.P
WN=dict(zip('zero one two three four five six seven eight nine ten eleven twelve thirteen'.split(),range(14)))
# (regex on entry, fact keys that may supply values, for that entry)
FE=[(r'17 non-U33 screen positives|screen positives: zero|Of 28,583|Among the 28,583|among 28,583|17/28,583',['RS','bases','T1','dist','RStaxa']),
(r'9 of 12 vertebrate',['elo']),(r'Six fish from four',['six','tiers']),(r'screen of 32,982|32,982 (metazoan|eligible|=)|32,623 \+ 358',['inc','conf','panelrec','panel','panelsp','conftaxa','panel360','panel26']),
(r'98 (records|further)',['excl']),(r'The set comprises',['conf','conftaxa','panelrec','panel','panelsp']),(r'other 26 of the 360|control panel of 360|Calibration on the 360',['panel360','panel26','inc','panel']),
(r'539 taxa',['multi']),(r'402 records',['dup','d303']),(r'16,902 (code-5|resolved)',['codes','Bivalvia','Mollusca']),(r'0 of 47|[Zz]ero of 47|all 47|47 (species|UAA|are|have)',['mtrf']),
(r'14 of 20|14/20',['panelpass']),(r'Anomia \(5/5\)',['panelpass']),(r'planting|Planting|counter alone|63 of 65|bootstrap|locked rule, which|At 5%|planted-null',['b_counter','b_counter+C','cp65']),
(r'Cnidaria \(0/892',['cov_Cnidaria','cov_Ctenophora','cov_Chaetognatha']),(r'Acanthocephala \(20',['cov_Acanthocephala','cov_Bryozoa','cov_Onychophora','cov_Gnathostomulida','cov_Nematoda']),(r'Rotifera \(8|unresolved rotifer',['cov_Rotifera']),
(r'Across the 32,623 confirmatory|UAA calls \(n|32,510 Stop',['conf','ucalls']),(r'Codetta returned|Codetta "\?"',['codetta','RS']),(r'27,618 NCBI',['RStaxa','bases','dist','T1']),(r'^The other 17|^14 C33 and 3',['bases']),(r'The locked screen-positive rate \(T1\)',['T1']),
(r'2,018 Stop genomes',['stopn33']),(r'review class holds|of the 2,020',['rev','stopn33']),(r'By NCBI mitochondrial translation|0 of 1,491 molluscs|Code 2: 13|Mollusca: 0',['codes','bases','Mollusca','Bivalvia','Echinodermata','Asteroidea']),(r'By phylum, they',['Chordata','Nematoda','Porifera']),(r'By provenance',['prov_RefSeq','prov_INSDC_only','prov_other_assembly','prov_DToL']),(r'Tightening the gene-overlap',['pass5','pass0']),
(r'^The locked screen-positive rate \(17|The 17 locked|17 locked confirm',['RS','bases']),(r'13 of the 17|14 hits|3 hits|6 Tier 1 \(|3 Tier 2|2 Tier 3|2 unconfirmed|1 dropped|^1 Tier 1|C33 \(14\)|Total: 13|five of the six|six Tier 1|Tier 1 fish|Six C33 fish',['tiers']),
(r'^Neither code-5|no confirmed genome',['Mollusca','Bivalvia']),(r'only 5 of the 23|all 23 primary',['power']),(r'0\.25–2\.81|minor U33 share',['u33best']),(r'left only 5 U33',['u33tier1']),
(r'None has an internal UAA|Each ends four|run past the stop|11 or 12 of 13|poly\(A\) ends confirming|UAG ends at most',['six','six2']),(r'78 ends|14 U33 UAA=Stop',['ends']),(r'tails place|the same holds for|tail reads each|17–110|132 of 134',['tails','six']),
(r'32 internal TAA|7 TAT|34 conserved',['roty']),(r'C33 calls cluster|9 of 239|7 of 219|Elopiformes|seven eel',['elo']),(r'77 of 77|323 echin',['Asteroidea','Echinodermata']),(r'182\.09',['ledger']),(r'all 32,623|across 32,623',['conf']),(r'28,566',['bases']),(r'17 screen|17 hits|the 17 |17/17',['bases']),(r'28,583',['RS']),
(r'31,564',['conftaxa'])]
INV={r['accession']:r['derived_from'] for r in csv.DictReader(open('results/refseq_mito_inventory.tsv'),delimiter='\t')}; BASE={r['acc'].split('.')[0] for r in facts.C}
facts.put('d303',[sum(INV.get(r['acc'],'').split('.')[0] in BASE for r in facts.C)],'results/refseq_mito_inventory.tsv (derived_from) + results/mitogenome_master.tsv'); facts.put('six2',['2','1'],'results/posthoc2/sixfish_terminal_uaa.tsv')
EXC=r'dossiers/README|title_audit|GATE_0'
def norm(v): return v.replace(',','')
def auto(text,ctx,vals,pref=None,bar=1):
    """assign each value to the best source line: shared words with ctx, co-occurring values, REPORT preferred"""
    w=W(ctx)-STOPW; res={}; score={}; cov=collections.defaultdict(set)
    for v,p in vals:
        for q in match.variants(v,p):
            for k in match.idx.get(q,()): cov[k].add(v)
    for k,c in cov.items():
        if re.search(EXC,L[k][0]): continue
        sh=len(w&L[k][4]); rare=sum(1 if len(norm(x))<3 and '.' not in x else 3 for x in c)
        s=rare*2+min(sh,10)*1.5+(3 if L[k][0]=='reports/REPORT.md' else 0)+(2 if pref and re.search(pref,L[k][0]+' '+L[k][2][:400]) else 0)-(0 if sh else 20); score[k]=(s,sh)
    for k,(s,sh) in sorted(score.items(),key=lambda x:-x[1][0])[:60]:
        if sh<bar: continue
        for v in cov[k]:
            if v not in res: res[v]=(L[k][0]+f':{L[k][1]}',L[k][2],sh)
    return res
