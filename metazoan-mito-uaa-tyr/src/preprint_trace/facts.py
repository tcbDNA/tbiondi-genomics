"""Recompute the manuscript's main counts from the raw tables (needs the local sweep FASTA and the NCBI taxonomy dump)."""
import csv,collections,math,json,re,os,hashlib
import paths
P=paths.ROOT; os.chdir(P)
T=lambda p:list(csv.DictReader(open(p),delimiter='\t'))
F={}   # name -> (set of allowed strings, source)
def put(k,vals,src): F[k]=([str(v) for v in (vals if isinstance(vals,(list,tuple,set)) else [vals])],src)
R=T('results/mitogenome_master.tsv'); MS='results/mitogenome_master.tsv'
C=[r for r in R if r['status']=='included' and r['role']=='confirmatory' and r['acc']!='AP002973.1']; S=[r for r in C if r['uaa']=='Stop']; RS=[r for r in S if r['n33']=='resolved']
inc=[r for r in R if r['status']=='included']; rc=collections.Counter(r['role'] for r in inc)
c=lambda n:f"{n:,}"
put('inc',c(len(inc)),MS); put('conf',c(len(C)),MS); put('conftaxa',c(len({r['taxid'] for r in C})),MS); put('panelrec',rc['panel']+rc['panel_species'],MS); put('panel',rc['panel'],MS); put('panelsp',rc['panel_species'],MS); put('excl',sum(r['status']!='included' for r in R),MS)
pan=T('results/gate2_panel.tsv'); accs={r['acc'] for r in inc}; put('panel360',len(pan),'results/gate2_panel.tsv'); put('panel26',sum(r['panel_id'] not in accs for r in pan),'results/gate2_panel.tsv + '+MS)
u=collections.Counter(r['uaa'] for r in C); put('ucalls',[c(u['Stop']),u['ambiguous'],u['insufficient'],u['weak'],u['Tyr']],MS)
put('RS',c(len(RS)),MS); put('RStaxa',c(len({r['taxid'] for r in RS})),MS); b=collections.Counter(r['base'] for r in RS); put('bases',[c(b['U']),b['C'],b['G'],b['C']+b['G']],MS); n3=collections.Counter(r['n33'] for r in S); put('stopn33',[c(n3['stem_audit_review']),c(n3['unresolved'])],MS)
put('rev',c(sum(r['n33']=='stem_audit_review' for r in C)),MS); cod=collections.Counter(r['codetta'] for r in RS); put('codetta',[c(cod['?']),cod['E'],cod['Q']],MS)
def wil(k,n,z=1.959964): p=k/n; d=1+z*z/n; q=(p+z*z/(2*n))/d; e=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return q-e,q+e
w=wil(17,len(RS)); put('T1',[f"{100*17/len(RS):.3f}",f"{100*w[0]:.3f}",f"{100*w[1]:.3f}",f"{100*17/len({r['taxid'] for r in RS}):.3f}"],MS)
h=collections.defaultdict(list)
for r in C: h[hashlib.md5(''.join(open('data/sweep/fa/'+str(sum(map(ord,r['acc']))%200).zfill(3)+f"/{r['acc']}.fna").read().split('\n',1)[1].split()).upper().encode()).hexdigest()].append(r)
rsa={r['acc'] for r in RS}; ds=sum(any(x['acc'] in rsa for x in v) for v in h.values()); put('dist',[c(ds),f"{100*17/ds:.3f}"],MS+' + sweep FASTA (md5 of sequence)')
put('dup',[sum(len(v)-1 for v in h.values()),sum(len(v)==2 and {x['provenance'] for x in v}=={'RefSeq','INSDC_only'} for v in h.values())],MS+' + sweep FASTA (md5 of sequence)'); DUPALL=sum(len(v) for v in h.values() if len(v)>1)
tx=collections.Counter(r['taxid'] for r in C); put('multi',[sum(k>1 for k in tx.values()),c(sum(k-1 for k in tx.values() if k>1))],MS)
code={}
for l in open('data/ncbi/taxonomy/nodes.dmp'):
    f=[x.strip() for x in l.rstrip('\n').rstrip('|').split('\t|\t')]; code[f[0]]=f[8]
cc=collections.Counter(code.get(r['taxid']) for r in RS); nn=collections.Counter(code.get(r['taxid']) for r in RS if r['base']!='U'); put('codes',[c(cc['2']),nn['2'],cc['4'],nn['4'],c(cc['5']),nn['5'],sum(v for k,v in cc.items() if k not in('2','4','5')),12,1],MS+' + NCBI nodes.dmp')
for k,key,val in (('Chordata','phylum','Chordata'),('Nematoda','phylum','Nematoda'),('Porifera','phylum','Porifera'),('Mollusca','phylum','Mollusca'),('Bivalvia','class','Bivalvia'),('Asteroidea','class','Asteroidea'),('Echinodermata','phylum','Echinodermata')):
    x=[r for r in RS if r[key]==val]; put(k,[c(len(x)),sum(r['base']!='U' for r in x)],MS)
for k in ('RefSeq','INSDC_only','other_assembly','DToL'): x=[r for r in RS if r['provenance']==k]; put('prov_'+k,[c(len(x)),sum(r['base']!='U' for r in x)],MS)
for k in ('pass5','pass0'): x=[r for r in RS if r[k]=='True']; put(k,[c(len(x)),sum(r['base']!='U' for r in x)],MS)
for p in ('Cnidaria','Ctenophora','Chaetognatha','Acanthocephala','Bryozoa','Onychophora','Gnathostomulida','Nematoda','Rotifera'): x=[r for r in C if r['phylum']==p]; put('cov_'+p,[sum(r['n33']=='resolved' for r in x),len(x),len(x)-sum(r['n33']=='resolved' for r in x)],MS)
B=T('results/posthoc2/fig3_bootstrap.tsv'); BS='results/posthoc2/fig3_bootstrap.tsv'
pc=lambda v:(f"{100*float(v):.1f}".rstrip('0').rstrip('.')) if f"{100*float(v):.1f}".endswith('.0') else f"{100*float(v):.1f}"
for nm in ('counter','counter+Codetta(rep 0)'):
    x=[r for r in B if r['series']==nm]; put('b_'+nm[:8]+('C' if '+' in nm else ''),[pc(int(r['called'])/int(r['observations'])) for r in x]+[r['called'] for r in x]+[r['observations'] for r in x]+[r['genomes'] for r in x]+[pc(r['boot_lo']) for r in x]+[pc(r['boot_hi']) for r in x]+['0','5','10','25','50','100','95','65','325','five','one'],BS)
put('cp65',[f"{100*(1-0.025**(1/65)):.1f}"],'computed: two-sided exact limit on 0 of 65')
E6=T('results/posthoc2/sixfish_terminal_uaa.tsv'); S6='results/posthoc2/sixfish_terminal_uaa.tsv'
put('six',[min(int(r['conserved_Tyr_sites']) for r in E6),max(int(r['conserved_Tyr_sites']) for r in E6),c(min(int(r['aligned_codons']) for r in E6)),c(max(int(r['aligned_codons']) for r in E6)),min(int(r['UAU'])+int(r['UAC']) for r in E6),max(int(r['UAU'])+int(r['UAC']) for r in E6),13]+sorted({int(r['encoded_UAA'])+int(r['incomplete']) for r in E6})+sorted({int(r['encoded_UAA']) for r in E6})+sorted({int(r['incomplete']) for r in E6})+[int(r['transcript_supported_UAA']) for r in E6[:2]],S6)
G6=T('results/posthoc2/sixfish_terminal_genes.tsv'); ends=collections.defaultdict(list)
for r in G6:
    if r['stop_class']: ends[(r['species'],r['gene'])].append(r['profile_check']=='agrees' and r['cterm']=='conserved')
put('ends',[len(ends),sum(any(v) for v in ends.values()),len(T('results/posthoc2/sixfish_terminal_refs.tsv'))],'results/posthoc2/sixfish_terminal_genes.tsv; sixfish_terminal_refs.tsv')
J=T('results/posthoc2/sixfish_polya_junctions.tsv'); put('tails',[int(r['tail_reads_UAA']) for r in J if r['verdict'].startswith('transcript-supp')]+[x for r in J if r['gene']=='COX2' and r['verdict']=='transcript ends with AGG' for x in (r['codon_counts'].split(':')[1].split(';')[0],r['tail_reads_near_end'])],'results/posthoc2/sixfish_polya_junctions.tsv')
RP=T('results/posthoc2/read_power.tsv'); C4=T('results/posthoc2/read_power_check4.tsv'); put('power',[sum(r['supports_1pct_claim']=='yes' for r in RP),len(RP),len(C4)],'results/posthoc2/read_power.tsv; read_power_check4.tsv')
t1={'NC_037360.1','CM165726.1','NC_013627.1','NC_013436.1','NC_013628.2','CM163399.1'}; best={}
for r in C4:
    if r['acc'] in t1: best[r['acc']]=min(best.get(r['acc'],9),float(r['upper95_U33_share']))
put('u33best',[f"{100*min(best.values()):.2f}",f"{100*max(best.values()):.2f}"],'results/posthoc2/read_power_check4.tsv'); put('u33tier1',[sum(int(r['U33_reads']) for r in C4 if r['acc'] in t1|{'OZ318138.1'})],'results/posthoc2/read_power_check4.tsv')
EL={r['group']:r for r in T('results/posthoc2/elopomorpha_c33.tsv')}; put('elo',[EL['Elopomorpha']['C33'],EL['Elopomorpha']['resolved_stop_genomes'],EL['Elopomorpha']['gains_min'],EL['Vertebrata']['C33'],c(int(EL['other vertebrates']['resolved_stop_genomes'])),EL['other vertebrates']['C33'],EL['other vertebrates']['gains_min'],EL['Elopomorpha: Albuliformes']['C33'],EL['Elopomorpha: Albuliformes']['gains_min'],EL['Elopomorpha: Anguilliformes']['C33'],EL['Elopomorpha: Anguilliformes']['resolved_stop_genomes'],EL['Elopomorpha: Anguilliformes']['gains_min'],EL['Elopomorpha: Elopiformes']['resolved_stop_genomes'],EL['Elopomorpha: Notacanthiformes']['resolved_stop_genomes'],EL['genus Anguilla']['resolved_stop_genomes'],int(EL['Elopomorpha subfamily: Ophichthinae']['resolved_stop_genomes'])-int(EL['Elopomorpha subfamily: Ophichthinae']['C33']),EL['Elopomorpha subfamily: Ophichthinae']['resolved_stop_genomes'],0],'results/posthoc2/elopomorpha_c33.tsv')
TI=T('results/cold_review/tiers.tsv'); tk='tier'; cnt=collections.Counter((r['base'][0],r['tier'].split()[0]+(' '+r['tier'].split()[1] if r['tier'].startswith('Tier') else '')) for r in TI); put('tiers',[cnt[('C','Tier 1')],cnt[('C','Tier 2')],cnt[('C','Tier 3')],cnt[('C','unconfirmed')],cnt[('C','dropped')],cnt[('G','Tier 1')],cnt[('G','Tier 2')],cnt[('G','unconfirmed')],sum(v for k,v in cnt.items() if k[0]=='C'),sum(v for k,v in cnt.items() if k[0]=='G'),len(TI),sum(v for k,v in cnt.items() if k[1].startswith('Tier')),sum(v for k,v in cnt.items() if k[1]=='unconfirmed')],'results/cold_review/tiers.tsv'); TIERS=cnt
MD=T('results/posthoc2/mtrf1a_stop_denominator.tsv'); K=[r for r in MD if r['counted']=='yes']; put('mtrf',[len(K),sum(r['set']=='Gate 2 background' for r in K),sum(r['set']=='Gate 4' for r in K),sum(r['p216']=='S' for r in K)],'results/posthoc2/mtrf1a_stop_denominator.tsv'); MDN=len(MD)
PW=collections.defaultdict(lambda:[0,0]); lin={r['panel_id']:r['lineage'] for r in pan}
for r in T('results/posthoc2/locus_gate_corroboration.tsv'):
    if r['set']=='panel20': PW[lin[r['acc']]][1]+=1; PW[lin[r['acc']]][0]+=r['locked_status'].startswith('pass')
put('panelpass',[x for v in PW.values() for x in v]+[sum(v[0] for v in PW.values()),sum(v[1] for v in PW.values())],'results/posthoc2/locus_gate_corroboration.tsv')
for r in T('results/posthoc2/conspecific_specimen_audit.tsv')[:1]: AUDK=list(r)
L=[x for x in csv.reader(open('data/manifests/download_ledger.tsv',newline=''),delimiter='\t')][1:]; L=[x for x in L if x]; ib=lambda x:int(x) if x.strip().isdigit() else 0
put('ledger',[f"{sum(ib(r[3]) for r in L)/1e9:.2f}",len(L),sum(bool(re.match(r'[0-9a-f]{64}',r[4])) for r in L),sum(not re.match(r'[0-9a-f]{64}',r[4]) for r in L)],'data/manifests/download_ledger.tsv'); NOSUM=sum(not re.match(r'[0-9a-f]{64}',r[4]) for r in L)
D=json.loads(next(l for l in open('results/gate3/run/out.jsonl') if 'CM003310.1' in l[:40]))['uaa']; put('roty',[D['n_taa_internal'],D['n_taa_genes'],D['n_taa_tyrcol'],D['tyrcols'],D['tat'],D['tac'],D['tag_tyrcol'],D['llr_Y']],'results/gate3/run/out.jsonl (sha256 in data/manifests/gate3_raw.sha256)')
if __name__=='__main__':
    for k,v in F.items(): print(k,v[0][:16],'|',v[1][:50])
    print('DUPALL',DUPALL,'MDN',MDN,'NOSUM',NOSUM,dict(TIERS),{k:tuple(v) for k,v in PW.items()})
    print(AUDK)
