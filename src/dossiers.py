"""Build the outlier dossiers (reports/dossiers/*.md). Everything numeric is computed here from committed outputs; the "Checks" section is curated text whose source file is named. 17 non-U33 Stop candidates, Rotylenchulus, Pratylenchus goodeyi."""
import sys,csv,json,os,collections
sys.path.insert(0,"src"); os.environ["GATE3_UNLOCK"]="1"; import gate3_classify as K, gate3_post as P
D="reports/dossiers/"; os.makedirs(D,exist_ok=True)
B="TCAG"; T5="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"; AA={a+b+c:T5[i] for i,(a,b,c) in enumerate((a,b,c) for a in B for b in B for c in B)}
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
MAST=lambda f:{r["acc"]:r for r in csv.DictReader(open(f),delimiter="\t")}
G3=MAST("results/mitogenome_master.tsv"); G3U=MAST("results/gate3_unverified/master.tsv")
FE={}
for l in open("data/sweep/features_fixed.jsonl"):
    x=json.loads(l); FE[x["acc"]]=x["f"]
WANT=["NC_037360.1","NC_015546.1","CM165726.1","NC_005796.1","NC_013627.1","NC_013436.1","PZ495457.1","NC_013628.2","CM163399.1","OZ407450.1","OZ318138.1","PV357903.1","PX804232.1","KY021065.1","NC_082198.1","KU726971.1","NC_047225.1","CM003310.1"]
R={}
for l in open("results/gate3/run/out.jsonl"):
    d=json.loads(l)
    if d["acc"] in WANT: R[d["acc"]]=d
for l in open("results/gate3_unverified/run/out.jsonl"):
    d=json.loads(l)
    if d["acc"]=="ON525593.1": R["ON525593.1"]=d
LOOP=collections.defaultdict(list)
for r in csv.DictReader(open("results/gate4/loop.tsv"),delimiter="\t"):
    if r["covered"]=="True" and r["label"]=="MTRF1L": LOOP[r["source"]].append(r)
MT={"CM165726.1":["PROT_Albula_glossodonta_BYU","GEN_A_glossodonta_hap1"],"CM163399.1":["GEN_Lycenchelys"],"OZ407450.1":["GEN_Cr_zetlandica"],"OZ318138.1":["GEN_Cr_infrequens"],"CM003310.1":["GEN_R_reniformis","TSA_R_reniformis"]}
STAT={"NC_037360.1":"supported (C33)","NC_015546.1":"supported (C33; lineage evidence only)","CM165726.1":"supported (C33; lineage evidence)","NC_005796.1":"supported (C33; lineage evidence only, single 2004 record)","NC_013627.1":"supported (C33)","NC_013436.1":"supported (C33)","PZ495457.1":"supported (C33; own-specimen reads only)","NC_013628.2":"supported (C33)","CM163399.1":"supported (C33)","OZ407450.1":"supported (G33; nuclear-copy caveat)","OZ318138.1":"supported (G33; nuclear-copy caveat)","PV357903.1":"UNCONFIRMED (C33, one specimen)","PX804232.1":"UNCONFIRMED (C33, one specimen)","KY021065.1":"UNCONFIRMED (C33, one record, no reads)","NC_082198.1":"UNCONFIRMED (C33, one isolate, no reads)","KU726971.1":"DROPPED (C33; second U33 tRNA-Tyr-like locus on chromosome II)","NC_047225.1":"DROPPED (G33; no mitochondrial reads, no independent record)","CM003310.1":"SENSITIVITY (UAA=Tyr, blocked by the trnY stem audit)","ON525593.1":"SENSITIVITY (supported weak; UNVERIFIED record, not confirmed)"}
CHK={
"NC_037360.1":["Conspecific and congeneric: all 9 confident *Microdevario* calls are C33 (`results/gate3/candidate_NC_037360.1.txt`, GATE_3_SWEEP). Independent *M. kubotai* records PV357791–PV357794 (one study), AP011403 (RefSeq source), AB937094.","Reads: amplicon run SRR4556446 10 with C33, 0 with U33; WGS run SRR39401262 572 with C33, 1 with U33 (`results/gate3/reads_microdevario.txt`).","NUMT: no nuclear assembly in NCBI, so not checked. mtRF1a: **untested (no nuclear assembly)**."],
"NC_015546.1":["Same genus as *M. kubotai* (both species C33). No separate read check was done for *M. nanus*. The record is the RefSeq copy of AP011402.","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"CM165726.1":["This record is the mitochondrial contig of the assembly `GCA_056783665.1` (VGP fAlbGlo1 hap1).","Lineage evidence: independent *A. vulpes* mitogenome OZ556087.1 is C33; *A. vulpes* own run SRR22396674 has 469 reads C33 and 0 U33 (`results/gate3/reads_albula.txt`). No read test was done on this *A. glossodonta* record itself.","NUMT: 0 nuclear hits ≥200 bp and ≥80% in hap1 and in `GCA_019788955.1` (`results/gate3/numt_albula_hap1.txt`, `numt_albula_BYU.txt`); hap2 not searched; a planted 10%-diverged 2 kb copy was found. AP002973 is a named extra with a doubtful species label (`notes_albula_records.md`)."],
"NC_005796.1":["Albuliformes is uniformly C33 in the sweep (4 genomes, 1 of them panel, 1 a named extra). **No conspecific, read or NUMT check was done for *Pterothrissus* itself**; it is one 2004 record."],
"NC_013627.1":["Independent records: PX939991.1, PX754814.1, PZ284732.1 (all C33; `conspecific_b.txt`).","Reads: the record's own 31-mer matched 0 reads in SRR36227421 because that fish differs by one base nearby; with its own sequence, 272 reads C33 and 0 U33 (`reads_b_NC_013627.1_redo.txt`). The reads are from a different specimen (PX754814) than the record (A404).","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"NC_013436.1":["Independent records: PZ743035.1, PZ285086.1, PZ285087.1 (C33) and OR482581.1 (filed as *Stemonidium*, sequence is *S. sector*). OR546184.1, filed as *S. sector* with U33, is a frogfish and is set aside.","Reads: SRR22396609, 363 C33, 0 U33 (`reads_b_NC_013436.1.txt`); this is the specimen of OR482581.","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"PZ495457.1":["Single record. Reads: SRR38368693 is the record's own specimen: 419 C33, 0 U33 (`reads_b_PZ495457.1.txt`). No independent specimen.","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"NC_013628.2":["Independent records: PZ285083.1, PZ285084.1, PZ285088.1, PZ285090.1 (all C33).","Reads: SRR38094758 85 C33, 0 U33 (1 read with A) from specimen PZ285090, not the record (`reads_b_NC_013628.2.txt`).","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"CM163399.1":["Four *L. verrillii* records are C33 (this one, QB104173.1, PX940805.1, PV339947.1); other *Lycenchelys* species in the search (7 records) are U33 (`conspecific_b.txt`).","Reads: SRR37037826 622 C33, 1 U33, from voucher PX940805's specimen (`reads_b_CM163399.1.txt`).","NUMT: 0 nuclear hits ≥200 bp, ≥80%; an 82 kb contig (JBQQXK010000041.1) that is ≥99% mitochondrial over most of its length was excluded (`numt_lycenchelys.txt`)."],
"OZ407450.1":["Reads: RNA ERR14857142 35 G33, 0 other; DNA ERR14842678 55 G33, 0 other (`reads_b_OZ407450.1.txt`).","Nuclear copies: 322 hits ≥200 bp, ≥80% in hap1, 38 cover the site and **all read G** (`numt_craniella_zetl.txt`). No U33 read was seen, but the source copy of the record cannot be proven.","Same Sanger Aquatic Symbiosis pipeline as *C. infrequens* (BioProjects PRJEB107262, PRJEB96775); the whole mitogenomes are 97.1% identical; the trnY is identical."],
"OZ318138.1":["Reads: DNA ERR14835550 40 G33, 0 U33; RNA SRR33975184 66 G33, 0 U33 (separate study; species label not independently checked); RNA ERR14792883 2 G33 (`reads_b_OZ318138.1.txt`; two lines restored from the log after an overwrite).","Nuclear copies: 217 hits ≥200 bp, ≥80% in hap1, 4 cover the site and **all read G** (`numt_craniella_infr.txt`).","Relatives: *Cinachyrella* (3 records), *C. villosa* and *C. arb* are U33 (`posthoc_tetillidae`); *C. arb* own reads 216 U33, 0 other. The second GUA-like hit in all Tetillidae lies inside ND2 (NC_010198), so there is no second trnY."],
"PV357903.1":["Reads: SRR25465297 88 C33, 0 U33 (own specimen). *Coloconger cadenati* is U33. No second specimen."],
"PX804232.1":["Reads: SRR36182985 200 C33, 0 U33 (own specimen). Other *Ophichthus* records in the search are U33 (25 of 26; this one is the C33). No second specimen."],
"KY021065.1":["One record; no SRA runs for the species; no conspecific."],
"NC_082198.1":["One isolate (Taizhou): NC_082198.1 is the RefSeq copy of OR241493.1. Other *Toxocara* records with a confident call are U33 (13 of 15 U33, 2 are this species and its copy). No SRA runs."],
"KU726971.1":["Reads: SRR15202050 176 C33 and 0 U33 (1 G); SRR15202047 271 C33, 0 U33 (`reads_b_KU726971.1.txt`). So the base is real in the DNA of two populations.","**Why dropped:** the genome is multi-chromosome; chromosome II (KU726972.1) carries another GUA-loop tRNA-Tyr-like locus with U33, score 11.1 against 10.0 here (the gate minimum, one caller). Other *Globodera* are U33. The C33 locus may not be tRNA-Tyr.","NUMT: 10 nuclear hits, none over N33 (`numt_globodera.txt`)."],
"NC_047225.1":["**Why dropped:** RAD-seq runs (SRR9586127, SRR9586188) hold no reads matching either allele or the control 31-mers, so they say nothing (`reads_b_NC_047225.1.txt`). The only other record is the RefSeq source MK855100.1 (identical). 6 other confident *Bufotes* calls are U33. One caller, 4/5 stem."],
"CM003310.1":["Independent genome reads SRR26324223 (PRJNA1026095): all 32 TAA sites at depth 70–474, TAA 10,613 reads, TAT+TAC 23 (`roty_SRR26324223.txt`).","RNA-seq ERR744555 (PRJEB8325): 16 of 32 TAA sites covered, 55 reads, 0 TAT/TAC (`roty_ERR744555.txt`).","49 COX1 fragments from many labs carry internal TAA; **transl_table 14 on such records is filing practice, not evidence**.","TAG sites: COX1 404, COX3 239, CYTB 62 (`roty_tag_sites.txt`): only COX3 239 is a conserved-Tyr column; only COX1 404 is seen in RNA-seq.","**Why blocked:** the locked code counts 3 of 5 anticodon-stem pairs (score 12.8, with one second caller), so the locus audit blocks a new-clade claim.","**Stem-count finding (REPORT §3.6; not adopted).** PREREG does not say how pairs are counted; the locked code counts at model columns, where the aligner deletes one stem column and gives 3 of 5. Counted beside the loop, the arm above is 5 of 5 pairs (one G·U), as in *P. goodeyi* under the locked code. With 5 of 5 the genome would be *resolved* and meet §3 as a new clade. See `src/stem_direct_check.py`."],
"ON525593.1":["Weak audit: supported weak (10 supported sites in 4 genes; replication from *Radopholus similis*, a panel genome in the same family).","No COI or other mitochondrial record for *P. goodeyi* in GenBank besides this one; no SRA run or assembly.","ARAGORN reports a UUA-anticodon locus (15202–15260, 59 nt, 90% AT), tRNAscan-SE reports none; settled as not evidence (29 of 444 other nematode genomes have one).","The record is UNVERIFIED and linear; 18 of 35 TAA sites have p(Tyr) < 0.1."],
}
def seqof(acc):
    p=(G3U[acc]["path"] if acc=="ON525593.1" else "data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna"); return "".join(l.strip() for l in open(p) if l[0]!=">").upper()
PAIR={("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
def trny(acc,d):
    t=d["trny"]; s=seqof(acc); L=[]
    if t.get("status")!="called": return [f"Not called: {t.get('reason')} (score {t.get('score')}, anticodon {t.get('anticodon')})."]
    a,b=sorted((int(t["start"]),int(t["end"]))); x=s[a-1:b]; x=rc(x) if t["strand"]=="-" else x; x=x.replace("T","U"); g=int(t["g34_pos"]); j=(b-g) if t["strand"]=="-" else (g-a)   # 0-based index of G34 in the tRNA
    L.append(f"Location {t['start']}–{t['end']} ({t['strand']}), CM score {t['score']}, reason `{t['reason']}`, anticodon {t['anticodon']}, loop `{t['loop']}`, N33 = **{t['N33']}**, stem pairs {t['stem_pairs']}/5.")
    L.append(""); L.append("```"); L.append(x)
    loop=x[j-2:j+5]; s5=x[j-7:j-2]; s3=x[j+5:j+10]; pr="".join("|" if (s5[i],s3[4-i]) in PAIR and not (s5[i],s3[4-i]) in {("G","U"),("U","G")} else ("o" if (s5[i],s3[4-i]) in PAIR else "x") for i in range(5)) if len(s5)==5 and len(s3)==5 else "n/a"
    L.append(f"anticodon arm: 5'-{s5}-[{loop}]-{s3}-3'   pairs (5' to 3' of the 5' strand): {pr}  (| Watson-Crick, o G·U, x mismatch)"); L.append("```")
    if loop!=t["loop"]: L.append(f"WARNING: loop recomputed from coordinates ({loop}) differs from the caller's ({t['loop']}).")
    s2=d.get("second") or {}; L.append(f"Second callers: ARAGORN {len(s2.get('aragorn',[]))} Tyr locus/loci, tRNAscan-SE {len(s2.get('trnascan',[]))}. Verified GUA arms from the CM search: {t.get('n_ok')}.")
    return L
def build(acc):
    d=R[acc]; m=(G3U if acc=="ON525593.1" else G3)[acc]; u=d.get("uaa") or {}; z=K.n33(d,FE.get(acc) if acc!="ON525593.1" else None); L=[]
    L.append(f"# {m['species']} `{acc}`"); L.append(""); L.append(f"**Status: {STAT[acc]}.**"); L.append("")
    L.append("## Record"); L.append(""); L.append(f"- {m['phylum']}, {m['class']}, {m['order']}, {m['family']}. Provenance {m['provenance']}. Role: {m['role']}."); L.append(f"- Length {m['length']} bp; topology in the record: {m['topology']}. Ambiguous bases: {m['nN']}.")
    L.append(f"- UAA class: **{K.uaa_class(d)[0]}**. Internal TAA {u.get('n_taa_internal')} in {u.get('n_taa_genes')} genes, {u.get('n_taa_tyrcol')} at conserved-Tyr columns of {u.get('tyrcols')}; internal TAG {u.get('n_tag_internal')}; LLR {u.get('llr_Y')}; Codetta TAA = {(d.get('codetta') or {}).get('TAA')}; aligned codons {u.get('aligned_codons')}, genes {u.get('n_genes_hit')}.")
    L.append(f"- Conserved-Tyr columns: TAT {u.get('tat')}, TAC {u.get('tac')}, TAA {u.get('n_taa_tyrcol')}, TAG {u.get('tag_tyrcol')}, other {u.get('tyrcol_other')}.")
    L.append(""); L.append("## tRNA-Tyr(GUA)"); L.append(""); L+=trny(acc,d); L.append(f"- Locus status (PREREG audit): **{z['status']}** ({z['reason']}); pass at 0/5/10 nt: {z['pass0']}/{z['pass5']}/{z['pass10']}.")
    L.append(""); L.append("## Protein-coding genes"); L.append("")
    r=dict(m,path=m.get("path") or None); 
    if not r["path"]: r.pop("path")
    cm,tag,fr=P.cols(dict(r)); ta=collections.Counter(g for (g,c),v in cm.items() if v=="TAA"); tg=collections.Counter(g for g,c in tag)
    f=d.get("first") or {}; tf=d.get("tagfirst") or {}; L.append("| Gene | aligned domain | reading frames | internal TAA | internal TAG | terminal codon (first in-frame rule) |"); L.append("|---|---|---|---|---|---|")
    genes=sorted({p[3] for p in d["pcg"]}|set(f))
    for g in genes: L.append(f"| {g} | {'yes' if any(p[3]==g for p in d['pcg']) else 'no'} | {len(fr.get(g,[]))} | {ta.get(g,0)} | {tg.get(g,0)} | {f.get(g,'')} |")
    L.append(""); L.append("Terminal codons: "+", ".join(f"{k} {v}" for k,v in sorted(collections.Counter(f.values()).items())))
    if acc in("CM003310.1","ON525593.1"):
        L.append(""); L.append("## Conserved-Tyr excerpts (each internal TAA; relatives = 10 closest Gate 3 Stop species)"); L.append("")
        S=[q for q in G3.values() if q["status"]=="included" and q["provenance"]!="named_extra" and q["uaa"]=="Stop"]; rr=dict(m,uaa="weak"); 
        if acc=="ON525593.1": rr["path"]=m["path"]
        rel=P.closest(rr,S); L.append("Relatives: "+", ".join(x["species"]+f" ({k})" for k,x in rel)); L.append(""); L.append("| gene, column | relatives' codons | relatives' residues |"); L.append("|---|---|---|")
        for (g,c) in sorted(k for k,v in cm.items() if v=="TAA"):
            v=[P.cols(x)[0].get((g,c)) for _,x in rel]; al=collections.Counter(AA.get(x,"?") if x else "-" for x in v); L.append(f"| {g} {c} | {' '.join(x or '---' for x in v)} | {', '.join(f'{k} {n}' for k,n in al.most_common())} |")
    L.append(""); L.append("## mtRF1a"); L.append("")
    if acc in MT:
        L.append("| source | length | 203 | 216 | loop 196–222 | label | best bacterial match |"); L.append("|---|---|---|---|---|---|---|")
        for s in MT[acc]:
            for q in sorted(LOOP[s],key=lambda q:-int(q["ident_to_human"]))[:2]: L.append(f"| {s} | {q['length']} | {q['p203']} | {q['p216']} | `{q['loop']}` | {q['label']} | {q['bact_best']} |")
        L.append(""); L.append("Paralog check: best human match MTRF1L by blastp (`results/gate4/labels.tsv`), placed in `results/gate4/tree.nwk`; not bacterial (E. coli RF1 identity ≤48%; for the contaminant that failed see *Hoplolaimus columbus* in GATE_4)." if acc!="CM003310.1" else "Paralog check: all *Rotylenchulus* models are MTRF1L by best human match (`labels.tsv`) and cluster with nematode RF1 in `results/gate4/tree.nwk`; E. coli RF1 identity 35–38%. The species has 8 genome and 40 TSA models, all Q203/S216.")
    else: L.append("Untested: no nuclear assembly for the species" if acc not in("KU726971.1","NC_082198.1","NC_047225.1","KY021065.1") else "Not run: candidate dropped or unconfirmed (Gate 4 scope).")
    L.append(""); L.append("## Checks and their sources"); L.append(""); L+=[f"- {x}" for x in CHK[acc]]
    return "\n".join(L)+"\n"
idx=["# Outlier dossiers","","Built by `src/dossiers.py` from committed outputs. 17 non-U33 Stop candidates (9 C33 supported, 4 C33 unconfirmed, 1 C33 dropped, 2 G33 supported, 1 G33 dropped), plus *Rotylenchulus* and *Pratylenchus goodeyi* as sensitivity cases. The status words are post-hoc; the locked rules count all 17.","","| Dossier | Species | Family | Base | Status |","|---|---|---|---|---|"]
for a in WANT+["ON525593.1"]:
    m=(G3U if a=="ON525593.1" else G3)[a]; fn=f"{a}_{m['species'].replace(' ','_')}.md"; open(D+fn,"w").write(build(a)); t=R[a]["trny"]; idx.append(f"| [{a}]({fn}) | *{m['species']}* | {m['family']} | {t.get('N33')} | {STAT[a]} |")
open(D+"README.md","w").write("\n".join(idx)+"\n"); print("dossiers written",len(WANT)+1)
