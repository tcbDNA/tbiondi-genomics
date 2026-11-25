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
TIER={r["acc"]:r for r in csv.DictReader(open("results/cold_review/tiers.tsv"),delimiter="\t")}
STAT={a:f"{r['tier']} ({r['base']}; NCBI mitochondrial code {r['ncbi_mt_code']}; post-hoc tier rule, REPORT 7.2.1)" for a,r in TIER.items()}
STAT["CM003310.1"]="UNCONFIRMED SINGLE-GENOME CANDIDATE (PREREG 3 label; blocked by the locked trnY stem count)"
STAT["ON525593.1"]="SENSITIVITY (supported weak; UNVERIFIED record, not confirmed)"
CHK={
"NC_037360.1":["**Distinct genomes:** the 9 confident *Microdevario* calls reported in Gate 3 are 7 distinct genomes from 3 sources: PV357791-PV357794 (one study, 4 isolates), AP011402 and AP011403 (with RefSeq copies NC_015546 and NC_037360, identical by checksum), AB937094.","Reads: amplicon run SRR4556446 10 C33, 0 U33; WGS run SRR39401262 572 C33, 1 U33 (`results/gate3/reads_microdevario.txt`).","NUMT: no nuclear assembly in NCBI, so not checked. mtRF1a: **untested (no nuclear assembly)**."],
"NC_015546.1":["Tier 3: no same-species evidence (the record is the RefSeq copy of AP011402, identical by checksum); the genus is uniformly C33 and *M. kubotai* is Tier 1. No separate read check was done for *M. nanus*.","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"CM165726.1":["**PREREG T6 read check (run after the cold review).** The record's BioProject PRJNA1427880 has no reads. Same specimen (fAlbGlo1, BioSample SAMN55933885), PacBio run SRR38145622 (PRJNA1454682): 32 reads aligned across N33, **32 C33, 0 U33** (unchanged down to 85% identity). Independent specimen (HI-2016, 2016/17; PRJNA668352), RNA-seq SRR13036352: **1,797 C33, 1 U33** at >=98%, 1,821 and 1 at >=90%. `results/cold_review/t6_reads.txt`. That fish is *A. glossodonta*, not *A. virgata*: the majority CYTB consensus of its reads is identical to CM165726.1 over 1,110 bp (COX1 differs at 1 of 1,515), while *A. virgata* CYTB (KJ910045.1, the only *A. virgata* record over 350 bp) is 82.7% identical; 279,991 of its reads carry *A. glossodonta*-diagnostic CYTB 25-mers and 149 *A. virgata*-diagnostic ones, 148 of which carry both (`src/albula_reads_id.py`, `src/albula_kmers.py`).","This record is the mitochondrial contig of `GCA_056783665.1` (VGP fAlbGlo1 hap1).","**Panel evidence, not counted:** *A. vulpes* 469 to 0 (SRR22396674) is a calibration-panel species. OZ556087, filed as *A. vulpes* (BioSample SAMEA123450433, Indian Ocean), is C33 but matches *A. glossodonta* (CYTB 99.8% to CM165726.1, 100% to *A. glossodonta* OL542768; the panel *A. vulpes* NC_082992 is 92.5% identical to CM165726.1), so its label is doubtful. It is counted nowhere.","NUMT: 0 nuclear hits >=200 bp and >=80% in hap1 and in `GCA_019788955.1` (`numt_albula_hap1.txt`, `numt_albula_BYU.txt`); hap2 not searched; a planted 10%-diverged 2 kb copy was found. AP002973 is a named extra with a doubtful species label.","C33 in *Albula* was known before the PREREG lock (GATE_2 A8)."],
"NC_005796.1":["**No SRA data** for *Pterothrissus*. Tier 3: no same-species evidence; every confident Albulidae call is C33 and *A. glossodonta* is Tier 1. No conspecific, read or NUMT check was done for this 2004 record."],
"NC_013627.1":["Independent records: PX939991.1, PX754814.1, PZ284732.1 (all C33; `conspecific_b.txt`).","Reads: the record's own 31-mer matched 0 reads in SRR36227421 because that fish differs by one base nearby; with its own sequence, 272 reads C33 and 0 U33 (`reads_b_NC_013627.1_redo.txt`). The reads are from a different specimen (PX754814) than the record (A404).","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"NC_013436.1":["Independent records: PZ743035.1, PZ285086.1, PZ285087.1 (C33) and OR482581.1 (filed as *Stemonidium*, sequence is *S. sector*). OR546184.1, filed as *S. sector* with U33, is a frogfish and is set aside.","Reads: SRR22396609, 363 C33, 0 U33 (`reads_b_NC_013436.1.txt`); this is the specimen of OR482581.","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"PZ495457.1":["Single record. Reads: SRR38368693 is the record's own specimen: 419 C33, 0 U33 (`reads_b_PZ495457.1.txt`). No independent specimen.","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"NC_013628.2":["Independent records: PZ285083.1, PZ285084.1, PZ285088.1, PZ285090.1 (all C33).","Reads: SRR38094758 85 C33, 0 U33 (1 read with A) from specimen PZ285090, not the record (`reads_b_NC_013628.2.txt`).","NUMT not checked. mtRF1a: **untested (no nuclear assembly)**."],
"CM163399.1":["Four *L. verrillii* records are C33 (this one, QB104173.1, PX940805.1, PV339947.1); other *Lycenchelys* species in the search (7 records) are U33 (`conspecific_b.txt`).","Reads: SRR37037826 622 C33, 1 U33, from voucher PX940805's specimen (`reads_b_CM163399.1.txt`).","NUMT: 0 nuclear hits >=200 bp and >=80%. Contig JBQQXK010000041.1, excluded, is **90,807 bp** (not 82 kb) of tandem mitochondrial copies about 19.6-19.8 kb apart; **all 5 tRNA-Tyr copies on it are C33** (6250, 25864, 45640, 65254, 85032)."],
"OZ407450.1":["Reads: own-specimen DNA ERR14842678 55 G33, 0 other; RNA ERR14857142 35 G33, 0 other (`reads_b_OZ407450.1.txt`). Tier 2.","Nuclear copies: 322 hits >=200 bp, >=80% in hap1, 38 cover the site and **all read G**. The source copy of the record cannot be proven.","Same Sanger Aquatic Symbiosis pipeline as *C. infrequens* (PRJEB107262, PRJEB96775); whole mitogenomes 97.1% identical; trnY identical."],
"OZ318138.1":["Independent specimen: RNA-seq SRR33975184 (PRJNA1276125), label confirmed by 1,086,464 *C. infrequens*-diagnostic 31-mers against 822 *C. zetlandica*; aligned reads across N33 **49-55 G33, 0 other**. Own DNA ERR14835550 40 G33, 0 U33. Tier 1.","Nuclear copies: 217 hits >=200 bp, >=80% in hap1, 4 cover the site and **all read G**.","Relatives: *Cinachyrella* (3 records), *C. villosa*, *C. arb* are U33; *C. arb* own reads 216 U33, 0 other."],
"PV357903.1":["Reads: SRR25465297 88 C33, 0 U33 (own specimen); Tier 2. *Coloconger cadenati* is U33. No second specimen.","mtRF1a: **skipped by the Gate 4 plan**; no nuclear assembly exists either."],
"PX804232.1":["Reads: SRR36182985 200 C33, 0 U33 (own specimen); Tier 2. Other *Ophichthus* records in the search are U33 (25 of 26). No second specimen.","mtRF1a: **skipped by the Gate 4 plan**; no nuclear assembly exists either."],
"KY021065.1":["One record; no SRA runs for the species; no conspecific. Unconfirmed under the tier rule; the gene tree joins it with *O. asakusae* (support 1.0), which the rule does not count.","mtRF1a: skipped by the Gate 4 plan."],
"NC_082198.1":["One isolate (Taizhou): NC_082198.1 is the RefSeq copy of OR241493.1. 13 of 15 other *Toxocara* records with a confident call are U33; the two C33 are this species and its copy. No SRA runs. Unconfirmed; code 5.","mtRF1a: skipped by the Gate 4 plan."],
"KU726971.1":["Reads: SRR15202050 176 C33, 0 U33 (1 G); SRR15202047 271 C33, 0 U33 (`reads_b_KU726971.1.txt`).","**Dropped, confirmed.** Chromosome II (KU726972.1) has a tRNA-Tyr(GUA) with U33 at 6699-6752, placed at exactly those coordinates by Infernal and ARAGORN, in a normal tRNA cluster. Here on chromosome I the two callers disagree (16851-16905 against 16830-16905; score 10.0).","mtRF1a: **skipped** (candidate dropped), although a nuclear assembly exists.","NUMT: 10 nuclear hits, none over N33."],
"NC_047225.1":["**Unconfirmed** (no longer \"dropped\": the tier rule has no positive evidence against it). RAD-seq runs (SRR9586127, SRR9586188) hold no reads matching either allele or the control 31-mers. The only other record is the RefSeq source MK855100.1 (identical). 6 other confident *Bufotes* calls are U33. One caller, 4/5 stem.","mtRF1a: skipped by the Gate 4 plan."],
"CM003310.1":["**Classification:** PREREG 3's \"unconfirmed single-genome candidate\"; the locked trnY stem count (3 of 5 at model columns; score 12.8) blocks it as a new clade. UAA evidence does not depend on the tRNA.","Independent genome reads SRR26324223 (PRJNA1026095): all 32 TAA sites at depth 70-474, TAA 10,613 reads, TAT+TAC 23. RNA-seq ERR744555 (PRJEB8325): 16 of 32 sites covered, 55 reads, 0 TAT/TAC.","**tRNA-Tyr boundaries.** The covariance model's call (22684-22632, 53 nt) has 9 nt before the anticodon arm and 27 after. ARAGORN's (22704-22643, 62 nt) has 29 and 16, like *P. goodeyi* (26/14) and *Radopholus* (27/14); it is 89% identical to *P. goodeyi* (51/57) and 86% to *Radopholus* (50/58), against 74% for *M. incognita* and 69% for *P. vulnus*. The model boundaries look wrong.","**Stem count (REPORT 7.2.2).** Counted beside the loop the arm pairs 5 of 5 (one G·U). GATE_2 A13 already recorded, in the lock commit, that a contiguous-arm recount disagreed in 4 of 305 calls. The locked count is kept.","**COX1 fragments:** 89 of 92 *Rotylenchulus* records have internal UAA (86 with two or more), in 5 named species and 7 locality-labelled isolates. **NCBI Taxonomy assigns mitochondrial code 14 to the genus**, so `transl_table` 14 is not independent evidence.","TAG sites (`roty_tag_sites.txt`; tylenchid Stop relatives rechecked): COX3 239 is Tyr in 18 of 18 readable; COX1 404 and CYTB 62 in 0 of 18. Only COX1 404 is seen in RNA-seq.","NUMTs in `GCA_001026735`: 3 nuclear-like segments (1,013, 658 and 269 bp), none over trnY; weak (N50 6 kb).","Termination (TSA, PRJNA286314): COX1 transcripts run past an in-frame UAA (15681-15683) and end near 15591, after a tandem UAG (15612-15617); COX3 transcripts pass its first UAA. Unresolved.","Shared internal UAA columns with *P. goodeyi*: 9 (7 at conserved-Tyr columns). COX1 324 is Tyr in 9 of 11 tylenchid Stop relatives."],
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
SECOND_TXT={"CM165726.1":"tRNA-Thr(UGU), opposite strand; one-base slip","NC_005796.1":"tRNA-Thr(UGU), opposite strand; one-base slip","NC_037360.1":"tRNA-Thr(UGU), opposite strand; one-base slip","NC_015546.1":"tRNA-Thr(UGU), opposite strand; one-base slip","CM163399.1":"tRNA-Thr(UGU), opposite strand; one-base slip","PV357903.1":"tRNA-Thr(UGU), opposite strand; one-base slip","NC_047225.1":"tRNA-Thr(UGU), opposite strand; one-base slip","OZ407450.1":"low-score hit inside the ND2 reading frame (455 aa), no second caller","OZ318138.1":"low-score hit inside the ND2 reading frame (475 aa), no second caller","NC_013627.1":"ARAGORN-only call 12891-12946, antisense inside ND5; no model hit","KY021065.1":"ARAGORN-only call 10915-10979, antisense inside ND4; no model hit","KU726971.1":"a real second tRNA-Tyr on chromosome II (KU726972.1, 6699-6752, U33), by Infernal and ARAGORN"}
import trny_extract as TX
def second(acc,d):
    L=["","### Second GUA loci",""]; t=d["trny"]; p=(G3U[acc]["path"] if acc=="ON525593.1" else "data/sweep/fa/"+str(sum(map(ord,acc))%200).zfill(3)+f"/{acc}.fna")
    h,e=TX.run(p,t.get("reason")=="sens") if t.get("status")=="called" else ([],None); ar=d["trna"]["aragorn"]
    for x in sorted((x for x in (h or []) if x["anticodon"]=="GUA" and x["loop_ok"] and x["stem_pairs"]>=3),key=lambda x:-x["score"]):
        u,v=sorted((x["start"],x["end"])); a0,b0=sorted((int(t["start"]),int(t["end"])))
        if u==a0 and v==b0: continue
        L.append(f"- Model hit {x['start']}-{x['end']} ({x['strand']}), score {x['score']}, loop `{x['loop']}`, N33 {x['N33']}; ARAGORN genes overlapping it: "+(", ".join(f"{q[3]}({q[4]})" for q in ar if q[0]<=v and q[1]>=u) or "none")+".")
    n=len([1 for q in ar if q[4]=="GTA"]); L.append(f"- ARAGORN calls {n} Tyr(GTA) gene(s) in the record; tRNAscan-SE {len((d.get('second') or {}).get('trnascan',[]))}.")
    if acc in SECOND_TXT: L.append(f"- **Interpretation:** {SECOND_TXT[acc]} (REPORT 7.2.7).")
    return L
def build(acc):
    d=R[acc]; m=(G3U if acc=="ON525593.1" else G3)[acc]; u=d.get("uaa") or {}; z=K.n33(d,FE.get(acc) if acc!="ON525593.1" else None); L=[]
    L.append(f"# {m['species']} `{acc}`"); L.append(""); L.append(f"**Status: {STAT[acc]}.**"); L.append("")
    L.append("## Record"); L.append(""); L.append(f"- {m['phylum']}, {m['class']}, {m['order']}, {m['family']}. Provenance {m['provenance']}. Role: {m['role']}."); L.append(f"- Length {m['length']} bp; topology in the record: {m['topology']}. Ambiguous bases: {m['nN']}.")
    L.append(f"- UAA class: **{K.uaa_class(d)[0]}**. Internal TAA {u.get('n_taa_internal')} in {u.get('n_taa_genes')} genes, {u.get('n_taa_tyrcol')} at conserved-Tyr columns of {u.get('tyrcols')}; internal TAG {u.get('n_tag_internal')}; LLR {u.get('llr_Y')}; Codetta TAA = {(d.get('codetta') or {}).get('TAA')}; aligned codons {u.get('aligned_codons')}, genes {u.get('n_genes_hit')}.")
    L.append(f"- Conserved-Tyr columns: TAT {u.get('tat')}, TAC {u.get('tac')}, TAA {u.get('n_taa_tyrcol')}, TAG {u.get('tag_tyrcol')}, other {u.get('tyrcol_other')}.")
    L.append(""); L.append("## tRNA-Tyr(GUA)"); L.append(""); L+=trny(acc,d); L.append(f"- Locus status (PREREG audit): **{z['status']}** ({z['reason']}); pass at 0/5/10 nt: {z['pass0']}/{z['pass5']}/{z['pass10']}.")
    L+=second(acc,d)
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
    else: L.append("Untested: no nuclear assembly for the species." if acc not in("KU726971.1","NC_082198.1","NC_047225.1","KY021065.1","PV357903.1","PX804232.1") else "Not run: skipped by the Gate 4 plan (see the checks below).")
    L.append(""); L.append("## Checks and their sources"); L.append(""); L+=[f"- {x}" for x in CHK[acc]]
    return "\n".join(L)+"\n"
idx=["# Outlier dossiers","","Built by `src/dossiers.py` from committed outputs. 17 non-U33 Stop candidates (post-hoc tier rule, REPORT 7.2.1: C33: Tier 1, 6; Tier 2, 3; Tier 3, 2; unconfirmed, 2; dropped, 1; G33: Tier 1, 1; Tier 2, 1; unconfirmed, 1), plus *Rotylenchulus* (unconfirmed single-genome candidate) and *Pratylenchus goodeyi* (sensitivity). The locked rules count all 17.","","| Dossier | Species | Family | Base | Status |","|---|---|---|---|---|"]
for a in WANT+["ON525593.1"]:
    m=(G3U if a=="ON525593.1" else G3)[a]; fn=f"{a}_{m['species'].replace(' ','_')}.md"; open(D+fn,"w").write(build(a)); t=R[a]["trny"]; idx.append(f"| [{a}]({fn}) | *{m['species']}* | {m['family']} | {t.get('N33')} | {STAT[a]} |")
open(D+"README.md","w").write("\n".join(idx)+"\n"); print("dossiers written",len(WANT)+1)
