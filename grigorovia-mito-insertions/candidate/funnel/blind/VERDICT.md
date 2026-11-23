# VERDICT: frameshifts in three budding-yeast mitochondrial genomes

## Short answer

| Genome | Frameshifts in the 7 core genes | Direction | Verdict |
|---|---|---|---|
| genome_1 | **0** | none | No frameshifts, so there is nothing to classify. All 7 genes read through in one frame. |
| genome_2 | **4** (cox1 ×1, cob ×3) | all **+1** | **Systematic.** Every site is the same motif. An extra A/T sits right before a **GGC** glycine codon, so the original frame reads AGG or TGG and the shifted frame reads GGC. All four are at glycines conserved in both references and in the other two genomes. GGC is used in frame nowhere else in genome_2. |
| genome_3 | **5** (cob, cox2, cox3 ×2, atp6) | all **+1** | **Systematic, with a different motif.** At every site the extra base is the **C of a CTG codon** in the original frame. Skipping it makes the next codon TGT or TGC (Cys). CTG is never used in frame anywhere else in genome_3. The pattern is looser than in genome_2: the sites are spread over 4 genes, the residue after the shift is not conserved, and 2 of the 5 exact positions rely on comparison with genomes 1 and 2 rather than on the references. |

How to read the coordinates: they are 1-based positions on the forward strand of each `genome_N.fna`. "Coordinate" is the position of the extra base. Codons are written 5'→3' on the coding strand, and genome_2's sites are all on the minus strand. "Original reading frame" means the frame being read just before the shift. Translation uses NCBI code 3 (TGA=W, ATA=M, CTN=T).

---

## 1. What was done (overview)

1. **Reference names.** `reference_proteins.faa` has three entries all called `refA_cob`. The first is plain cob (385 aa). The second and third are cob exons fused to two different intron ORFs (517 and 423 aa). They were renamed `refA_cob_v1/_v2/_v3`. Each genome was searched with all 16 references (refA and refB versions of cox1, cox2, cox3, cob, atp6, atp8 and atp9, plus the two fusions).
2. **miniprot `-T 3`** was run against each genome as given (linear).
3. **Circular handling.** These are circular mtDNAs, and two genes cross the sequence end: **genome_1 cox1** and **genome_3 cob**. In the linear run genome_1 cox1 aligned only from aa 235 to 534. Genome_3 cob gave a spurious "+2" call, or two +1 calls, next to the break. I therefore re-ran on genomes with their first 8 kb appended to the end ("circ"). I also re-ran on genomes rotated to start at their midpoint, as an independent check. All final calls come from the circ runs, with coordinates folded back onto the original sequence.
4. **Frameshift penalty sweep.** The default is `-F 23`; I also ran `-F 6`, `-F 12`, `-F 15` and `-F 35`.
5. **Event extraction.** Frameshifts were taken from the `cg:Z:` CIGAR. `F` means extra genome nucleotides (a +1 for one extra base). `G` means a residue aligned to only 1–2 nt (a −1 for a 2-nt codon).
6. **Placement windows.** miniprot's exact frameshift position is only one of several equally good choices. For each event I tested every placement within a few codons and recorded all placements that give the best identity to each reference.
7. **Cross-genome comparison.** At each site I printed the codons that all three genomes align to the same reference residues. This gives a nucleotide-level view of what the "extra" base is relative to related genomes.
8. **Pattern tests.** I counted in-frame codon usage for all 7 genes in each genome. I located every GGC and CTG trinucleotide relative to the reading frame, checked start and stop codons at the gene ends, and checked how far the unshifted frame runs before hitting a stop.

All the exact commands are in section 7.

## 2. Genes found (circ run, best reference per gene)

All 7 core genes were found in full length in all three genomes. Every gene ends at TAA within 0–2 codons of the end of the alignment, so no frameshift is hidden at a C-terminus.

| Gene | genome_1 (32,558 bp, GC 23.1%) | genome_2 (29,214 bp, GC 23.8%) | genome_3 (31,230 bp, GC 22.6%) |
|---|---|---|---|
| cox1 | + 29788→2286 (crosses the end; 2 introns, 2087 and 1341 nt) | − 23355–24942 (no intron) | + 14822–18210 (2 introns, about 377 and 1422 nt) |
| atp8 | + 2611–2754 | − 1758–1901 | + 18680–18823 |
| atp6 | + 2992–3768 | − 936–1709 | + 19221–20007 |
| cox2 | + 7690–8436 | + 2624–3376 | + 26528–27278 |
| cob | − 23227–24381 (no intron) | − 25750–26904 (no intron) | + 29429→1036 (crosses the end; 1 intron, 1688 nt; exon 2 starts at about 91) |
| atp9 | − 24884–25108 | + 27697–27924 | + 12026–12253 |
| cox3 | + 26691–27497 | − 22079–22885 | + 9449–10257 |

---

## 3. genome_1: no frameshifts

* At default settings, and at `-F 12`, `-F 15` and `-F 35`, no alignment of any core gene contains an F or G operation. The circular, rotated and linear runs all agree.
* At `-F 6`, which is far below the default, the only calls are at alignment ends: the last residue of atp9 (`75M2F1M`) and the divergent N-terminus of cox2 (`8M1G3M1G1G…`, refA only). Both are low-similarity terminal regions, only one reference supports each, and they disappear at `-F 12`. I treat them as artifacts.
* The motif codons that shift in genomes 2 and 3 are absent in frame here as well: AGG 0, GGC 0, CTG 0, TGG 2.

## 4. genome_2: 4 frameshifts, all +1

All four sites are on the minus strand. They are stable across `-F 12` to `-F 35`, and refA and refB call them at identical positions.

| # | Gene | Coordinate of extra base | Dir. | Extra base | Codons in original frame (shift codon in [ ]) | Read after the shift | Ref residue right after the shift (refA / refB) |
|---|---|---|---|---|---|---|---|
| 1 | cox1 | **24711** (24712 equally good, since AA; miniprot: 24711) | +1 | A | TTA ATT GGA **[AGG]** CTT TGG TAA | **GGC** TTT GGT → G F G | **G78 / G78** |
| 2 | cob | **26496** (unique; miniprot: 26496) | +1 | A | TGT GTA TAT **[AGG]** CCA AAT GTC | **GGC** CAA ATG → G Q M | **G137 / G137** |
| 3 | cob | **26477** (unique; miniprot: 26477) | +1 | T | TCA CAT TGA **[TGG]** CGC ACT AGT (frame as set by shift #2) | **GGC** GCA CTA → G A T | **G143 / G143** |
| 4 | cob | **26263** (26264 and 26266 equally good; **miniprot called 26266**) | +1 | A | GGA ATT ACT **[AGG]** CAA TTT AGA | **GGC** AAT TTA → G N L | **G214 / G214**. With miniprot's placement at 26266 (ATT [A]CT AGG → CTA GGC) the residue after is **T213 / T213**. ACT, ACA and CTA are all Thr, so this choice cannot be resolved from the protein. |

If the ribosome stayed in the original frame, it would hit TAA or TAG after 4, 10, 7 and 5 codons at sites 1–4.

**Why this looks systematic.**
* **Same direction:** 4 of 4 are +1, and each is a single extra base.
* **Same motif:** the original-frame codon at the shift is AGG (3 sites) or TGG (1 site), followed by C. Skipping its first base always gives **GGC**.
* **Same residue:** all four GGCs encode a **glycine conserved in refA, refB, genome_1 and genome_3**. At these four positions (cox1 78, cob 137, 143, 214) genome_1 uses GGA/GGT/GGT/GGT and genome_3 uses GGT/GGA/GGA/GGA. Only genome_2 has "(A/T)+GGC" there.
* **Exclusive:** GGC occurs as an in-frame codon in genome_2 **only at these 4 sites**. The only other GGC in the 7 genes (26614) is out of frame. AGG and TGG never occur as in-frame codons in genome_2. So every in-frame (A/T)GGC is a frameshift site, and every frameshift site is one.
* **Not explained by chance:** apart from the sites, (A/T)GGC occurs about once per 5,400 coding positions. Four random insertion sites all landing on it is vanishingly unlikely (on the order of 10⁻¹³).
* **Not a homopolymer pattern:** only site 1 is in a 2-base run (AA). This is not what homopolymer sequencing errors look like.
* **Clustered:** three of the four sites are within 77 codons of cob (aa 137, 143, 214).

This is what one would expect from a recurring "+1 at AGG/TGG-C" signal, such as ribosomes skipping a codon they cannot decode well. It is equally what a deliberately placed motif would look like. Sequence data alone cannot tell programmed translational frameshifting apart from an engineered or systematic artifact, but it does rule out random scattered errors.

## 5. genome_3: 5 frameshifts, all +1

All five sites are on the plus strand. They are stable across `-F 12` to `-F 35`.

| # | Gene | Coordinate of extra base | Dir. | Extra base | Codons in original frame (shift codon in [ ]) | Read after the shift | Ref residue right after the shift (refA / refB) |
|---|---|---|---|---|---|---|---|
| 1 | cob (exon 2) | **121** (120–124 equally good; miniprot: 121 with both refs) | +1 | C | TTA AGA TAC **[CTG]** TCA TGC TAA | **TGT** CAT GCT → C H A | **L81 / M81** |
| 2 | cox2 | **27197** (unique; miniprot: 27197 with both refs) | +1 | C | AGT GAA TTG **[CTG]** CGG ATC CGG | **TGC** GGA TCC → C G S | **C225 / C226** |
| 3 | cox3 | **9815** (9814 equally good for refA; refB window 9811–9818; miniprot: refA 9815, refB 9818) | +1 | C | TTA GGA GCT **[CTG]** CTG ACC TCC | **TGC** TGA CCT → C W P | **C123 / V125**. With refB's placement at 9818 the residue after is W126. |
| 4 | cox3 | **10197** (preferred by comparison with genomes 1 and 2). **miniprot called 10200** with both refs. Windows that fit the references equally well: refA 10199–10201, refB 10196–10201. | +1 | C | GTA TTA TAT **[CTG]** CCA TGT CTT | **TGC** CAT GTC → C H V | **T250 / G252**. With miniprot's placement (TAT CTG [C]CA → CAT) the residue after is **H251 / H253**. |
| 5 | atp6 | **19803** (preferred by comparison with genomes 1 and 2). **miniprot called 19800** (an A; refB only). Window that fits the reference equally well: 19800–19801. | +1 | C (A under miniprot's placement) | GTT AAT AAT **[CTG]** TGC AGG TCA | **TGT** GCA GGT → C A G | **L192 / L192**. With miniprot's placement (AAT [A]AT CTG → ATC TGT) the residue after is **I191 / I191**. |

Notes on these sites:
* refA_atp6 failed to map to genome_3 under any setting I tried (default, `-M 0`, `-k 5`, `-n 1`, `-S`; reason unclear). Site 5 therefore rests on the refB alignment. refA has the same residues around it (refA 185–198 `LRLGSNILAGHLLM`, refB `LRLSSNILAGHLLM`).
* In the linear-file run, cob site 1 appeared as a 2-nt "+2" event at 120–121 (refA) or as two +1 events at 99 and 121 (refB). That happened because exon 1 sits at the far end of the file. The circular run resolves it to a single +1 at 121. That site lies 10 codons into exon 2, after 6 codons that match refB exactly (G W L L R Y), so it is not a splice-junction artifact.
* If the ribosome stayed in the original frame, it would hit a stop after 4, 11, 11, 24 and 16 codons at sites 1–5.

**Why the exact position at sites 4 and 5 departs from miniprot's.** The references have T (= CTG) at cox3 aa250 and I at atp6 aa191, so miniprot shifts one codon later or earlier there. The two related genomes settle it at nucleotide level, using the cross-genome codon table in `work/run/compare_sites.txt`:
* **atp6:** genome_1 and genome_2 read `AAT AAT TGT GCA GGT`. genome_3 reads `AAT AAT C TGT GCA GGT`. Removing the C at 19803 reproduces genome_1 exactly. Any other choice leaves mismatches.
* **cox3 #2:** genome_1 reads `TAT TGT CAT GTA`. genome_3 minus the C at 10197 reads `TAT TGC CAT GTC`, with only synonymous differences. miniprot's choice (10200) leaves 4 mismatches.

**Why this looks systematic, and why less tightly than genome_2.**
* **Same direction:** 5 of 5 are +1, each a single extra base.
* **Same base:** the extra base is **C in all 5** under the placements above. Even taking miniprot's raw positions, it is C in 4 of 5, in a genome whose genes are only about 12% C.
* **Same codon:** at every site the original-frame codon at the shift is **CTG**, and the extra base is its C. Using miniprot's raw positions, CTG is still the shift codon at 3 sites (cob, cox2, cox3 #1). At the other two it is the codon right before (cox3 #2) or right after (atp6) miniprot's extra base. After the skip the ribosome always reads **TGT or TGC (Cys)**.
* **Exclusive:** CTG is **never used as an in-frame codon** in genome_3 except at these sites. About 51 other CTG trinucleotides occur inside the genes, all of them out of frame. Other CTN codons are read normally (CTA occurs 13 times in frame), so the effect is specific to CTG. Genomes 1 and 2 also never use CTG in frame.
* **Not explained by chance:** CTG starts about 1% of coding positions. Five random insertions all landing on a CTG C within their windows is on the order of 10⁻⁷ or less.
* **Relative to the other genomes:** the extra C is an insertion compared with genome_1 and genome_2. Removing it restores their codons at atp6 exactly, and with only synonymous differences at cox3 #2.
* **Looser than genome_2 in three ways:**
  * The sites are spread over 4 genes at unrelated positions: cob aa81, cox2 aa225, cox3 aa123 and aa251, atp6 aa192.
  * The Cys that follows the shift is **not** a conserved residue. The references have L/M, C/C, C/V, T/G and L at those positions, and genome_1 has I, C, V, C and C.
  * At 2 of the 5 sites the exact placement comes from the other genomes rather than from the references.
  * Two of the five sites (cox2, cox3 #1) also sit in unusually GC-rich 81-nt windows (42% GC against 23% genome-wide). The other three do not (25–36%).

Overall, genome_3's frameshifts recur at one specific sequence motif (an original-frame CTG whose C is skipped), so I classify them as **systematic**. The evidence is weaker than for genome_2 because the residue after the shift is not conserved and the sites are scattered across genes rather than clustered.

The two genomes use different motifs: AGG/TGG→GGC in genome_2 and CTG→TGY in genome_3. None of the nine sites falls at a homologous position in another genome, and genome_1 reads straight through at every one of them.

---

## 6. Negative results, sensitivity and caveats

* **Penalty sweep.** The 4 calls in genome_2 and the 5 in genome_3 are identical at `-F 12`, `-F 15`, `-F 23` (default) and `-F 35`. At `-F 6` (and for a few refB alignments in genome_3 at `-F 12`), extra F/G calls appear only at alignment ends, in divergent N- or C-terminal segments, or at group-I intron junctions. Each is supported by one reference only, and none survives the default setting. Examples: N-termini of cox2, cox3 and atp6, the C-termini of cob and cox3, a cox1 C-terminal stretch in genome_2, and the junctions of cob intron 1 and the cox1 introns in genome_3. I did not count any of them.
* **No −1 frameshifts** were found in any core gene. The only G call, `refA_cob_v2` in genome_2 at 23632–23633, lies in the cob intron-ORF part of that fusion reference, where identity is about 28%. It is not credible.
* **Intron ORFs.** The intron-ORF parts of refA_cob_v2 and v3 were mapped separately without splicing. They have no homolog in genome_1 or genome_2. In genome_3 they have only a weak, frameshift-free hit inside a cox1 intron (15734–16759, about 24% nt-level matches). So frameshifts were only scanned reliably in the 7 core genes.
* **Starts and ends.** Every gene starts at or within a few codons of an in-frame ATG or ATA, with no stop in between, and ends at TAA. No terminal frameshift was found.
* **Exact positions are uncertain.** Section 4 and 5 give the window of equally good positions for each site. Counts and directions do not depend on that choice. The residue "immediately after" does, so both versions are given where they differ.
* **What cannot be detected:** a frameshift exactly at a splice junction or inside an intron (yeast group I/II introns do not follow miniprot's GT–AG model, so junction placement is uncertain by a few nt), and compensating +1/−1 pairs a codon or two apart. The `-F 6` run would expose such pairs, and it showed none away from the terminal artifacts above.
* **Mechanism is not determined.** The data show recurring patterns, which argues against random sequencing errors. They cannot say whether the shifts are programmed translational frameshifts, consistent assembly or sequencing artifacts, or engineered insertions.

---

## 7. Exact commands

All work was done inside `<blind folder>`, using `/opt/homebrew/bin/miniprot` (version `0.18-r281`) and `python3 -I`. Helper scripts are in `work/scripts/`.

### 7a. Final consolidated pipeline (`work/reproduce.sh`; outputs in `work/run/`, which the tables above cite)

Run as `bash work/reproduce.sh`. The full file:

```bash
#!/bin/bash
# Full pipeline for VERDICT.md. Run from anywhere: bash work/reproduce.sh
set -euo pipefail
D=<blind folder>
cd "$D"
MP=/opt/homebrew/bin/miniprot
S=work/scripts
O=work/run
mkdir -p $O
glen() { case $1 in 1) echo 32558;; 2) echo 29214;; 3) echo 31230;; esac; }   # genome lengths (composition.txt)

# 0. inputs
$MP --version > $O/miniprot_version.txt
python3 -I $S/composition.py genome_1.fna genome_2.fna genome_3.fna > $O/composition.txt
python3 -I $S/prep_refs.py reference_proteins.faa > $O/refs.faa      # refA_cob x3 -> refA_cob_v1/v2/v3

for g in 1 2 3; do
  # 1. linear genomes, default settings, yeast mito code
  $MP -T 3 -t 1 -u genome_$g.fna $O/refs.faa > $O/g${g}_linear.paf 2> $O/g${g}_linear.log
  # 2. circular handling: append first 8 kb (genes crossing the sequence end become contiguous)
  python3 -I $S/circularize.py genome_$g.fna 8000 > $O/genome_${g}_circ.fna
  $MP -T 3 -t 1 -u $O/genome_${g}_circ.fna $O/refs.faa > $O/g${g}_circ.paf 2> $O/g${g}_circ.log
  $MP -T 3 -t 1 --gff --aln $O/genome_${g}_circ.fna $O/refs.faa > $O/g${g}_circ.gff_aln 2>/dev/null
  # 3. independent check of the end-crossing genes: genome rotated to start at its midpoint
  python3 -I $S/rotate.py genome_$g.fna > $O/genome_${g}_rot.fna
  $MP -T 3 -t 1 $O/genome_${g}_rot.fna $O/refs.faa > $O/g${g}_rot.paf 2>/dev/null
  # 4. sensitivity sweep of the frameshift penalty (default -F 23)
  for F in 6 12 15 35; do
    $MP -T 3 -t 1 -u -F $F $O/genome_${g}_circ.fna $O/refs.faa > $O/g${g}_circ_F$F.paf 2>/dev/null
  done
done

for g in 1 2 3; do
  L=$(glen $g)
  echo "### genome_$g" 
  # 5. frameshift events (cg ops F/G) from the circular default run
  python3 -I $S/fs_events.py $O/genome_${g}_circ.fna $O/refs.faa $O/g${g}_circ.paf $L
  # 6. placement windows (all equally good positions for each +1 event)
  python3 -I $S/fs_window.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf
  # 7. in-frame codon usage per genome, start/stop check
  python3 -I $S/codons.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf
  python3 -I $S/termini.py $O/genome_${g}_circ.fna $O/refs.faa $O/g${g}_circ.paf | sort -u -k2,2
  # 8. where GGC and CTG fall relative to the reading frame
  for m in GGC CTG; do python3 -I $S/motif_scan.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf $m; done
  python3 -I $S/motif_null.py $O/genome_${g}_circ.fna $L $O/refs.faa $O/g${g}_circ.paf AGGC,TGGC,GGC,CTGT,CTGC,CTG x
done > $O/analysis.txt

# 9. residue-level alignment windows around every frameshift
for g in 2 3; do python3 -I $S/aln_windows.py $O/g${g}_circ.gff_aln 40; done > $O/aln_windows.txt

# 10. codons of all three genomes at each site (reference-anchored)
SP="$O/g1_circ.paf:$O/genome_1_circ.fna:32558 $O/g2_circ.paf:$O/genome_2_circ.fna:29214 $O/g3_circ.paf:$O/genome_3_circ.fna:31230"
for x in "refA_cox1 70 86" "refA_cob_v1 130 150" "refA_cob_v1 206 220" "refA_cob_v1 74 88" "refA_cox2 218 232" "refA_cox3 116 130" "refA_cox3 244 258" "refB_atp6 184 198"; do
  echo "### $x"; python3 -I $S/compare_sites.py $O/refs.faa $x $SP
done > $O/compare_sites.txt

# 11. codons around each site (canonical placement), original-frame run-on to a stop, local GC
{
  echo "## genome_2"; python3 -I $S/site_codons.py genome_2.fna 24711 - 26496 - 26477 - 26263 - 26266 -
  python3 -I $S/stop_distance.py genome_2.fna 24711 - 26496 - 26477 - 26263 -
  python3 -I $S/local_context.py genome_2.fna 24711 - 26496 - 26477 - 26263 -
  echo "## genome_3"; python3 -I $S/site_codons.py genome_3.fna 121 + 27197 + 9815 + 10197 + 10200 + 19803 + 19800 +
  python3 -I $S/stop_distance.py genome_3.fna 121 + 27197 + 9815 + 10197 + 19803 +
  python3 -I $S/local_context.py genome_3.fna 121 + 27197 + 9815 + 10197 + 19803 +
} > $O/sites.txt

# 12. refA_atp6 did not map to genome_3: retry with more sensitive settings
awk '/^>refA_atp6$/{p=1;print;next} /^>/{p=0} p' $O/refs.faa > $O/refA_atp6.faa
for opt in "" "-M 0" "-k 5" "-M 0 -k 5 -n 1" "-S"; do
  echo "opts=[$opt]"; $MP -T 3 -t 1 -u $opt $O/genome_3_circ.fna $O/refA_atp6.faa 2>/dev/null | cut -f1-12
done > $O/refA_atp6_retry.txt

# 13. cob intron-ORF parts of refA_cob_v2/v3 mapped on their own (no splicing)
python3 -I - "$O/refs.faa" > $O/intron_orfs.faa <<'PY'
import sys
seqs = {}; n = None
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith('>'): n = l[1:]; seqs[n] = ''
    else: seqs[n] += l
print('>refA_cob_v2_intronORF_172-517\n' + seqs['refA_cob_v2'][171:])
print('>refA_cob_v3_intronORF_143-423\n' + seqs['refA_cob_v3'][142:])
PY
for g in 1 2 3; do $MP -T 3 -t 1 -u -S $O/genome_${g}_circ.fna $O/intron_orfs.faa 2>/dev/null | cut -f1-12,14-19; done > $O/intron_orfs.paf
echo done
```

What each helper script does (all in `work/scripts/`; each has a usage line at the top):

| Script | Purpose |
|---|---|
| `prep_refs.py` | Rename the duplicate `refA_cob` headers to `_v1/_v2/_v3`. |
| `circularize.py` | Append the first N bases of a genome to its end. |
| `rotate.py` | Rotate a genome to start at its midpoint. |
| `composition.py` | Genome length and base composition. |
| `fs_events.py` | Walk the `cg:Z:` CIGAR (M/I/D/F/G/N/U/V) and report each F/G event: forward coordinate, extra bases, reference residue after, codon context. |
| `fs_window.py` | For each F event, try every placement of the extra base(s) in an ungapped window bounded by neighbouring events. Report all placements with maximal identity to that reference (code 3 translation). |
| `aln_windows.py` | Cut miniprot `--aln` blocks around frameshift marks (`!`). |
| `compare_sites.py` | Codons of all three genomes aligned to the same reference residues. |
| `codons.py` | In-frame codon counts for the 7 genes (M/D codons from the best alignment per gene), plus gene coordinates. |
| `termini.py` | Start codon, last codon, and distance to the next TAA/TAG. |
| `motif_scan.py` | Every occurrence of a motif in the genes, with its codon position in the aligned frame. |
| `motif_null.py` | Coding-sequence composition and per-position motif frequency. |
| `site_codons.py` | Original-frame and shifted-frame codons around a site. |
| `stop_distance.py` | How far the unshifted frame runs before a stop. |
| `local_context.py` | 81-nt coding-strand context and GC fraction around a site. |

### 7b. Exploratory commands run before the consolidated pipeline (outputs in `work/`)

These used the same tools and parameters. They are listed in the order run. In the loops, `g` = 1, 2, 3 and `$L` = that genome's length (32558, 29214, 31230). Ad-hoc inline Python one-liners from this stage (composition, renaming, rotating, circularizing, context printing, residue look-ups) did the same jobs as the scripts above and were later saved as those scripts.

```bash
cd <blind folder>
# linear runs
for g in 1 2 3; do /opt/homebrew/bin/miniprot -T 3 -t 1 genome_$g.fna work/refs_renamed.faa > work/g${g}.paf 2> work/g${g}.paf.log; done
for g in 1 2 3; do /opt/homebrew/bin/miniprot -T 3 -t 1 --gff --aln genome_$g.fna work/refs_renamed.faa > work/g${g}.gffaln 2>/dev/null; done
python3 -I work/scripts/fs_events.py genome_$g.fna work/refs_renamed.faa work/g$g.paf          # g = 1,2,3
# rotated genomes (midpoint)
/opt/homebrew/bin/miniprot -T 3 -t 1 work/genome_${g}_rot.fna work/refs_renamed.faa > work/g${g}_rot.paf
# circularized genomes (+8 kb), default and -F 12
/opt/homebrew/bin/miniprot -T 3 -t 1 -u work/genome_${g}_circ.fna work/refs_renamed.faa > work/g${g}_circ.paf
/opt/homebrew/bin/miniprot -T 3 -t 1 -u -F 12 work/genome_${g}_circ.fna work/refs_renamed.faa > work/g${g}_circ_F12.paf
/opt/homebrew/bin/miniprot -T 3 -t 1 --aln work/genome_${g}_circ.fna work/refs_renamed.faa > work/g${g}_circ.aln
python3 -I work/scripts/aln_windows.py work/g2_circ.aln 40 ; python3 -I work/scripts/aln_windows.py work/g3_circ.aln 40
python3 -I work/scripts/fs_window.py work/genome_${g}_circ.fna $L work/refs_renamed.faa work/g${g}_circ.paf
python3 -I work/scripts/codons.py work/genome_${g}_circ.fna $L work/refs_renamed.faa work/g${g}_circ.paf
python3 -I work/scripts/compare_sites.py work/refs_renamed.faa <REF> <START> <END> work/g1_circ.paf:work/genome_1_circ.fna:32558 work/g2_circ.paf:work/genome_2_circ.fna:29214 work/g3_circ.paf:work/genome_3_circ.fna:31230
# penalty sweep
for F in 6 15 35; do /opt/homebrew/bin/miniprot -T 3 -t 1 -F $F work/genome_${g}_circ.fna work/refs_renamed.faa > work/g${g}_circ_F$F.paf; done
# intron-ORF parts alone, no splicing
/opt/homebrew/bin/miniprot -T 3 -t 1 -u -S work/genome_${g}_circ.fna work/intron_orfs.faa
# start/stop checks and motif scans
python3 -I work/scripts/termini.py work/genome_${g}_circ.fna work/refs_renamed.faa work/g${g}_circ.paf
python3 -I work/scripts/motif_scan.py work/genome_${g}_circ.fna $L work/refs_renamed.faa work/g${g}_circ.paf GGC   # and CTG
python3 -I work/scripts/motif_null.py work/genome_${g}_circ.fna $L work/refs_renamed.faa work/g${g}_circ.paf <MOTIFS> x
# refA_atp6 retries on genome_3 (whole circ genome and the extracted 19001-20300 region)
/opt/homebrew/bin/miniprot -T 3 -t 1 -u [ "" | -M 0 | -k 5 | "-M 0 -k 5 -n 1" | -S ] work/genome_3_circ.fna work/refA_atp6.faa
/opt/homebrew/bin/miniprot -T 3 -t 1 -u [ "" | "-k 5 -n 1 -m 0" | "-k 4 -n 1 -l 4" ] work/g3_atp6_region.fna work/refA_atp6.faa
python3 -I work/scripts/site_codons.py genome_2.fna 24711 - 26496 - 26477 - 26263 - 26266 -
python3 -I work/scripts/site_codons.py genome_3.fna 121 + 27197 + 9815 + 10197 + 10200 + 19803 + 19800 +
```

Two bugs were fixed during the exploratory stage, and all tables above use the fixed versions:
* An early translation table in `fs_window.py` and `fs_events.py` read CTN as Leu. It now uses code 3, where CTN is Thr.
* An ad-hoc motif listing counted `I` operations as consuming genome bases. This was replaced by `motif_scan.py`.
