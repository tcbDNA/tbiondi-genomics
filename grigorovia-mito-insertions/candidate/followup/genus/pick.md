# Mitochondrial contigs picked from the SPAdes assemblies

Rule: a contig is mitochondrial if it has a DIAMOND blastx hit (code 3, e <= 1e-10) to an S. cerevisiae mitochondrial protein, or a minimap2 (asm20) alignment of >= 300 nt to the congener mitogenome OR703946, or GC < 25% (yeast mtDNA is AT-rich; nuclear contigs here are 30-40%), and k-mer coverage >= 3x the assembly's length-weighted median (mtDNA is multi-copy).

## G. humatica

- assembly: 6556 contigs, 13,388,792 nt; length-weighted median k-mer coverage 18.4
- mitochondrial contigs kept: 5, total 62,112 nt

  - NODE_74_length_53413_cov_109.893787: 53,413 nt, coverage 110, genes atp6,atp8,atp9,cob,cox1,cox2,cox3, aligned to OR703946 3,683 nt, GC 0.182
  - NODE_202_length_7239_cov_79.569673: 7,239 nt, coverage 80, genes VAR1, aligned to OR703946 2,462 nt, GC 0.163
  - NODE_514_length_563_cov_75.705761: 563 nt, coverage 76, genes -, aligned to OR703946 0 nt, GC 0.249
  - NODE_662_length_472_cov_153.959494: 472 nt, coverage 154, genes -, aligned to OR703946 0 nt, GC 0.246
  - NODE_959_length_425_cov_79.399425: 425 nt, coverage 79, genes -, aligned to OR703946 0 nt, GC 0.212

- matching but low coverage, >= 500 nt (possible nuclear copies, not used): 31
  - NODE_90_length_44593_cov_17.621439: 44,593 nt, coverage 17.6, genes -
  - NODE_117_length_33800_cov_16.010349: 33,800 nt, coverage 16.0, genes -
  - NODE_182_length_11583_cov_17.894664: 11,583 nt, coverage 17.9, genes -
  - NODE_245_length_3401_cov_33.549940: 3,401 nt, coverage 33.5, genes -
  - NODE_350_length_1048_cov_26.020597: 1,048 nt, coverage 26.0, genes -
  - NODE_355_length_984_cov_1.100331: 984 nt, coverage 1.1, genes -
  - NODE_361_length_948_cov_1.939150: 948 nt, coverage 1.9, genes cox3
  - NODE_369_length_895_cov_1.037897: 895 nt, coverage 1.0, genes cob
  - NODE_375_length_879_cov_2.375312: 879 nt, coverage 2.4, genes cox2
  - NODE_378_length_877_cov_1.712500: 877 nt, coverage 1.7, genes -

## G. transvaalensis

- assembly: 15060 contigs, 14,056,744 nt; length-weighted median k-mer coverage 7.6
- mitochondrial contigs kept: 27, total 33,261 nt

  - NODE_741_length_5483_cov_78.715825: 5,483 nt, coverage 79, genes cox3, aligned to OR703946 265 nt, GC 0.251
  - NODE_749_length_5432_cov_230.272267: 5,432 nt, coverage 230, genes VAR1, aligned to OR703946 58 nt, GC 0.292
  - NODE_899_length_4188_cov_112.521154: 4,188 nt, coverage 113, genes cob, aligned to OR703946 345 nt, GC 0.275
  - NODE_1265_length_2119_cov_101.535149: 2,119 nt, coverage 102, genes cox1, aligned to OR703946 0 nt, GC 0.269
  - NODE_1298_length_2008_cov_64.216344: 2,008 nt, coverage 64, genes cox1, aligned to OR703946 0 nt, GC 0.196
  - NODE_1463_length_1388_cov_50.504267: 1,388 nt, coverage 51, genes -, aligned to OR703946 467 nt, GC 0.218
  - NODE_1467_length_1377_cov_43.723787: 1,377 nt, coverage 44, genes -, aligned to OR703946 0 nt, GC 0.203
  - NODE_1496_length_1291_cov_61.166946: 1,291 nt, coverage 61, genes atp6, aligned to OR703946 235 nt, GC 0.232
  - NODE_1533_length_1156_cov_83.543992: 1,156 nt, coverage 84, genes cox2, aligned to OR703946 288 nt, GC 0.226
  - NODE_1552_length_1092_cov_128.586103: 1,092 nt, coverage 129, genes cox1, aligned to OR703946 131 nt, GC 0.281
  - NODE_1559_length_1075_cov_128.731557: 1,075 nt, coverage 129, genes atp9, aligned to OR703946 89 nt, GC 0.284
  - NODE_1591_length_1004_cov_32.670718: 1,004 nt, coverage 33, genes -, aligned to OR703946 0 nt, GC 0.214
  - NODE_1651_length_862_cov_27.431193: 862 nt, coverage 27, genes -, aligned to OR703946 48 nt, GC 0.188
  - NODE_1689_length_759_cov_37.775758: 759 nt, coverage 38, genes cox1, aligned to OR703946 45 nt, GC 0.241
  - NODE_1737_length_636_cov_67.638734: 636 nt, coverage 68, genes -, aligned to OR703946 138 nt, GC 0.233
  - NODE_1768_length_596_cov_59.173038: 596 nt, coverage 59, genes -, aligned to OR703946 1,509 nt, GC 0.273
  - NODE_1787_length_563_cov_77.250000: 563 nt, coverage 77, genes -, aligned to OR703946 0 nt, GC 0.243
  - NODE_2036_length_375_cov_26.202899: 375 nt, coverage 26, genes atp8, aligned to OR703946 0 nt, GC 0.232
  - NODE_2182_length_313_cov_40.504673: 313 nt, coverage 41, genes -, aligned to OR703946 0 nt, GC 0.211
  - NODE_2410_length_264_cov_26.896970: 264 nt, coverage 27, genes cox2, aligned to OR703946 0 nt, GC 0.216
  - NODE_2425_length_262_cov_42.251534: 262 nt, coverage 42, genes atp6, aligned to OR703946 54 nt, GC 0.206
  - NODE_2909_length_195_cov_25.135417: 195 nt, coverage 25, genes cox1, aligned to OR703946 0 nt, GC 0.451
  - NODE_4395_length_191_cov_23.413043: 191 nt, coverage 23, genes cox1, aligned to OR703946 0 nt, GC 0.414
  - NODE_4396_length_191_cov_23.358696: 191 nt, coverage 23, genes cob, aligned to OR703946 0 nt, GC 0.450
  - NODE_6366_length_171_cov_22.833333: 171 nt, coverage 23, genes cox1, aligned to OR703946 0 nt, GC 0.392
  - NODE_6984_length_152_cov_25.018868: 152 nt, coverage 25, genes cob, aligned to OR703946 0 nt, GC 0.362
  - NODE_9082_length_118_cov_23.842105: 118 nt, coverage 24, genes cob, aligned to OR703946 0 nt, GC 0.263

- matching but low coverage, >= 500 nt (possible nuclear copies, not used): 5
  - NODE_166_length_16922_cov_9.356714: 16,922 nt, coverage 9.4, genes -
  - NODE_500_length_8641_cov_12.322524: 8,641 nt, coverage 12.3, genes -
  - NODE_1755_length_611_cov_1.783203: 611 nt, coverage 1.8, genes -
  - NODE_1807_length_536_cov_1.736842: 536 nt, coverage 1.7, genes -
  - NODE_1835_length_508_cov_1.755501: 508 nt, coverage 1.8, genes -

## G. transvaalensis_NRRL_Y-17245

- assembly: 1068 contigs, 12,204,016 nt; length-weighted median k-mer coverage 3.8
- mitochondrial contigs kept: 0, total 0 nt


- matching but low coverage, >= 500 nt (possible nuclear copies, not used): 4
  - NODE_26_length_89661_cov_3.823216: 89,661 nt, coverage 3.8, genes cox2
  - NODE_87_length_45739_cov_3.736302: 45,739 nt, coverage 3.7, genes -
  - NODE_172_length_20835_cov_4.068070: 20,835 nt, coverage 4.1, genes -
  - NODE_429_length_1024_cov_2.693770: 1,024 nt, coverage 2.7, genes -

