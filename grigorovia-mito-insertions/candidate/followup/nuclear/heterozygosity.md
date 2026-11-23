# G. jiainica nuclear heterozygosity and heterozygous variants in the panel genes

Genome-wide (25 largest contigs, 3,654,777 nt; sites with 20-120 reads at MAPQ >= 20, base quality >= 20): 2,987,281 sites.
- heterozygous-like (second base in 25-75% of reads): 35,687 = 11.95 per kb
- reads disagree with the assembly (> 90% another base): 20 = 0.007 per kb
- alternative-base fraction, sites with any: {0.0: 209807, 0.1: 11266, 0.2: 3917, 0.3: 8308, 0.4: 12132, 0.5: 9959, 0.6: 3768, 0.7: 960, 0.8: 185, 0.9: 47, 1.0: 4}

## Panel genes (locus from the best DIAMOND hit of the congener protein); pileup view

| gene | locus | heterozygous SNVs | indels in >= 25% of reads (not multiple of 3) | in-frame indels | stop-gains (single-base view) |
|---|---|---|---|---|---|
| ATP22 | JAJLYE010000001.1:178497-181055 | 39 | - | 179306+TCA(5/19) | - |
| CBP1 | JAJLYE010000102.1:10573-11652 | 0 | - | - | - |
| CBP3 | JAJLYE010000011.1:18036-18827 | 9 | - | - | - |
| CBP6 | JAJLYE010000229.1:48002-48460 | 0 | - | - | - |
| CBS2 | JAJLYE010000024.1:4204-5331 | 3 | - | - | - |
| COR1 | JAJLYE010000008.1:115588-117017 | 41 | 115929-TC(17/26) | - | - |
| COX10 | JAJLYE010000023.1:104210-105244 | 3 | - | - | - |
| COX12 | JAJLYE010000002.1:231819-232058 | 0 | - | - | - |
| COX13 | JAJLYE010000165.1:3021-3404 | 3 | - | - | - |
| COX15 | JAJLYE010000132.1:2350-3723 | 12 | - | - | - |
| COX4 | JAJLYE010000128.1:1704-2135 | 2 | - | - | - |
| COX5A | JAJLYE010000231.1:10327-10809 | 4 | - | - | - |
| COX5B | JAJLYE010000231.1:10327-10770 | 3 | - | - | - |
| COX6 | JAJLYE010000251.1:19135-19494 | 6 | - | - | - |
| COX7 | JAJLYE010000063.1:17583-17765 | 3 | - | - | - |
| COX8 | JAJLYE010000035.1:19592-19816 | 0 | - | - | - |
| COX9 | JAJLYE010000082.1:25266-25421 | 3 | - | - | - |
| CYC1 | JAJLYE010000282.1:11489-11794 | 0 | - | - | - |
| CYC3 | JAJLYE010000072.1:36412-37242 | 9 | - | 36536-TAATAA(10/23) | - |
| CYC7 | JAJLYE010000282.1:11456-11794 | 0 | - | - | - |
| CYT1 | JAJLYE010000239.1:10664-11485 | 11 | - | - | - |
| MEF1 | JAJLYE010000030.1:17443-19686 | 28 | - | - | - |
| MRP4 | JAJLYE010000001.1:244650-245693 | 8 | - | - | - |
| MRPL10 | JAJLYE010000071.1:36865-37710 | 2 | - | - | - |
| MRPL15 | JAJLYE010000245.1:35515-36114 | 0 | - | - | - |
| MRPL16 | JAJLYE010000041.1:47790-48380 | 6 | - | - | - |
| MRPL3 | JAJLYE010000062.1:47871-48851 | 33 | - | - | - |
| MRPL4 | JAJLYE010000267.1:26367-26852 | 0 | - | - | - |
| MRPL7 | JAJLYE010000024.1:88871-89671 | 17 | - | - | - |
| MRPS17 | JAJLYE010000042.1:2098-2850 | 5 | - | - | - |
| MRPS5 | JAJLYE010000025.1:85863-86660 | 7 | - | - | - |
| MRPS9 | JAJLYE010000003.1:175033-175812 | 13 | - | - | - |
| MSS51 | JAJLYE010000210.1:42831-44156 | 0 | - | - | - |
| PET111 | JAJLYE010000063.1:17979-20309 | 71 | 20114-ATTATTA(7/25) | - | - |
| PET309 | JAJLYE010000030.1:20767-23227 | 40 | 22944-ATTGAAGAAGGAAAAA(7/25), 22956+TT(10/25) | 22963+ATT(12/22) | 22025C>G(27/54) |
| PET494 | JAJLYE010000929.1:27776-29473 | 23 | - | 28639-ACGTTC(18/45), 28951+ATC(13/38) | - |
| QCR10 | JAJLYE010000002.1:52753-52977 | 1 | - | - | - |
| QCR2 | JAJLYE010000060.1:45570-46259 | 1 | - | - | - |
| QCR6 | JAJLYE010000320.1:5727-5939 | 0 | - | - | - |
| QCR7 | JAJLYE010000084.1:20241-20621 | 4 | - | - | - |
| QCR8 | JAJLYE010000004.1:70109-70384 | 2 | - | - | - |
| QCR9 | JAJLYE010000070.1:27873-28073 | 0 | - | - | - |
| RIP1 | JAJLYE010000014.1:67141-67767 | 0 | - | - | - |
| TUF1 | JAJLYE010000269.1:17593-18837 | 12 | - | - | - |

## Flagged spots: reads grouped by the exact sequence they carry (alignment-based; code/explore/read_haplotypes.py)

| gene | window | read versions (reads) | length vs assembly | translation in the gene frame | reads matching the assembly |
|---|---|---|---|---|---|
| ATP22 | JAJLYE010000001.1:179294-179334 | 7 / 6 (of 14) | +3 / +0 | SSSSSSSSSSAQNN (7); SSSSSSSSSAQNN (6) | 6 |
| COR1 | JAJLYE010000008.1:115917-115957 | 18 / 5 (of 25) | -2 / +19 | IDKFDKIDSLFK (18); IDNANNKINKFDKIESLFK (5) | 0 |
| CYC3 | JAJLYE010000072.1:36524-36564 | 11 / 6 (of 23) | -6 / +0 | KDNSNNNNLMI (11); KDNSNNNNNNLMI (6) | 6 |
| PET111 | JAJLYE010000063.1:20102-20142 | 15 / 8 (of 26) | +0 / -6 | VNNNNNTNNNNNN (15); VNNNNNDINNN (8) | 15 |
| PET309 | JAJLYE010000030.1:22013-22053 | 18 / 10 (of 36) | +0 / +0 | QGNRNECKLILNE (18); QGNSNECKLILNE (10) | 10 |
| PET309 | JAJLYE010000030.1:22932-22972 | 12 / 9 (of 23) | -16 / +5 | DEEMMEEG (12); DEEMIEEGLKIIEEG (9) | 0 |
| PET494 | JAJLYE010000929.1:28627-28667 | 20 / 14 (of 38) | +0 / -6 | RRQEDINERLKIS (20); RRQEDINLKIS (14) | 20 |
| PET494 | JAJLYE010000929.1:28939-28979 | 16 / 12 (of 29) | +3 / +0 | RKANNNDDDDTNNV (16); RKANNNDDDTNNV (12) | 12 |
