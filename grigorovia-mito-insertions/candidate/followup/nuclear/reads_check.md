# G. jiainica: assembly against read consensus

Reads mapped: 4467459 + 0 primary mapped (86.02% : N/A)

| gene | miniprot, assembly | miniprot, read consensus | DIAMOND, assembly (cov, id, frameshifts, stops) | DIAMOND, read consensus |
|---|---|---|---|---|
| COR1 | broken 1.00 fs1 st0 | present 1.00 fs0 st0 | 0.99, 61.6%, 1, 0 (JAJLYE010000008.1:115588-117017) | 0.99, 59.1%, 0, 0 (JAJLYE010000008.1:115573-117000) |
| COX9 | absent | absent | 0.91, 65.4%, 0, 0 (JAJLYE010000082.1:25266-25421) | 0.91, 65.4%, 0, 0 (JAJLYE010000082.1:25268-25423) |
| PET111 | absent | absent | 0.98, 31.3%, 0, 0 (JAJLYE010000063.1:20309-17979) | 0.98, 30.2%, 0, 0 (JAJLYE010000063.1:20297-17967) |
| PET122 | absent | absent | none | none |
| PET309 | broken 0.90 fs1 st0 | absent | 0.90, 36.5%, 1, 0 (JAJLYE010000030.1:20767-23227) | 0.90, 34.2%, 0, 0 (JAJLYE010000030.1:20733-23222) |

## Reads hitting the congener proteins (DIAMOND --sensitive, E <= 1e-5; both read files)

| gene | protein length | reads | mean identity |
|---|---|---|---|
| CBS2 | 423 | 225 | 52.5% |
| COR1 | 489 | 400 | 65.5% |
| COX12 | 81 | 89 | 87.5% |
| COX8 | 75 | 48 | 69.8% |
| COX9 | 57 | 56 | 68.8% |
| MSS51 | 429 | 741 | 70.2% |
| PET111 | 816 | 176 | 41.9% |
| PET122 | 260 | 296 | 51.8% |
| PET309 | 954 | 248 | 46.0% |
| QCR9 | 67 | 49 | 77.2% |

**PET122: missing from the assembly; assembled from the 203 read pairs that hit it (SPAdes 4.3.0):** contig 1657 nt, aligns to congener residues 2-254 of 260 at 46.8% identity, frameshifts 0, stops 0.
