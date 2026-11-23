# Pass 1: S. cerevisiae queries (miniprot 0.18, standard code). Call = default run, or the sensitive run (*) where the default found nothing

Rules (PREREG): present >= 70% coverage, no frameshift or stop; broken = frameshift or stop; absent < 30%; else partial.

| gene | G. jiainica | G. yakushimaensis | G. humatica | G. transvaalensis | M. servazzii (control) | K. africana (control) |
|---|---|---|---|---|---|---|
| CYC1 | present 0.94 | present 0.95 | present 0.95 | present 0.95 | present 1.00 | present 0.94 |
| QCR6 | present 0.82 | present 0.86 | present 0.84 | present 0.86 | present 1.00 | partial 0.67 |
| QCR7 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 1.00 |
| COX5A | present 0.84 | present 1.00 | present 0.90 | present 0.87 | present 1.00 | present 1.00 |
| TUF1 | present 1.00 | present 0.98 | present 0.92 | present 0.92 | present 1.00 | present 1.00 |
| COX4 | present 0.93 | present 0.95 | present 1.00 | present 0.92 | present 1.00 | present 1.00 |
| COX8 | present 0.94 | present 0.94 | absent* | absent* | present 1.00 | present* 1.00 |
| CYT1 | present 0.86 | present 0.99 | present 0.88 | present 0.88 | present 0.95 | present 0.94 |
| COX9 | absent* | present* 0.97 | present* 0.90 | present* 0.90 | present 0.98 | present 1.00 |
| COR1 | broken 1.00 fs1 st0 | present 1.00 | present 0.94 | present 1.00 | present 1.00 | present 1.00 |
| QCR2 | present 0.93 | present 1.00 | present 0.96 | present 0.93 | present 1.00 | present 1.00 |
| RIP1 | present 0.94 | present 0.99 | present 0.97 | present 0.98 | present 1.00 | present 0.92 |
| QCR8 | present 1.00 | present 1.00 | present 0.99 | present 0.99 | present 1.00 | present 1.00 |
| COX7 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 1.00 |
| CBP3 | present 0.82 | present 0.80 | present 0.80 | present 0.81 | present 0.87 | present 0.79 |
| COX10 | present 0.75 | present 0.75 | present 0.80 | present 0.78 | present 0.97 | present 0.90 |
| QCR9 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | absent* | present 0.98 |
| MSS51 | present 0.96 | present 0.97 | present 0.94 | present 0.97 | present 0.99 | present 0.92 |
| COX13 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 0.87 | present 1.00 |
| MRPL4 | present 0.84 | broken 0.92 fs0 st1 | present 0.87 | present 0.83 | present 0.96 | present 0.97 |
| QCR10 | absent* | present 0.97 | present* 0.97 | present* 0.97 | present* 0.95 | present 1.00 |
| COX15 | present 0.83 | present 0.96 | present 0.85 | present 0.92 | present 0.83 | present 0.92 |
| COX12 | present 0.99 | present 0.96 | present 0.96 | present 0.96 | present 0.95 | present 0.98 |
| CYC7 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 1.00 | present 0.90 |
| COX5B | present 1.00 | present 1.00 | present 0.89 | present 0.85 | present 1.00 | present 1.00 |
| COX6 | present 0.92 | present 0.82 | present 0.84 | present 0.90 | present 1.00 | present 1.00 |
| CYC3 | present 0.96 | present 0.96 | present 0.95 | present 0.95 | present 1.00 | present 1.00 |
| CBP6 | present* 0.95 | present 0.95 | present 0.95 | present 1.00 | present 1.00 | present* 1.00 |
| PET122 | absent* | absent* | absent* | absent* | absent* | present* 0.82 |
| PET54 | absent* | absent* | absent* | absent* | absent* | present* 0.89 |
| CBS1 | absent* | absent* | absent* | absent* | absent* | absent* |
| AEP2 | absent* | absent* | absent* | absent* | absent* | absent* |
| MEF1 | present 0.96 | present 0.98 | present 0.96 | present 0.96 | present 0.97 | present 0.96 |
| PET309 | absent* | absent* | absent* | absent* | present* 0.94 | present 1.00 |
| MRPS5 | present 1.00 | present 0.94 | partial 0.62 | present 0.94 | present 1.00 | present 1.00 |
| MRPL3 | present 0.83 | present 0.89 | present 0.84 | present 0.89 | present 0.95 | present 0.84 |
| MRPL7 | present 0.92 | present 0.92 | present 0.92 | present 0.92 | present 0.98 | present 0.98 |
| MRPL10 | present 0.87 | present 0.87 | present 0.93 | present 0.84 | present 0.87 | present 0.88 |
| MRPL15 | absent* | present* 0.91 | present 0.81 | present* 0.78 | absent* | present 0.90 |
| MRPL16 | present 1.00 | present 0.98 | present 0.85 | present 0.85 | present 0.85 | present 0.85 |
| CBP2 | absent* | absent* | absent* | absent* | absent* | present 0.96 |
| CBP1 | absent* | absent* | absent* | absent* | present* 1.00 | present* 0.98 |
| PET494 | absent* | absent* | absent* | absent* | absent* | present 0.89 |
| PET111 | absent* | absent* | absent* | absent* | absent* | absent* |
| CBS2 | absent* | absent* | absent* | absent* | absent* | present 1.00 |
| MRP4 | present 0.85 | present 0.93 | present 0.87 | present 0.87 | present 1.00 | present 1.00 |
| MRPS9 | present 0.87 | present 0.87 | present 0.87 | present 0.87 | present 0.87 | present 1.00 |
| ATP22 | absent* | absent* | absent* | absent* | absent* | present 0.96 |
| MRPS17 | present 1.00 | present 1.00 | present 1.00 | present 0.99 | present 0.98 | present 0.98 |
