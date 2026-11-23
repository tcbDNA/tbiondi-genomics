# metazoan-uaa-tyr-check

## Resource rule (set by the user, 2026-11-23; replaces the 20 GB / 1 hour approval gates)
- Downloads and long-running jobs no longer need the user's approval.
- Keep at least 30 GB of free disk at all times. Before any large download or job, estimate its size; if it would leave less than 30 GB free, stop and ask.
- Log every download (URL/accession, size, checksum) in the download ledger (`data/manifests/download_ledger.tsv`), and report cumulative totals in each gate report.
- Run long jobs in the background with progress logged to a file, so `/status` can see them.
- Keep the other projects' saved state intact. Do not delete their files to make room.

## Pre-registration
- No trnY, UAA or mtRF1a caller runs on a non-control genome until `PREREG.md` is committed. The user reviews the revised `PREREG.md` before it is committed.
