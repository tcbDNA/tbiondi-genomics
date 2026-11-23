#!/usr/bin/env bash
# Fetch every public input by accession into data/raw/ (about 1 GB in total, plus 0.3 GB of reads if wanted).
#
#   bash scripts/setup_data.sh            # genomes, taxonomy, fungal TPA records, held-out record
#   bash scripts/setup_data.sh --reads    # also the Illumina run behind the Grigorovia jiainica mitogenome
#   bash scripts/setup_data.sh --genus    # also the Grigorovia follow-up inputs (1.6 GB of reads, Rfam models; the
#                                         # nuclear assemblies, UniProt panel and BUSCO data are fetched by the scripts)
#
# Notes:
# - The scan used RefSeq mitochondrion release 238 (2026). NCBI's FTP serves only the current release. To rebuild
#   release 238 exactly, this script fetches the 18,220 RefSeq records listed in protocol/accessions.tsv by
#   accession.version (scripts/fetch_by_accession.py, NCBI E-utilities, no API key, a few hours at the polite rate).
# - NCBI Taxonomy changes over time. The rows used are in protocol/taxonomy_extract/ (see step 2); the codes and pools are
#   also recorded per genome in protocol/accessions.tsv (genetic_code) and protocol/taxonomy_assignments.tsv.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw/refseq_mito data/raw/fungal_tpa data/raw/taxdump data/raw/OR703946

# 1. RefSeq mitochondrion records used by the scan (release 238), by accession.version.
uv run python scripts/fetch_by_accession.py \
    --accessions <(awk -F'\t' 'NR>1 && $2 ~ /^refseq/ {print $1}' protocol/accessions.tsv) \
    --gbff data/raw/refseq_mito/mitochondrion.1.genomic.gbff.gz \
    --fasta data/raw/refseq_mito/mitochondrion.1.1.genomic.fna.gz

# 2. NCBI taxonomy. The analysis reads protocol/taxonomy_extract/: the rows it uses, cut from the dump downloaded on
#    2026-11-21 (taxdump.tar.gz md5 6d18d6c9d60bbc20759869f6b3592eff), which reproduces the published values. That
#    dump is not archived by NCBI. For an independent check, the closest dated archive is fetched and verified here;
#    it differs for 4 of 20,915 genomes (1 taxon added after 2026-11-01; 3 family/order changes), so
#    'prep.py --taxonomy data/raw/taxdump' stops at those 4 unless --use-archived-taxonomy is added.
curl -fL -o data/raw/taxdump/taxdmp_2026-11-01.zip https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2026-11-01.zip
python3 -c 'import hashlib,sys; h=hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest(); sys.exit(0 if h=="42a27451954fd10ddcb50f0549ab7869" else "taxdmp_2026-11-01.zip: md5 "+h+", expected 42a27451954fd10ddcb50f0549ab7869")' data/raw/taxdump/taxdmp_2026-11-01.zip
unzip -o -q data/raw/taxdump/taxdmp_2026-11-01.zip names.dmp nodes.dmp merged.dmp -d data/raw/taxdump

# 3. The 2,695 fungal mitogenomes of Sarhan et al. 2025 (GenBank TPA BK072095-BK074789) and their accession table
#    (figshare, doi:10.6084/m9.figshare.28750034; save it as data/raw/fungal_tpa/All_fungal_mitochondrial_genomes_accessions.tsv).
uv run python code/audit/fetch_fungal_tpa.py

# 4. Held-out relative: Grigorovia yakushimaensis mitochondrial sequence OR703946.1 (GenBank, UNVERIFIED).
curl -fL -o data/raw/OR703946/OR703946.1.gb \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=OR703946.1&rettype=gbwithparts&retmode=text"

# 5. Optional: reads behind the Grigorovia jiainica mitogenome (SRA SRR16974510, Illumina HiSeq 2500 paired, 282 MB).
if [[ "${1:-}" == "--reads" ]]; then
    mkdir -p data/reads/SRR16974510
    for r in 1 2; do
        curl -fL -C - --retry 5 -o data/reads/SRR16974510/SRR16974510_$r.fastq.gz \
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/010/SRR16974510/SRR16974510_$r.fastq.gz"
    done
fi
# 6. Optional: inputs of the Grigorovia follow-up tests (candidate/followup/), from NCBI's AWS mirror of SRA.
if [[ "${1:-}" == "--genus" || "${2:-}" == "--genus" ]]; then
    mkdir -p data/rfam
    for m in RF00028 RF00029; do curl -fL -o data/rfam/$m.cm https://rfam.org/family/$m/cm; done
    for r in SRR16974597 SRR16974369 ERR9706995; do
        d=data/reads/genus/$r; mkdir -p $d
        curl -fL --retry 5 -o $d/$r.sra "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r"
        fasterq-dump --split-files --threads 6 -O $d $d/$r.sra
    done
fi
echo "done"
