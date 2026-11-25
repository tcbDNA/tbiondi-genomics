"""Download ledger. RECONSTRUCTED on 2026-11-23 from files on disk and session logs: no ledger was kept during Gates 0-2.
Rows with sha256=NA were deleted before hashing or not retained. API responses and aborted transfers were never metered."""
import os,hashlib,csv,glob
def sha(f):
    h=hashlib.sha256()
    with open(f,"rb") as fh:
        for b in iter(lambda: fh.read(1<<22),b""): h.update(b)
    return h.hexdigest()
R=[]
def add(gate,src,f,how="on_disk"): R.append([gate,src,f,os.path.getsize(f),sha(f),how])
def gone(gate,src,f,size,note): R.append([gate,src,f,size,"NA",note])
SP="data/sources/fulltext/"
for f,u in (("lin2026.pdf","https://www.biorxiv.org/content/10.64898/2026.08.14.744788v1.full.pdf (web-fetch tool)"),("zhao2026.pdf","user browser download (10.1186/s12915-026-02767-5)"),("zhao2026_MOESM1.docx","user browser download"),("zhao2026_MOESM4.docx","user browser download"),("zhao2026_MOESM9.pdf","static-content.springer.com MOESM9"),("PMC6319601.xml","Europe PMC"),("PMC2761399.xml","Europe PMC"),("PMC3121625.xml","Europe PMC"),("matsushima2026.xml","arphahub 192470 xml")):
    if os.path.exists(SP+f): add("0/2",u,SP+f)
add("0","http://eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz","data/pfam/Pfam-A_enone.tar.gz")
add("1","ftp.ncbi refseq/release/mitochondrion r238","data/ncbi/refseq_mito/mitochondrion.1.genomic.gbff.gz"); add("1","ftp.ncbi pub/taxonomy/new_taxdump","data/ncbi/taxonomy/new_taxdump.tar.gz")
add("1","data.orthodb.org odb13v0_species","data/orthodb/odb13v0_species.tab.gz"); add("1","ENA portal read_run RNA-Seq taxids","data/ena/invert_rnaseq.tsv")
for f in sorted(glob.glob("data/ncbi/organelle_fasta/*.fna")): add("1","NCBI efetch nuccore (organelle FASTA)",f)
FG={"mitogenomes.zip":"68295238","Hsq.chr.pep":"68294806","Pmu.chr.pep":"68294893","Ach.chr.pep":"68294896","Pvi.chr.pep":"68294980"}
for f,i in FG.items(): add("2",f"https://ndownloader.figshare.com/files/{i} (10.6084/m9.figshare.33250323, CC BY 4.0)","data/figshare/"+f)
add("2","https://zenodo.org/records/4284483 refseq89m (CC BY 4.0)","data/mitos_ref/refseq89m.tar.bz2")
AC={"AnoChi":"GCA_062075034.1","AnoEph":"GCA_990049579.1","EniAen":"GCA_056601665.1","HetSqu":"GCA_964263325.1","PodMac":"GCA_032359645.1","RadSim":"GCA_013357305.1"}
for n,a in AC.items(): add("2","NCBI datasets "+a,f"data/nuclear/{n}.zip")
for f in sorted(glob.glob("data/mtrf/bg/*.faa")): add("2","NCBI datasets protein "+os.path.basename(f)[:-4]+" (unpacked size; zip not kept)",f)
for f in sorted(glob.glob("env/pkgs/*.conda")+glob.glob("env/pkgs/*.tar.bz2")): add("0/2","conda-forge/bioconda",f)
add("2","bioconda mitos 2.1.10","env/tools/mitos_pkg/mitos-2.1.10-pyhdfd78af_0.conda")
S="https://sra-pub-run-odp.s3.amazonaws.com/sra/"
gone("2",S+"SRR11101524","SRR11101524.sra",1526990622,"deleted after search; size matched server; not hashed, not vdb-validated")
gone("2",S+"SRR11101525","SRR11101525.sra (1st)",1402354936,"CORRUPT (52.7 MB oversize, two writers); result discarded; deleted")
gone("2",S+"SRR11101525","SRR11101525.sra (2nd)",1349688568,"deleted after search; size matched server; vdb-validate rc 0")
gone("2",S+"ERR15194044","ERR15194044.sra",2336610385,"deleted after search; size matched server; not vdb-validated")
gone("2","ENA SRR1695471.fastq.gz","SRR1695471.fastq.gz",9975193,"deleted after search; size matched ENA fastq_bytes")
w=csv.writer(open("data/manifests/download_ledger.tsv","w"),delimiter="\t"); w.writerow("gate source_url_or_accession file bytes sha256 status".split()); w.writerows(R)
import collections; t=collections.Counter()
for r in R: t[r[0]]+=r[3]
print(len(R),"rows; total GB",round(sum(r[3] for r in R)/1e9,2),{k:round(v/1e9,2) for k,v in t.items()},"unhashed rows",sum(r[4]=="NA" for r in R))
