library("seqinr")
library("kableExtra")
library("R.utils")

#downloading E.coli coding DNA sequence

URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile = "ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")

#downloading Saprospirales bacterium (GCA_003448025) coding DNA sequence
library("seqinr")
library("kableExtra")
library("R.utils")
URL1="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_58_collection/saprospirales_bacterium_gca_003448025/cds/"
download.file(URL1,destfile="Sbac_cds.fa.gz")
gunzip("Sbac_cds.fa.gz")



#listing the length of cds for E.coli and S.bacterium
cds <- seqinr::read.fasta("ecoli_cds.fa")
str(head(cds))
seqinr:: read.fasta("Sbac_cds.fa")



