library("seqinr")
library("kableExtra")
library("R.utils")

#QUESTION 1
#downloading E.coli coding DNA sequence

URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile = "ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")


#downloading Saprospirales bacterium (GCA_003448025) coding DNA sequence

URL1="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_58_collection/saprospirales_bacterium_gca_003448025/cds/Saprospirales_bacterium_gca_003448025.ASM344802v1.cds.all.fa.gz"
download.file(URL1,destfile="Sbac_cds.fa.gz")
gunzip("Sbac_cds.fa.gz")



#listing the length of cds for E.coli and S.bacterium
cds <- seqinr::read.fasta("ecoli_cds.fa")
str(head(cds))
cds2 <- seqinr:: read.fasta("Sbac_cds.fa")

CS_E.coli <- length(cds)   #coding seq for E.coli calculating the no. of coding sequences
CS_S <- length(cds2) # coding seq for S.bacterium calculating the no. of coding sequences

#QUESTION 2
head(summary(cds))
str(summary(cds))

head( summary(cds)[,1])
len <- sapply(X=cds,FUN=length) #calculating the
median(len)
mean(len)

head(summary(cds2))
str(summary(cds2))
head( summary(cds2)[,1])
len2 <- sapply(X=cds,FUN=length)
median(len2)
mean(len2)













