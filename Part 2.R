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
boxplot(len,len2,main= "Length of Coding Sequences", at=c(1,2),
        names= c("E.Coli","S.bacterium"), ylab="Length")

#Question 3
#calculating dna composition frquency for E.coli
GC(cds[[1]])
count(cds[[1]],1)
count(cds[[1]],2)
count(cds[[1]],3)
summary(cds[1:3])
sum(sapply(cds[1:3],length))
length(unlist(cds[1:3]))
dna1 <- unlist(cds)
GC(dna1)
dna_composition <- count(dna1,1)
#Calculating dna composition frequency for S.bacterium

GC(cds2[[1]])
count(cds2[[1]],1)
count(cds2[[1]],2)
count(cds2[[1]],3)
summary(cds2[1:3])
sum(sapply(cds2[1:3],length))
length(unlist(cds2[1:3]))
dna2 <- unlist(cds2)
GC(dna2)
dna_composition2 <- count(dna2,1)


barplot(dna_composition,xlab="nucleotides",ylab="frequency", main="E.Coli CDS composition")

barplot(dna_composition2,xlab="nucleotides",ylab="frequency", main="S.bacterium CDS composition")



#Calculating protein sequence for E.coli
translate(cds[[1]])
prot <- lapply(cds,translate)
prot[[1]]
aa <-unique(prot[[2]])
aa <- aa[aa !="*"]

protein <- unlist(prot)
protein
protein_composition <- count(protein,alphabet=aa,wordsize=1)
protein_composition

protseq<- count(protein,alphabet=aa,wordsize=1)

#Calculating protein sequence for S.bacterium
translate(cds2[[1]])
prot2 <- lapply(cds2,translate)
prot2[[1]]
aa2 <-unique(prot[[2]])
aa2 <- aa2[aa2 !="*"]

protein2 <- unlist(prot2)
protein2
protein_composition2 <- count(protein2,alphabet=aa2,wordsize=1)
protein_composition2

protseq2<- count(protein2,alphabet=aa2,wordsize=1)



barplot(protseq,xlab="Protien",ylab="frequency", main="E.Coli CDS composition")
barplot(protseq2,xlab="Protien2",ylab="frequency", main="S.bacterium CDS composition")




#Question4
#Condon usage for E.coli

uco(cds[[2]])
codon_usage <- uco(cds[[2]],index = "rscu",as.data.frame =TRUE)

codon_usage
codon_usage%>%
  kbl()%>%
  kable_paper(full_width = F, html_font = "Arial")

plot(codon_usage)

#Condon usage for S.bacterium
uco(cds2[[2]])
codon_usage2 <- uco(cds2[[2]],index = "rscu",as.data.frame =TRUE)

codon_usage2
codon_usage2%>%
  kbl()%>%
  kable_paper(full_width = F, html_font = "Arial")
plot(codon_usage2)


#Question5
#Kmer profiling for E.coli
protein <- unlist(prot)

Kmer1 <- count(protein,wordsize=3,alphabet=aa)
Kmer2 <- count(protein,wordsize=4,alphabet=aa)
Kmer3 <- count(protein,wordsize=5,alphabet=aa)


str(Kmer1)
head(Kmer1)
aa[aa !="*"]

myfreq1 <- count(protein,wordsize=3,alphabet=aa,freq=TRUE)
myfreq2 <- count(protein,wordsize=4,alphabet=aa,freq=TRUE)
myfreq3 <- count(protein,wordsize=5,alphabet=aa,freq=TRUE)
head(myfreq1,n=10) #Kmer profiling for Kmers length 3
tail(myfreq1,n=10)

head(myfreq2,n=10) #Kmer profiling for Kmers length 4
tail(myfreq2,n=10)

head(myfreq3,n=10) #Kmer profiling for Kmers length 5
tail(myfreq3,n=10)

#Kmer profiling for S.bacterium
protein2 <- unlist(prot2)

Kmer4<- count(protein2,wordsize =3,alphabet=aa2)
Kmer5<- count(protein2,wordsize = 4,alphabet=aa2)
Kmer6<- count(protein2,wordsize = 5,alphabet=aa2)

str(Kmer4)
head(Kmer4)
aa[aa !="*"]

myfreq4<- count(protein2,wordsize =3,alphabet=aa2,freq = TRUE)
myfreq5 <- count(protein2,wordsize = 4,alphabet=aa2,freq=TRUE)
myfreq6 <- count(protein2,wordsize = 5,alphabet=aa2,freq=TRUE)
head(myfreq4,n=10) #Kmer profiling for Kmers length 3
tail(myfreq4,n=10)

head(myfreq5,n=10) #Kmer profiling for Kmers length 4
tail(myfreq5,n=10)

head(myfreq6,n=10) #Kmer profiling for Kmers length 5
tail(myfreq6,n=10)

plot(myfreq1,main="Kmer Profiling for E.Coli",xlab = "Kmers of Length 3",ylab="Frequency")
plot(myfreq2,main="Kmer Profiling for E.Coli",xlab = "Kmers of Length 4",ylab="Frequency")
plot(myfreq3,main="Kmer Profiling for E.Coli",xlab = "Kmers of Length 5",ylab="Frequency")
plot(myfreq4,main="Kmer Profiling for S.bacterium",xlab = "Kmers of Length 3",ylab="Frequency")
plot(myfreq5,main="Kmer Profiling for S.bacterium",xlab = "Kmers of Length 4",ylab="Frequency")
plot(myfreq6,main="Kmer Profiling for S.bacterium",xlab = "Kmers of Length 5",ylab="Frequency")

