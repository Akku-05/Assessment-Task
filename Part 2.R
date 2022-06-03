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

# Calculating the Coding DNA present in both the organisms
len1 <- as.numeric(summary(cds)[,1])
len2 <- as.numeric(summary(cds2)[,1])

sum(len1)
sum(len2)


length(unlist(cds[1:3]))
length(unlist(cds2[1:3]))



dna1 <- unlist(cds)
dna2 <- unlist(cds2)
head(dna1)
head(dna2)

#QUESTION 2
# Mean and Median for E.coli
head(summary(cds))
str(summary(cds))

head( summary(cds)[,1])
len3 <- sapply(X=cds,FUN=length) #calculating the
median(len3)
mean(len3)
#Mean and median for S.bacterium
head(summary(cds2))
str(summary(cds2))
head( summary(cds2)[,1])
len4 <- sapply(X=cds2,FUN=length)
median(len4)
mean(len4)

#Boxplot of coding sequences
boxplot(len3,len4,main= "Length of Coding Sequences", at=c(1,2),
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

# Bar plot fpr DNA composition
barplot(dna_composition,xlab="nucleotides",ylab="frequency", main="E.Coli DNA composition")

barplot(dna_composition2,xlab="nucleotides",ylab="frequency", main="S.bacterium DNA composition")



#Calculating protein sequence for E.coli
translate(cds[[1]])
prot <- lapply(cds,translate)
prot[[1]]
aa <-unique(prot[[2]])
aa <- aa[aa !="*"]

protein <- unlist(prot)

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

protein_composition2 <- count(protein2,alphabet=aa2,wordsize=1)


protseq2<- count(protein2,alphabet=aa2,wordsize=1)

#Barplots for Protein Composition

barplot(protseq,xlab="Protein composition",ylab="frequency", main="E.Coli Protein composition")
barplot(protseq2,xlab="Protein composition",ylab="frequency", main="S.bacterium Protein composition")





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



head(Kmer1)
aa[aa !="*"]

myfreq1 <- count(protein,wordsize=3,alphabet=aa,freq=TRUE)
sfreq1 <- sort(myfreq1)

myfreq2 <- count(protein,wordsize=4,alphabet=aa,freq=TRUE)
sfreq2 <- sort(myfreq2)

myfreq3 <- count(protein,wordsize=5,alphabet=aa,freq=TRUE)
sfreq3 <- sort(myfreq3)



#10 overexpressed kmers of length 3 for E.coli
head(sfreq1,n=10) #Kmer profiling for Kmers length 3

#10 underexpressed kmers of length 3 for E.coli
tail(sfreq1,n=10)

#10 overexpressed kmers of length 4 for E.coli
head(sfreq2,n=10) #Kmer profiling for Kmers length 4

#10 underexpressed kmers of length 5 for E.coli
tail(sfreq2,n=10)

#10 overexpressed kmers of length 3 for E.coli
head(sfreq3,n=10) #Kmer profiling for Kmers length 5

#10 underexpressed kmers of length 5 for E.coli
tail(sfreq3,n=10)

#Kmer profiling for S.bacterium
protein2 <- unlist(prot2)

Kmer4<- count(protein2,wordsize =3,alphabet=aa2)
Kmer5<- count(protein2,wordsize = 4,alphabet=aa2)
Kmer6<- count(protein2,wordsize = 5,alphabet=aa2)


head(Kmer4)
aa[aa !="*"]

myfreq4<- count(protein2,wordsize =3,alphabet=aa2,freq = TRUE)
sfreq4 <- sort(myfreq4)

myfreq5 <- count(protein2,wordsize = 4,alphabet=aa2,freq=TRUE)
sfreq5 <- sort(myfreq5)

myfreq6 <- count(protein2,wordsize = 5,alphabet=aa2,freq=TRUE)
sfreq6 <- sort(myfreq6)


#10 under represented kmers of length 3 for S.bacterium
head(myfreq4,n=10) #Kmer profiling for Kmers length 3

#10 over represented kmers of length 3 for S.bacterium
tail(sfreq4,n=10)

#10 under represented kmers of length 4 for S.bacterium
head(sfreq5,n=10) #Kmer profiling for Kmers length 4

#10 over represented kmers of length 4 for S.bacterium
tail(sfreq5,n=10)

#10 under represented kmers of length 5 for S.bacterium
head(sfreq6,n=10) #Kmer profiling for Kmers length 5

#10 over represented kmers of length 4 for S.bacterium
tail(sfreq6,n=10)

plot(sfreq1,main="Frequency of expression of kmers for E.coli",xlab = "Kmers of Length 3",ylab="Frequency")

plot(sfreq4,main="Frequency of expression of kmers for S.bacterium",xlab = "Kmers of Length 3",ylab="Frequency")
