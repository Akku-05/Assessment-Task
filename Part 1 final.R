#dowloading growth data and gene expression data

URL= "https://raw.githubusercontent.com/markziemann/SLE712_files/master/assessment_task3/bioinfo_asst3_part1_files/gene_expression.tsv"
gene= "gene_expression.tsv"
download.file(URL,destfile=gene)

URL= "https://raw.githubusercontent.com/markziemann/SLE712_files/master/assessment_task3/bioinfo_asst3_part1_files/growth_data.csv"
growth= "growth_data.csv"
download.file(URL,destfile=growth)

GE <- read.delim("gene_expression.tsv")
GD <-read.csv("growth_data.csv")
str(GE)
#Changing row names to gene expression name.
row.names(GE) 
colnames(GE)
df <- data.frame(GE)
rownames(df) <- df$Name_Description
df$Name_Description=NULL

# subsetting first 6 genes

df[1:6,]

#formatting table
library(kableExtra)
df[1:6,]

head(df) %>%
  kbl(caption = "gene expression data table") %>%
  kable_paper(full_width = F, html_font = "Arial")

#question 2

colMeans(df[1:6,1:3])

rowMeans(df[1:6,1:3])

df$mean <- rowMeans(df)
A <- df[1:6,]
GeneMean <- rowMeans(A)
head(GeneMean) %>%
  kbl(caption = "gene expression mean table") %>%
  kable_paper(full_width = F, html_font = "Arial")

#question 3

B <-order(df$mean, decreasing=TRUE)
HighestG <-B[1:10]

#question 4
df$mean < 10
Mean_filtered <- df[which(df$mean < 10),]
subset(df,df$mean < 10)

#question 5

hist(df$mean, xlab = "mean", ylab = "expression", main = "Gene Expression")

#question 6 (already downloaded file)

GD
colnames(GD)

#question 7
NE<-subset(GD,Site =="northeast")
SW <- subset(GD,Site =="southwest")
NE2005 <- NE$Circumf_2005_cm
NE2020 <- NE$Circumf_2020_cm
SW2005 <- SW$Circumf_2005_cm
SW2020 <- SW$Circumf_2020_cm
mean(NE2005)
mean(NE2020)
mean(SW2005)
mean(SW2020)
a <- sd(NE2005)
b <-sd(NE2020)
c <-sd(SW2005)
d <-sd(SW2020)
boxplot(a,b,c,d, main= "Growth SD", at=c(1,2,3,4),
        names= c("NE2005","NE2020","SW2005","
                 SW2020"), ylab="nshvd")


# question 8
boxplot(a,b,c,d, main= "Growth SD", at=c(1,2,3,4),
        names= c("NE2005","NE2020","SW2005","
                 SW2020"), ylab="nshvd")
grid()

#question9
NE2010 <- NE$Circumf_2010_cm

SW2010 <- SW$Circumf_2010_cm
growthrateNE <- ((NE2020-NE2010)/NE2010) #growth rate = ((present-past)growth/past growth)
growthrateSW <- ((SW2020-SW2010)/SW2010)
SWMean10 <- mean(growthrateSW)
NEMean10 <-mean(growthrateNE)
NEMean10      #Mean for northeast region for past 10 years
SWMean10     #Mean for southwest region for past 10 years
MeanGrowth <- mean(growthrateNE+growthrateSW)
MeanGrowth      #Mean growth rate of both the site for past 10 years


#Question 10

res <- t.test (growthrateNE,growthrateSW)
p.value <-res$p.value
HEADER <- paste("P Value:", signif(p.value,3))

mtext(HEADER)

p.value2 <-wilcox.test(growthrateNE,growthrateSW)
HEADER <- paste("P Value:", signif(p.value,3))

mtext(HEADER)

