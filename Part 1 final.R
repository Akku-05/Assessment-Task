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

