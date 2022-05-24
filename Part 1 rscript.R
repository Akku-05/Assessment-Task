url <- "https://raw.githubusercontent.com/markziemann/SLE712_files/master/assessment_task3/bioinfo_asst3_part1_files/gene_expression.tsv"
destfile <- "gene_expression.tsv"

download.file(url,destfile)

library (readr)
urlfile <- 

urlfile <-"https://github.com/markziemann/SLE712_files/blob/master/assessment_task3/bioinfo_asst3_part1_files/growth_data.csv"
destfile1 <- "growth_data.csv"

download.file(urlfile,destfile1) 

GE<-read.table(file="gene_expression.tsv", header= FALSE, fill= TRUE)
str(GE)
