#!/usr/bin/env Rscript
#Control batch effects for prep, chips, etc

#Install dependencies
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

#Get arguments
#args[1] = dataframe
#args[2] = sample info
#args[3] = output file name
#args[4] = DE equation
args = commandArgs(trailingOnly=TRUE)

#import counts_data
count_table <- read.table(args[1],sep='\t',header=TRUE,row.names=1)

#create conditions dataframe
sample_table <- read.table(text=readLines(args[2], warn = FALSE), header=TRUE, sep='\t')
names(sample_table) <- tolower(names(sample_table))

#run DESeq2 analysis on data
dds <- DESeqDataSetFromMatrix(countData = count_table, colData = sample_table, design = as.formula(paste('~', toLower(toString(args[4]), sep=''))))
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]

#write output to new file
write.table(as.data.frame(resOrdered),file=args[3],sep='\t',col.names=T,row.names=T)
