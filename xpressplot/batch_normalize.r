#!/usr/bin/env Rscript

# Control batch effects for prep, chips, etc

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sva", version = "3.8")
BiocManager::install("bladderbatch", version = "3.8")
BiocManager::install("limma", version = "3.8")

source("http://bioconductor.org/biocLite.R")
biocLite("pamr")

library(sva)
library(bladderbatch)
library(pamr)
library(limma)

# Get arguments
# args[1] = dataframe
# args[2] = batch info
# args[5] = output file name
args = commandArgs(trailingOnly=TRUE)

# Import expression matrix and sample info
data = read.csv(args[1], sep='\t', row.names=1)
data_m = as.matrix(data)
info = read.csv(args[2], sep='\t')

# Identify batch effect to normalize on
names(info) <- tolower(names(info))
batch = info$batch

# Prepare model
modcombat = model.matrix(~1, data=info)

# Perform combat normalization on known batch parameter
combat_data = ComBat(dat=data_m, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
write.table(as.data.frame(combat_data), file=args[3], sep='\t', col.names=T, row.names=T)
