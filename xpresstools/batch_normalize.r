#Control batch effects for prep, chips, etc

#Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sva", version = "3.8")
BiocManager::install("bladderbatch", version = "3.8")
BiocManager::install("limma", version = "3.8")

source("http://bioconductor.org/biocLite.R")
biocLite("pamr")

#Initialize workspace
library(sva)
library(bladderbatch)
library(pamr)
library(limma)

#args[1] = dataframe
#args[2] = batch info
#args[3] = dataframe sep
#args[4] = batch sep
#args[5] = output file name

#Import expression matrix and sample info
data = read.csv(args[1], sep=args[3], row.names=1)
data_m = as.matrix(data)
info = read.csv(args[2],sep=args[4])

#Identify batch effect to normalize on
batch = info$Chip

#Prepare model
modcombat = model.matrix(~1, data=info)

#Perform combat normalization on known batch parameter
combat_data = ComBat(dat=data_m, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
write.csv(combat_data, file = args[5])
