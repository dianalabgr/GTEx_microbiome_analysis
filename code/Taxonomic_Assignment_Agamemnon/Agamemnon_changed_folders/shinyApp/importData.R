path <- paste(getwd(), "/filipos_results", sep = "") 
setwd(path)

library("Biobase")
library("metagenomeSeq")

i = list.files()[1]
x = read.csv(i, sep = "\t", stringsAsFactors = FALSE)
pd = x[, -ncol(x)]

pd = data.frame(pd)
rownames(pd) <- pd$TaxID
pd = AnnotatedDataFrame(pd)

mats = sapply(list.files(), function(i)
{
			x = read.csv(i, sep="\t", stringsAsFactors = FALSE)[, 10]
})

rownames(mats) = rownames(pd)
colnames(mats) = gsub("\\.tab", "", colnames(mats))

phenotypes = read.csv("../phenoData/phenotypes_filip.tab", sep = " ", stringsAsFactors = FALSE)
rownames(phenotypes) <- phenotypes$External_ID

hmp_metag = newMRexperiment(mats, featureData = pd)
phenoData(hmp_metag) = AnnotatedDataFrame(phenotypes)
hmp_metag = filterData(hmp_metag, present = 1, depth = 1)

