Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
set.seed(220)
setwd("/home/bsb/breast/cox/")
library(data.table)
ex=fread('GSE202203_TPM_Raw_gene_3207.tsv.gz')
genename=ex$X
save(genename,file='genename.RData')
