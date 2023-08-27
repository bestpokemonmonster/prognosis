Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("/home/bsb/breast/bulk/TCGA")
set.seed(2023)
library(data.table)
library(dplyr)
library(tidyverse)
###load rnaseq data
count=fread("TCGA-BRCA.htseq_counts.tsv.gz",header = T, sep = '\t',data.table = F)
View(count)
dim(count)
####probe data
pro=fread("gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
View(pro)
pro=pro[,c(1,2)]
View(pro)
###merge
countpro=merge(pro,count,by.y ="Ensembl_ID",by.x = "id" )
View(countpro)
exp=countpro[,-c(1)]
View(exp)
#limma to get average
library(limma)
exp=as.data.frame(avereps(exp[,-1],ID=exp$gene))
dim(exp)
#get primary count
explog=2^exp-1
#get colname
metadata <- data.frame(TCGA_id =colnames(explog))
table(substring(metadata$TCGA_id,14,15))
sample <- ifelse(substring(metadata$TCGA_id,14,15)=="11","normal","cancer")
#factor needed for dds
metadata$sample <- as.factor(sample)
view(metadata)
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=round(explog), 
                             colData=metadata, 
                             design=~sample)
dds <- DESeq(dds)

#if need express matrix to analyse PCA£¬WGCNA,CLUSTERING£¬require vst normalization
vsd <- vst(dds, blind = FALSE)
#see sample distribution
plotPCA(vsd, "sample")

#acquire normalized data,this step will filter out genes not fit the regulations
mRNA_exprSet <- as.data.frame(assay(vsd))
save(dds,file="mRNA_exprSet_dds_sample.RData")
save(mRNA_exprSet,file='mRNA_exprset_norm.RData')
write.csv(mRNA_exprSet, file = "mRNA_exprset_norm.csv",sep=",", row.names = T,quote = F)

#analyse differently express
expr_for_diff <- results(dds, tidy=TRUE)
#ensure which upregulate and downregulate
#single gene analyse
plot_data<-data.frame(mRNA_exprSet['CSAG1',])
plot_data<-t(plot_data)
plot_data<-data.frame(plot_data)
plot_data <- cbind(plot_data,sample=metadata$sample)
#visualization
library(ggplot2)
p<-ggplot(plot_data,aes(x=sample,y=CSAG1,fill=sample))+
  geom_boxplot()+
  theme_classic()+
  ggpubr::stat_compare_means(color="red")+
  ggsci::scale_fill_jco()+
  theme(legend.position = "none")+
  ylab("CSAG1 counts")
p
normalmarker=subset(expr_for_diff,log2FoldChange>0)
cancermarker=subset(expr_for_diff,log2FoldChange<0)
save(normalmarker,file='normalmarker.RData')
save(cancermarker,file='cancermarker.RData')










