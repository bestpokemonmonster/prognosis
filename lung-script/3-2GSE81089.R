Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("/home/bsb/lung/GSE81089/")
set.seed(2023)
gset = getGEO('GSE81089', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#express
library(data.table)
exp=fread('GSE81089_readcounts_featurecounts.tsv.gz')
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
x1=exp$Ensembl_gene_id
x1=as.character(x1)
x2=AnnotationDbi::select(org.Hs.eg.db, keys=x1, columns=c("SYMBOL"), keytype="ENSEMBL")
View(x2)
x3=subset(x2,SYMBOL!='NA')
exp2=merge(x3,exp,by.x=1,by.y=1)
exp2=exp2[,-1]
exp2=distinct(exp2,SYMBOL,.keep_all=T)
rownames(exp2)=exp2$SYMBOL
exp3=exp2[,-1]
name=data.frame(a=colnames(exp3))
name$status <- ifelse(substring(name$a,5)=="N","normal","cancer")
table(name$status)
normal=subset(name,status=='normal')
cancer=subset(name,status=='cancer')
normalname=normal[,1]
cancername=cancer[,1]
exp4=exp3[,c(normalname,cancername)]
colnames(exp4)
head(exp4)[1:5]
condition <- factor(c(rep("normal",19),rep("cancer",199)), levels = c("normal","cancer"))
condition
coldata<-data.frame(row.names=colnames(exp4),condition)
coldata
#dds
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=exp4,colData=coldata, design= ~ condition)
dds <- DESeq(dds)
dds


vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, "condition")
mRNA_exprSet <- as.data.frame(assay(vsd))

res = results(dds, contrast=c("condition", "normal", "cancer"))
res = res[order(res$pvalue),]
head(res)
summary(res)
res=as.data.frame(res)

can=exp3[,cancername]
nor=exp3[,normalname]
c1=(can['SLC6A4',])
mean(as.numeric(c1))
n1=(nor['SLC6A4',])
mean(as.numeric(n1))

plot_data<-data.frame(mRNA_exprSet["MAGEA11",])
plot_data<-t(plot_data)
plot_data<-data.frame(plot_data)
plot_data <- merge(plot_data,name,by.x=0,by.y=1)

library(ggplot2)
p<-ggplot(plot_data,aes(x=status,y=MAGEA11,fill=status))+
  geom_boxplot()+
  theme_classic()+
  ggpubr::stat_compare_means(color="red")+
  ggsci::scale_fill_jco()+
  theme(legend.position = "none")+
  ylab("counts")
print(p)

normalmarker3=subset(res,log2FoldChange>0)                 
cancermarker3=subset(res,log2FoldChange<0)
save(normalmarker3,file='normalmarker3.RData')
save(cancermarker3,file='cancermarker3.RData')