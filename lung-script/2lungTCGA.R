Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("/home/bsb/lung/TCGAlung/")
set.seed(2023)
library(data.table)
library(dplyr)
library(tidyverse)
### rnaseq data
count=fread("./LUAD/TCGA-LUAD.htseq_counts.tsv.gz",header = T, sep = '\t',data.table = F)
View(count)
dim(count)
count2=fread("./LUSC//TCGA-LUSC.htseq_counts.tsv.gz",header = T, sep = '\t',data.table = F)
View(count2)
dim(count2)
##probe data
pro=fread("./LUAD/gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
View(pro)
pro=pro[,c(1,2)]
pro2=fread("./LUSC/gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
View(pro2)
pro2=pro2[,c(1,2)]
###merge
countpro=merge(pro,count,by.y ="Ensembl_ID",by.x = "id" )
View(countpro)
exp=countpro[,-c(1)]
View(exp)
countpro2=merge(pro2,count2,by.y ="Ensembl_ID",by.x = "id" )
View(countpro2)
exp2=countpro2[,-c(1)]
View(exp2)
#limma average
library(limma)
exp=as.data.frame(avereps(exp[,-1],ID=exp$gene))
dim(exp)
exp2=as.data.frame(avereps(exp2[,-1],ID=exp2$gene))
dim(exp2)
#get parimary count
explog=2^exp-1
explog2=2^exp2-1
#see ditribution 
df1=explog[,1:20]
df2=explog2[,1:20]
df3=cbind(df1,df2)
group_list=c(rep('control',20),rep('cancer',20))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="control")
boxplot(df3,outline=FALSE, notch=T,col=group_list, las=2)
explog3=merge(explog,explog2,by.x=0,by.y=0)
rownames(explog3)=explog3[,1]
explog3=explog3[,-1]
#read colname
metadata <- data.frame(TCGA_id =colnames(explog3))
table(substring(metadata$TCGA_id,14,15))
sample <- ifelse(substring(metadata$TCGA_id,14,15)=="11","normal","cancer")
#factor needed for dds
metadata$sample <- as.factor(sample)
view(metadata)
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=round(explog3), 
                             colData=metadata, 
                             design=~sample)
dds <- DESeq(dds)


vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, "sample")

mRNA_exprSet <- as.data.frame(assay(vsd))
save(dds,file="mRNA_exprSet_dds_sample.RData")
save(mRNA_exprSet,file='mRNA_exprset_norm.RData')
write.csv(mRNA_exprSet, file = "mRNA_exprset_norm.csv",sep=",", row.names = T,quote = F)


expr_for_diff <- results(dds, tidy=TRUE)
#single gene analyse
plot_data<-data.frame(mRNA_exprSet['ITLN1',])
plot_data<-t(plot_data)
plot_data<-data.frame(plot_data)
plot_data <- cbind(plot_data,sample=metadata$sample)
#
library(ggplot2)
p<-ggplot(plot_data,aes(x=sample,y=ITLN1,fill=sample))+
  geom_boxplot()+
  theme_classic()+
  ggpubr::stat_compare_means(color="red")+
  ggsci::scale_fill_jco()+
  theme(legend.position = "none")+
  ylab("ITLN1 counts")
p
normalmarker=subset(expr_for_diff,log2FoldChange>0)
cancermarker=subset(expr_for_diff,log2FoldChange<0)
save(normalmarker,file='normalmarker.RData')
save(cancermarker,file='cancermarker.RData')