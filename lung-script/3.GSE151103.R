Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("/home/bsb/lung/GSE151103/")
set.seed(2023)
gset = getGEO('GSE151103', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#matrix
exp <- exprs(gset)
View(exp)
# remove na
cleaned_data <- exp[complete.cases(exp), , drop = FALSE]
print(cleaned_data)
#print place of na
missing_coords <- which(is.na(cleaned_data), arr.ind = TRUE)
print(missing_coords)
#see min 
min_val <- min(cleaned_data)
#pdata
pdata=pData(gset)
View(pdata)
#gpl
gpl=fData(gset)
View(gpl)
gpl=gpl[,c(1,10)]
#merge
exp.pl=merge(gpl,cleaned_data,by.x=1,by.y=0)
x=exp.pl$gene_assignment
a1=strsplit(x,split = " // ",fixed = T) 
gene.all = sapply(a1,function(x){x[2]})
exp.pl$gene_assignment=gene.all
exp.pl=exp.pl[,-c(1)]
colnames(exp.pl)[1]='Symbol'
missing_coords <- which(is.na(exp.pl), arr.ind = TRUE)
print(missing_coords)
exp.pl=subset(exp.pl,Symbol!='NA')
library(limma)
exp.pl=as.data.frame(avereps(exp.pl[,-1],ID=exp.pl$Symbol))
##group data
View(pdata)
pdata=pdata[,c(2,12)]
table(pdata$characteristics_ch1.2)
colnames(pdata)[2]='group'
normal=filter(pdata,group=='tissue: Normal')
tumor=filter(pdata,group=='tissue: Tumor')
exp2=exp.pl[,c(normal[,c(1)],tumor[,c(1)])]
colnames(exp2)
##
group_list=c(rep('control',172),rep('cancer',188))
##turn to factor
group_list=factor(group_list)
## limit order
group_list <- relevel(group_list, ref="control")
boxplot(exp2,outline=FALSE, notch=T,col=group_list, las=2)
exp3=normalizeBetweenArrays(exp2)
boxplot(exp3,outline=FALSE, notch=T,col=group_list, las=2) 
##make matrix for differently analyse
design=model.matrix(~ group_list) 
View(design)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp2)
##lmFit()£º
fit=lmFit(exp3,design)
##eBayes()
fit=eBayes(fit) 
##topTable()
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) 
View(allDiff)
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
write.csv(allDiff,file = "allDiff.csv")
#single gene
example=exp3[c('SLC6A4','SPP1'),]
View(example)
mean(as.matrix(example[1,c(1:172)]))
mean(as.matrix(example[1,c(173:360)]))
mean(as.matrix(example[2,c(1:172)]))
mean(as.matrix(example[2,c(173:360)]))
tumormarker2=subset(allDiff,logFC>0)
normalmarker2=subset(allDiff,logFC<0)
save(normalmarker2,file='normalmarker2.RData')
save(tumormarker2,file='cancermarker2.RData')