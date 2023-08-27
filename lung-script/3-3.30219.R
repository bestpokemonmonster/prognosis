Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("/home/bsb/lung/GSE30219/")
set.seed(2023)
gset = getGEO('GSE30219', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#matrix
exp <- exprs(gset)
View(exp)
min(exp)
#pdata
pdata=pData(gset)
View(pdata)
pdata=pdata[,c(2,8)]
#gpl
gpl=fData(gset)
View(gpl)
gpl=gpl[,c(1,11)]
#merge
exgpl=merge(gpl,exp,by.x = 1,by.y=0)
x=exgpl$`Gene Symbol`
a1=strsplit(x,split = " /// ",fixed = T) 
gene.all = sapply(a1,function(x){x[1]})
exgpl$`Gene Symbol`=gene.all
exgpl=exgpl[,-c(1)]
colnames(exgpl)[1]='Symbol'
missing_coords <- which(is.na(exgpl), arr.ind = TRUE)
print(missing_coords)
exgpl=subset(exgpl,Symbol!='NA')
library(limma)
exgpl=as.data.frame(avereps(exgpl[,-1],ID=exgpl$Symbol))
#
table(pdata$source_name_ch1)
tumour=subset(pdata,source_name_ch1=='Lung Tumour')
normal=subset(pdata,source_name_ch1=='Non Tumoral Lung')
exp=exgpl[,c(rownames(normal),rownames(tumour))]

##
group_list=c(rep('control',14),rep('cancer',293))
##
group_list=factor(group_list)
## 
group_list <- relevel(group_list, ref="control")
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
##
design=model.matrix(~ group_list) 
View(design)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp)
##lmFit()£º
fit=lmFit(exp,design)
##eBayes()
fit=eBayes(fit) 
##topTable()
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) 
View(allDiff)
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
#csv
write.csv(allDiff,file = "allDiff.csv")
#
example=exp[c('SCGB1A1','MMP12'),]
View(example)
mean(as.matrix(example[1,c(1:14)]))
mean(as.matrix(example[1,c(15:307)]))
mean(as.matrix(example[2,c(1:14)]))
mean(as.matrix(example[2,c(15:307)]))
tumormarker4=subset(allDiff,logFC>0)
normalmarker4=subset(allDiff,logFC<0)
save(normalmarker4,file='normalmarker4.RData')
save(tumormarker4,file='cancermarker4.RData')







