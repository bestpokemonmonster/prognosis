Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("/home/bsb/breast/bulk/GSE109169/")
set.seed(2023)
gset = getGEO('GSE109169', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#express matrix
exp <- exprs(gset)
View(exp)
boxplot(exp)
#Clinical information
pdata=pData(gset)
View(pdata)
#Platform Information
gpl=fData(gset)
View(gpl)
gpl=gpl[,c(1,10)]
#merge
exp=as.data.frame(exp)
exp.pl=merge(gpl,exp,by.x=1,by.y=0)
x=exp.pl$gene_assignment
a1=strsplit(x,split = " // ",fixed = T) #T means exact match
gene.all = sapply(a1,function(x){x[2]})
exp.pl$`Gene Symbol`=gene.all
exp.pl=exp.pl[,c(53,1:52)]
exp.pl=exp.pl[,-c(2,3)]
library(limma)
exp.pl=as.data.frame(avereps(exp.pl[,-1],ID=exp.pl$`Gene Symbol`))
##View grouping information
View(pdata)
pdata=pdata[,c(2,8)]
normal=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)
cancer=normal+1
exp2=exp.pl[,c(normal,cancer)]
View(exp2)
##create the grouped information, front is control, the latter is treat
group_list=c(rep('control',25),rep('cancer',25))
##turn to factor
group_list=factor(group_list)
##limit the order
group_list <- relevel(group_list, ref="control")
boxplot(exp2,outline=FALSE, notch=T,col=group_list, las=2)
exp3=normalizeBetweenArrays(exp2)
boxplot(exp3,outline=FALSE, notch=T,col=group_list, las=2)
library(RColorBrewer)
colors<-brewer.pal(25,"Set3") #produce set of color
colors
boxplot(exp3,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,10)) 
##Build a matrix for differently express analysis
design=model.matrix(~ group_list) #~followed factor£¬can¡®t be vector¡£0~ or ~£¬not add 0 shows first colum as control£¬compared with the second col
View(design)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp2)
##lmFit()£ºLinear fitting model construction
fit=lmFit(exp3,design)
##eBayes()
fit=eBayes(fit) 
##topTable
##coef=2 mean design¡¯s second col--tumour£¬compared with normal
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) #fdr mean adj.p
View(allDiff)
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
#csv
write.csv(allDiff,file = "allDiff.csv")
#single gene analyse
example=exp2[c('CST1','FABP4'),]
View(example)
tumormarker2=subset(allDiff,logFC>0)
normalmarker2=subset(allDiff,logFC<0)
save(normalmarker2,file='normalmarker2.RData')
save(tumormarker2,file='cancermarker2.RData')
