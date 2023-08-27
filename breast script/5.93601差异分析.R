Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("/home/bsb/breast/bulk/GSE93601/")
set.seed(2023)
gset = getGEO('GSE93601', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#express matrix
exp <- exprs(gset)
View(exp)
#Clinical information
pdata=pData(gset)
View(pdata)
#Platform Information
gpl=fData(gset)
View(gpl)
gpl=gpl[,c(1,2)]
gpl=gpl[gpl[,"ENTREZ_GENE_ID"]!="",]
### turn entrez id to gene symbol
library(org.Hs.eg.db)
###see all ID
keytypes(org.Hs.eg.db)
x1=gpl$ENTREZ_GENE_ID
x1=as.character(x1)
x2=AnnotationDbi::select(org.Hs.eg.db, keys=x1, columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID")#keys correspond keytype
View(x2)
gpl=merge(gpl,x2,by.x=2,by.y=1)
gpl2=gpl[,c(2,3)]
##merge
exp=as.data.frame(exp)
exp.pl=merge(gpl2,exp,by.x=1,by.y=0)
exp.pl=exp.pl[,-c(1)]
library(limma)
exp.pl=as.data.frame(avereps(exp.pl[,-1],ID=exp.pl$SYMBOL))
##group info
View(pdata)
pdata=pdata[,c(2,10)]
normal=filter(pdata,characteristics_ch1=="tissue: breast tumor-adjacent normal")
cancer=filter(pdata,characteristics_ch1=="tissue: breast tumor")
exp2=exp.pl[,c(normal[,c(1)],cancer[,c(1)])]
##create the grouped information, front is control, the latter is treat
group_list=c(rep('control',508),rep('T',602))
##turn to factor
group_list=factor(group_list)
##limit the order
group_list <- relevel(group_list, ref="control")
exp3=exp2[,c(1:10,1000:1010)]
boxplot(exp3,outline=FALSE, notch=T,col=group_list, las=2)
exp4=normalizeBetweenArrays(exp2)
exp5=exp4[,c(1:10,1000:1010)]
boxplot(exp5,outline=FALSE, notch=T,col=group_list, las=2)
library(RColorBrewer)
colors<-brewer.pal(1100,"Set3") #produce set of color
colors
boxplot(exp5,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,10)) 
##Build a matrix for differently express analysis
design=model.matrix(~ group_list) #~followed factor£¬can¡®t be vector¡£0~ or ~£¬not add 0 shows first colum as control£¬compared with the second col
View(design)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp2)
##lmFit()£ºLinear fitting model construction
fit=lmFit(exp4,design)
##eBayes()
fit=eBayes(fit) 
##topTable()
##
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) 
View(allDiff)
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
#csv
write.csv(allDiff,file = "allDiff.csv")
#single gene analyse
example=exp2[c('FN1','FABP4'),]
View(example)
mean(as.matrix(example[1,c(1:508)]))
mean(as.matrix(example[1,c(509:1100)]))
mean(as.matrix(example[2,c(1:508)]))
mean(as.matrix(example[2,c(509:1100)]))
tumormarker3=subset(allDiff,logFC>0)
normalmarker3=subset(allDiff,logFC<0)
save(normalmarker3,file='normalmarker3.RData')
save(tumormarker3,file='cancermarker3.RData')
