Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(GEOquery)
library(dplyr)
library(tidyverse)
setwd("/home/bsb/breast/bulk/GSE65194/")
set.seed(2023)
gset = getGEO('GSE65194', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#express matrix
exp <- exprs(gset)
View(exp)
min_val <- min(exp)
exp=as.data.frame(exp)
exp=exp+1.23
min_val <- min(exp)
#Clinical information
pdata=pData(gset)
View(pdata)
#Platform Information
gpl=fData(gset)
View(gpl)
gpl=gpl[,c(1,11)]
gpl=gpl[gpl[,"Gene Symbol"]!="",]
#merge
exp.pl=merge(gpl,exp,by.x=1,by.y=0)
x=exp.pl$`Gene Symbol`
a1=strsplit(x,split = " /// ",fixed = T) #exact match
gene.all = sapply(a1,function(x){x[1]})
exp.pl$`Gene Symbol`=gene.all
exp.pl=exp.pl[,-c(1)]
library(limma)
exp.pl=as.data.frame(avereps(exp.pl[,-1],ID=exp.pl$`Gene Symbol`))
##group info
View(pdata)
table(pdata$`sample_group:ch1`)
pdata=pdata[,c(2,40)]
colnames(pdata)[2]='group'
normal=filter(pdata,group=='Healthy')
TNBC=filter(pdata,group=='TNBC')
CellLine=filter(pdata,group=='CellLine')
luminalA=filter(pdata,group=='Luminal A')
luminalB=filter(pdata,group=='Luminal B')
Her2=filter(pdata,group=='Her2')
exp2=exp.pl[,c(normal[,c(1)],TNBC[,c(1)],CellLine[,c(1)],luminalA[,c(1)],luminalB[,c(1)],Her2[,c(1)])]
##control£¬latter is treat
group_list=c(rep('control',11),rep('cancer',167))
##turn as foctor
group_list=factor(group_list)
## limit order
group_list <- relevel(group_list, ref="control")
boxplot(exp2,outline=FALSE, notch=T,col=group_list, las=2)
library(RColorBrewer)
colors<-brewer.pal(25,"Set3") #produce a set of color
colors
boxplot(exp2,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,10)) 
##Build a matrix for differently express analysis
design=model.matrix(~ group_list)  #~followed factor£¬can¡®t be vector¡£0~ or ~£¬not add 0 shows first colum as control£¬compared with the second col
View(design)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp2)
##lmFit()£ºLinear fitting model construction
fit=lmFit(exp2,design)
##eBayes()
fit=eBayes(fit) 
##topTable()
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf)
View(allDiff)
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
#csv
write.csv(allDiff,file = "allDiff.csv")
#single gene analyse
example=exp2[c('S100P','ADIPOQ'),]
View(example)
mean(as.matrix(example[1,c(1:11)]))
mean(as.matrix(example[1,c(11:178)]))
mean(as.matrix(example[2,c(1:11)]))
mean(as.matrix(example[2,c(11:178)]))
tumormarker4=subset(allDiff,logFC>0)
normalmarker4=subset(allDiff,logFC<0)
save(normalmarker4,file='normalmarker4.RData')
save(tumormarker4,file='cancermarker4.RData')
