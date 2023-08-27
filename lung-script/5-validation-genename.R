#validation
setwd("/home/bsb/lung/GSE72094/")
library(GEOquery)
gset= getGEO('GSE72094', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#
pdata=pData(gset)
View(pdata[,50:70])
pdata=pdata[,c(2,20,19)]
a1=strsplit(pdata$characteristics_ch1.10,split = ": ",fixed = T) 
t1 = sapply(a1,function(x){x[2]})
pdata$fultime=t1
b1=strsplit(pdata$characteristics_ch1.9,split = ": ",fixed = T) 
s1= sapply(b1,function(x){x[2]})
pdata$fulstate=s1
pdata=pdata[,-c(2:3)]
pdata=subset(pdata,fulstate!='NA'&fultime!='NA')
pdata$fultime=as.numeric(pdata$fultime)
pdata[,"fultime"]=pdata[,"fultime"]/365
table(pdata$fulstate)
frame=data.frame(fulstate=c('Alive','Dead'),state=c(0,1))
pdata2=merge(pdata,frame,by.x=3,by.y=1)
pdata2=pdata2[,c(2,3,4,1)]
table(pdata2$state)
pdata2=pdata2[,-4]
#
exprs=exprs(gset)
gpl=fData(gset)
gpl=gpl[,c(1,4)]
exgpl=merge(gpl,exprs,by.x=1,by.y=0)
exgpl=exgpl[,-1]
genename=exgpl$GeneSymbol
save(genename,file='genename.RData')

