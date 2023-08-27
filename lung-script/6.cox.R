Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("/home/bsb/lung/cox/")
library(data.table)
set.seed(148)
#get gene list
load('../robust/cancersign.RData')
load('../robust/normalsign.RData')
load('../sclung/df_new.RData')
library(dplyr)
cancersign=cancersign[c(1:100),]
normalsign=normalsign[c(1:100),]
bulksign=rbind(cancersign,normalsign)
scsign=subset(df_new,power>0.71)
gene=bulksign[,1]
gene2=rownames(scsign)
gene3=c(gene,gene2)
gene4=unique(gene3)
#
load('../TCGAlung/mRNA_exprset_norm.RData')
gene5=rownames(mRNA_exprSet)
load('../GSE72094/genename.RData')
gene6=intersect(gene5,genename)
gene7=intersect(gene4,gene6)
exp=mRNA_exprSet[gene7,]
exp=as.data.frame(t(exp))
#survive data
surviveluad=fread('../TCGAlung/LUAD/TCGA-LUAD.survival.tsv')
survivelusc=fread('../TCGAlung/LUSC/TCGA-LUSC.survival.tsv')
survive=rbind(surviveluad,survivelusc)
#
name=survive[,1]
table(substring(name$sample,14,15))
survive$status <- ifelse(substring(name$sample,14,15)=="11","normal","cancer")
survive=subset(survive,status=='cancer')
name=survive[,1]
table(substring(name$sample,14,15))
#remove repeated patient
colnames(survive)[3]='patient'
library(dplyr)
survive2=distinct(survive,patient,.keep_all = T)
survive3=survive2[,c(1,4,2)]
#merge
survive3=as.data.frame(survive3)
exp2=merge(survive3,exp,by.x=1,by.y=0)
exp2[,"OS.time"]=exp2[,"OS.time"]/365
#group
library(caret)
group=createDataPartition(y=exp2[,2],p=0.7,list=F)
train=exp2[group,]
test=exp2[-group,]
#
prop.table(table(train$OS))
prop.table(table(test$OS))
prop.table(table(train$OS.time))
prop.table(table(test$OS.time))
#cox
library(survival)
coxscore=data.frame()
for (i in colnames(train)[4:ncol(train)]) {
cox<-coxph(Surv(OS.time,OS)~train[,i],data=train)
coxsummary=summary(cox)
coxscore=rbind(coxscore,cbind(gene=i,HR=coxsummary$coefficients[,'exp(coef)'],
z=coxsummary$coefficients[,'z'],
pvalue=coxsummary$coefficients[,'Pr(>|z|)']))
}
coxscoreall=coxscore
coxscore=subset(coxscore,pvalue<0.02)
intersect(coxscore$gene,gene)
intersect(coxscore$gene,gene2)
save(coxscore,file='coxscore.RData')
save(train,file='train.RData')
save(test,file='test.RData')
#multicox
coxgene=coxscore[,1]
colname=colnames(train)[1:3]
colname2=c(colname,coxgene)
multiexp=train[,colname2]
rownames(multiexp)=multiexp$sample
multiexp=multiexp[,-1]
#multicox
cox2<-coxph(Surv(OS.time,OS)~.,data=multiexp)
cox2=step(cox2,direction = 'both')
cox2
cox2$assign
b=cox2$assign
b=as.character(names(b))
b
length(b)

intersect(gene,b)
intersect(gene2,b)
can=cancersign[,1]
intersect(can,b)
nor=normalsign[,1]
intersect(nor,b)
sccan=subset(df_new,avg_log2FC<0)
sccan2=rownames(sccan)
intersect(sccan2,b)
scnor=subset(df_new,avg_log2FC>0)
scnor2=rownames(scnor)
intersect(scnor2,b)

#validation
rownames(test)=test[,1]
test=test[,-1]
riskscore=predict(cox2,type='risk',newdata = test)
risk=as.vector(ifelse(riskscore>median(riskscore),'hign','low'))
write.table(cbind(id=rownames(cbind(test[,1:2],riskscore,risk)),cbind(test[,1:2],riskscore,risk)),file="risk0.txt",sep="\t",quote=F,row.names=F)
#survive
risk2=read.table("risk0.txt",header=T,sep="\t")
diff=survdiff(Surv(OS.time, OS) ~risk,data = risk2)
pValue=1-pchisq(diff$chisq,df=1)
diff
fit <- survfit(Surv(OS.time, OS) ~ risk, data = risk2)
summary(fit)    
pdf(file="survival.pdf")
plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (year)",ylab="survival rate",
     main=paste("survival curve (p=",'2.11e-4',")",sep=""),mark.time=T,
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
legend("topright", c("high risk", "low risk"), lty = 2:3, col=c("red","blue"))
dev.off()
#roc
library(survivalROC)
risk3=read.table("risk0.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="ROC.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=risk3$OS.time, status=risk3$OS, marker = risk3$riskscore, 
                predict.time =3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
abline(0,1)
dev.off()
m1=round(median(risk2$riskscore),3)



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
b=as.data.frame(b)
exgpl2=merge(b,exgpl,by.x=1,by.y=1)
library(limma)
exgpl3=as.data.frame(avereps(exgpl2[,-1],ID=exgpl2$b))
#
exgpl4=as.data.frame(t(exgpl3))
test2=merge(pdata2,exgpl4,by.x=1,by.y=0)
rownames(test2)=test2[,1]
test2=test2[,-1]
colnames(test2)[1:2]=c('fultime','fulstate')
#os
riskscore2=predict(cox2,type='risk',newdata = test2)
risk4=as.vector(ifelse(riskscore2>median(riskscore2),'hign','low'))
write.table(cbind(id=rownames(cbind(test2[,1:2],riskscore2,risk4)),cbind(test2[,1:2],riskscore2,risk4)),file="risk.txt",sep="\t",quote=F,row.names=F)
#survive
risk5=read.table("risk.txt",header=T,sep="\t")
diff2=survdiff(Surv(fultime, fulstate) ~risk4,data = risk5)
pValue2=1-pchisq(diff2$chisq,df=1)
diff2
fit2 <- survfit(Surv(fultime, fulstate) ~ risk4, data = risk5)
summary(fit2)    
pdf(file="../cox/survival2.pdf")
plot(fit2, lty = 2:3,col=c("red","blue"),xlab="time (year)",ylab="survival rate",
     main=paste("survival curve (p=", '1.79e-4',")",sep=""),mark.time=T,
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
legend("topright", c("high risk", "low risk"), lty = 2:3, col=c("red","blue"))
dev.off()
#roc
library(survivalROC)
risk6=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="../cox/ROC2.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc2=survivalROC(Stime=risk6$fultime, status=risk6$fulstate, marker = risk6$riskscore2, 
                 predict.time =5, method="KM")
plot(roc2$FP, roc2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc2$AUC,3),")"),
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
abline(0,1)
dev.off()
m2=round(median(risk5$riskscore2),3)