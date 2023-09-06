Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
set.seed(220)
setwd("/home/bsb/breast/cox/")
library(data.table)
#gene list
load('../bulk/RobustRank/cancersign.RData')
load('../bulk/RobustRank/normalsign.RData')
load('../scRNA/SCTB/df_new.RData')
cancersign=cancersign[c(1:100),]
normalsign=normalsign[c(1:100),]
bulksign=rbind(cancersign,normalsign)
scsign=df_new
gene=bulksign[,1]
gene2=rownames(scsign)
gene3=c(gene,gene2)
intersect(gene,gene2)
gene4=unique(gene3)
#
c1=cancersign[,1]
n1=normalsign[,1]
sc1=rownames(scsign)
write.csv(c1,'breastcancergene.csv')
write.csv(n1,'breastnormalgene.csv')
write.csv(sc1,'breastscgene.csv')
#
load('../bulk/TCGA/mRNA_exprset_norm.RData')
gene5=rownames(mRNA_exprSet)
load('../cox/genename.RData')
gene6=intersect(gene5,genename)
gene7=intersect(gene4,gene6)
exp=mRNA_exprSet[gene7,]
exp=as.data.frame(t(exp))
#survival
survive=fread('./TCGA-BRCA.survival.tsv')
#see sample
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
#proportion
prop.table(table(train$OS))
prop.table(table(test$OS))
prop.table(table(train$OS.time))
prop.table(table(test$OS.time))
#cox of one gene
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
coxscore=subset(coxscore,pvalue<0.05)
save(coxscore,file='coxscore.RData')
save(train,file='train.RData')
save(test,file='test.RData')
#make multi-cox file
coxgene=coxscore[,1]
colname=colnames(train)[1:3]
colname2=c(colname,coxgene)
multiexp=train[,colname2]
rownames(multiexp)=multiexp$sample
multiexp=multiexp[,-1]
#multi-cox
cox2<-coxph(Surv(OS.time,OS)~.,data=multiexp)
cox2=step(cox2,direction = 'both')
cox2
cox2$assign
b=cox2$assign
b=as.character(names(b))
b
intersect(gene3,b)
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
coxscore2=subset(coxscore,gene%in%b)
#validation of tcga
rownames(test)=test[,1]
test=test[,-1]
riskscore=predict(cox2,type='risk',newdata = test)
risk=as.vector(ifelse(riskscore>median(riskscore),'hign','low'))
write.table(cbind(id=rownames(cbind(test[,1:2],riskscore,risk)),cbind(test[,1:2],riskscore,risk)),file="risk0.txt",sep="\t",quote=F,row.names=F)
#survive
risk2=read.table("risk0.txt",header=T,sep="\t")
diff=survdiff(Surv(OS.time, OS) ~risk,data = risk2)
pValue=1-pchisq(diff$chisq,df=1)
p1=signif(pValue,3)
diff
fit <- survfit(Surv(OS.time, OS) ~ risk, data = risk2)
summary(fit)    
pdf(file="16-survival.pdf")
plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (year)",ylab="survival rate",
     main=paste("survival curve (p=", '4.94e-3' ,")",sep=""),mark.time=T,
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
legend("topright", c("high risk", "low risk"), lty = 2:3, col=c("red","blue"))
dev.off()
#roc
library(survivalROC)
risk3=read.table("risk0.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="16-ROC.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=risk3$OS.time, status=risk3$OS, marker = risk3$riskscore, 
                predict.time =5, method="KM")
roc$AUC
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
abline(0,1)
dev.off()
m1=round(median(risk2$riskscore),3)


#validation
setwd("/home/bsb/breast/cox/")
library(GEOquery)
gset= getGEO('GSE202203', destdir=".",getGPL =T)
show(gset)
gset=gset[[1]]
#pdata
pdata=pData(gset)
View(pdata)
table(pdata$source_name_ch1)
pdata=pdata[,c(1,26,27,29,30)]
#express
ex=fread('GSE202203_TPM_Raw_gene_3207.tsv.gz')
a=as.data.frame(b)
ex2=merge(a,ex,by.x=1,by.y=1)
rownames(ex2)=ex2[,1]
ex2=ex2[,-c(1)]
min(ex2)
ex2=log(ex2+0.1)
ex3=as.data.frame(t(ex2))
test2=merge(pdata,ex3,by.x=1,by.y=0)
rownames(test2)=test2[,1]
test2=test2[,-1]
test2os=test2[,-(3:4)]
test2rfs=test2[,-(1:2)]
colnames(test2os)[1:2]=c('fultime','fulstate')
colnames(test2rfs)[1:2]=c('rfstime','rfs')
#os change form
a1=strsplit(test2os$fultime,split = ": ",fixed = T) 
t1 = sapply(a1,function(x){x[2]})
test2os$fultime=t1
b1=strsplit(test2os$fulstate,split = ": ",fixed = T) 
s1= sapply(b1,function(x){x[2]})
test2os$fulstate=s1
test2os=subset(test2os,fulstate!='NA')
#rfs
a2=strsplit(test2rfs$rfstime,split = ": ",fixed = T) 
t2 = sapply(a2,function(x){x[2]})
test2rfs$rfstime=t2
b2=strsplit(test2rfs$rfs,split = ": ",fixed = T) 
s2= sapply(b2,function(x){x[2]})
test2rfs$rfs=s2
test2rfs=subset(test2rfs,rfstime!='NA')
#OS validation
riskscore2=predict(cox2,type='risk',newdata = test2os)
risk4=as.vector(ifelse(riskscore2>median(riskscore2),'hign','low'))
write.table(cbind(id=rownames(cbind(test2os[,1:2],riskscore2,risk4)),cbind(test2os[,1:2],riskscore2,risk4)),file="risk.txt",sep="\t",quote=F,row.names=F)
#survive
risk5=read.table("risk.txt",header=T,sep="\t")
diff2=survdiff(Surv(fultime, fulstate) ~risk4,data = risk5)
pValue2=1-pchisq(diff2$chisq,df=1)
pValue2
p2=signif(pValue2,3)
diff2
fit2 <- survfit(Surv(fultime, fulstate) ~ risk4, data = risk5)
summary(fit2)    
pdf(file="16-survival2.pdf")
plot(fit2, lty = 2:3,col=c("red","blue"),xlab="time (year)",ylab="survival rate",
     main=paste("survival curve (p<", 2e-16 ,")",sep=""),mark.time=T,
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
legend("topright", c("high risk", "low risk"), lty = 2:3, col=c("red","blue"))
dev.off()
#roc
risk6=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="16-ROC2.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc2=survivalROC(Stime=risk6$fultime, status=risk6$fulstate, marker = risk6$riskscore2, 
                 predict.time =5, method="KM")
roc2$AUC
plot(roc2$FP, roc2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc2$AUC,3),")"),
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
abline(0,1)
dev.off()
m2=round(median(risk5$riskscore2),3)
#RFS
riskscore3=predict(cox2,type='risk',newdata = test2rfs)
risk7=as.vector(ifelse(riskscore3>median(riskscore3),'hign','low'))
write.table(cbind(id=rownames(cbind(test2rfs[,1:2],riskscore3,risk7)),cbind(test2rfs[,1:2],riskscore3,risk7)),file="risk2.txt",sep="\t",quote=F,row.names=F)
#survive
risk8=read.table("risk2.txt",header=T,sep="\t")
diff3=survdiff(Surv(rfstime, rfs) ~risk7,data = risk8)
pValue3=1-pchisq(diff3$chisq,df=1)
p3=signif(pValue3,3)
diff3
fit3 <- survfit(Surv(rfstime, rfs) ~ risk7, data = risk8)
summary(fit3)    
pdf(file="16-survival3.pdf")
plot(fit3, lty = 2:3,col=c("red","blue"),xlab="time (year)",ylab="survival rate",
     main=paste("survival curve (p=",p3 ,")",sep=""),mark.time=T,
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
legend("topright", c("high risk", "low risk"), lty = 2:3, col=c("red","blue"))
dev.off()
#roc
risk9=read.table("risk2.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="16-ROC3.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc3=survivalROC(Stime=risk9$rfstime, status=risk9$rfs, marker = risk9$riskscore3, 
                 predict.time =5, method="KM")
roc3$AUC
plot(roc3$FP, roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",'0.680',")"),
     lwd = 2, cex.main=1, cex.lab=1, cex.axis=1, font=1)
abline(0,1)
dev.off()
m3=round(median(risk8$riskscore),3)
