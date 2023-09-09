Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd('/home/bsb/breast/bulk/')
set.seed(2023)
library(RobustRankAggreg)
#load data
load('TCGA/cancermarker.RData')
load('GSE109169/cancermarker2.RData')
load('GSE93601/cancermarker3.RData')
load('GSE65194/cancermarker4.RData')
load('TCGA/normalmarker.RData')
load('GSE109169/normalmarker2.RData')
load('GSE93601/normalmarker3.RData')
load('GSE65194/normalmarker4.RData')
#filte 
cancermarker=subset(cancermarker,padj<0.05)
normalmarker=subset(normalmarker,padj<0.05)
tumormarker2=subset(tumormarker2,adj.P.Val<0.05)
normalmarker2=subset(normalmarker2,adj.P.Val<0.05)
tumormarker3=subset(tumormarker3,adj.P.Val<0.05)
normalmarker3=subset(normalmarker3,adj.P.Val<0.05)
tumormarker4=subset(tumormarker4,adj.P.Val<0.05)
normalmarker4=subset(normalmarker4,adj.P.Val<0.05)
#order
library(dplyr)
cancer1=arrange(cancermarker,log2FoldChange)
cancer2=arrange(tumormarker2,desc(logFC))
cancer3=arrange(tumormarker3,desc(logFC))
cancer4=arrange(tumormarker4,desc(logFC))
normal1=arrange(normalmarker,desc(log2FoldChange))
normal2=arrange(normalmarker2,logFC)
normal3=arrange(normalmarker3,logFC)
normal4=arrange(normalmarker4,logFC)
#order
c1=cancer1$row
c2=rownames(cancer2)
c3=rownames(cancer3)
c4=rownames(cancer4)
cancerlist=list(c1,c2,c3,c4)
cancersign=aggregateRanks(glist = cancerlist)
n1=normal1$row
n2=rownames(normal2)
n3=rownames(normal3)
n4=rownames(normal4)
normallist=list(n1,n2,n3,n4)
normalsign=aggregateRanks(glist = normallist)
cancersign=subset(cancersign,Score<0.05)
normalsign=subset(normalsign,Score<0.05)
save(cancersign,file='./RobustRank/cancersign.RData')
save(normalsign,file='./RobustRank/normalsign.RData')

#venn
library(gplots)
venn(cancerlist)
venn(normallist)
library(VennDiagram)
?venn.diagram
pdf('15-breastcancer.pdf')
venn.plot <- venn.diagram(
  x = list(
    TCGA_BRCA= cancer1$row,
    GSE109169=rownames(cancer2),
    GSE93601=rownames(cancer3),
    GSE65194=rownames(cancer4)
  ),
  filename = NULL,
  col = "black",
  lty = "dotted", 
  lwd = 3, 
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("blue", 'blue', "blue",'blue', "blue", "blue",
                        "blue", "blue", "yellow", "blue",
                        "blue", "red", "red", "red", "red"),
                        cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1,
  cat.fontface = "bold",
  cat.fontfamily = "serif"
)
grid.newpage()
grid.draw(venn.plot)
dev.off()
pdf('15-breastnormal.pdf')
venn.plot2 <- venn.diagram(
  x = list(
    TCGA_BRCA= normal1$row,
    GSE109169=rownames(normal2),
    GSE93601=rownames(normal3),
    GSE65194=rownames(normal4)
  ),
  filename = NULL,
  col = "black",
  lty = "dotted", 
  lwd = 3, 
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("blue", 'blue', "blue",'blue', "blue", "blue",
                "blue", "blue", "yellow", "blue",
                "blue", "red", "red", "red", "red"),
                        cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1,
  cat.fontface = "bold",
  cat.fontfamily = "serif"
)
grid.newpage()
grid.draw(venn.plot2)
dev.off()

C100=cancersign[1:100,]
C100$rank=1:100
N100=normalsign[1:100,]
N100$rank=1:100
write.csv(C100,'C100.csv')
write.csv(N100,'N100.csv')
