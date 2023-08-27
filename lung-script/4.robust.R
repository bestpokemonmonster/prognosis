Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd('/home/bsb/lung/robust/')
set.seed(2023)
library(RobustRankAggreg)
#read
load('../TCGAlung/cancermarker.RData')
load('../TCGAlung/normalmarker.RData')
load('../GSE151103/normalmarker2.RData')
load('../GSE151103/cancermarker2.RData')
load('../GSE81089/cancermarker3.RData')
load('../GSE81089/normalmarker3.RData')
load('../GSE30219/normalmarker4.RData')
load('../GSE30219/cancermarker4.RData')
#order
library(dplyr)
cancer1=arrange(cancermarker,log2FoldChange)
cancer2=arrange(tumormarker2,desc(logFC))
cancer3=arrange(cancermarker3,log2FoldChange)
cancer4=arrange(tumormarker4,desc(logFC))
normal1=arrange(normalmarker,desc(log2FoldChange))
normal2=arrange(normalmarker2,logFC)
normal3=arrange(normalmarker3,desc(log2FoldChange))
normal4=arrange(normalmarker4,logFC)


#filte
cancer1=subset(cancer1,padj<0.05)
normal1=subset(normal1,padj<0.05)
cancer2=subset(cancer2,adj.P.Val<0.05)
normal2=subset(normal2,adj.P.Val<0.05)
cancer3=subset(cancer3,padj<0.05)
normal3=subset(normal3,padj<0.05)
cancer4=subset(cancer4,adj.P.Val<0.05)
normal4=subset(normal4,adj.P.Val<0.05)
#
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
save(cancersign,file='./cancersign.RData')
save(normalsign,file='./normalsign.RData')

#venn
library(gplots)
venn(cancerlist)
venn(normallist)
library(VennDiagram)
pdf('11-lungcancer.pdf')
venn.plot <- venn.diagram(
  x = list(
    TCGA=cancer1$row,
    GSE151103=rownames(cancer2),
    GSE81089=rownames(cancer3)
  ),
  filename = NULL,
  col = "black",
  lty = "dotted", 
  lwd = 3, 
  fill = c("cornflowerblue", "green", "yellow"),
  alpha = 0.50,
  label.col = c("orange", "red", "darkorchid4", "black",
                "white", "darkblue", "blue"
               ),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "serif"
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

pdf('11-lungnormal.pdf')
venn.plot <- venn.diagram(
  x = list(
    TCGA=normal1$row,
    GSE151103=rownames(normal2),
    GSE81089=rownames(normal3)
  ),
  filename = NULL,
  col = "black",
  lty = "dotted", 
  lwd = 3, 
  fill = c("cornflowerblue", "green", "yellow"),
  alpha = 0.50,
  label.col = c("orange", "red", "darkorchid4", "black",
                "white", "darkblue", "blue"
  ),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "serif"
)
grid.newpage()
grid.draw(venn.plot)
dev.off()