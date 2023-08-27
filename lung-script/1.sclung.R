Sys.setenv(LANGUAGE="en")
options(stringsAsFactors=FALSE)
rm(list=ls())
setwd('/home/bsb/lung/sclung/')
set.seed(2023)
#load library
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
#read data
P6=read.table('./data/GSM4453581_P6_exp.txt',header=TRUE,row.names=1)
P14=read.table('./data/GSM4453589_P14_exp.txt',header=TRUE,row.names=1)
P18=read.table('./data/GSM4453593_P18_exp.txt',header=TRUE,row.names=1)
P19=read.table('./data/GSM4453594_P19_exp.txt',header=TRUE,row.names=1)
P40=read.table('./data/GSM4453615_P40_exp.txt',header=TRUE,row.names=1)
P41=read.table('./data/GSM4453616_P41_exp.txt',header=TRUE,row.names=1)
P6<-CreateSeuratObject(counts=P6,project="P6",min.cells=3,min.features=200)
P14<-CreateSeuratObject(counts=P14,project="P14",min.cells=3,min.features=200)
P18<-CreateSeuratObject(counts=P18,project="P18",min.cells=3,min.features=200)
P19<-CreateSeuratObject(counts=P19,project="P19",min.cells=3,min.features=200)
P40<-CreateSeuratObject(counts=P40,project="P40",min.cells=3,min.features=200)
P41<-CreateSeuratObject(counts=P41,project="P41",min.cells=3,min.features=200)
#harmony
scRNA_harmony<-merge(P6,y=c(P14,P18,P19,P40,P41))
scRNA_harmony<-NormalizeData(scRNA_harmony)%>%FindVariableFeatures(nfeatures=2000)%>%
              ScaleData()%>%RunPCA(npcs=30)
scRNA_harmony<-RunHarmony(scRNA_harmony,group.by.vars="orig.ident")
scRNA_harmony<-FindNeighbors(scRNA_harmony,reduction="harmony",dims=1:30)%>%FindClusters(resolution=0.5)
scRNA_harmony<-RunUMAP(scRNA_harmony,reduction="harmony",dims=1:30)
plot1=DimPlot(scRNA_harmony,reduction="umap",label=T)
plot2=DimPlot(scRNA_harmony,reduction="umap",group.by='orig.ident')
#combinate
plotc<-plot1+plot2
pdf('1-lung.pdf')
plotc
dev.off()

#annotation
Marker9=c('CAPS','SNTN','EPCAM','CLDN5','VWF','PECAM1',
         'GATA2','TPSAB1','TPSB2','COL1A1','COL1A2','DCN',
         'CD2','CD3D','CD3E','CD3G','CD79A','CD79B',
         'CD14','LYZ','CSF3R','S100A8','S100A9')
pdf('2-lung.pdf',width = 11,height = 9)
DotPlot(scRNA_harmony,features=Marker9)+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()
#'CAPS','SNTN','EPCAM'-Epithelial cells,cancer£»'CLDN5','VWF','PECAM1'-Endothelial cells£»
#'GATA2','TPSAB1','TPSB2'-Mast£»'COL1A1','COL1A2','DCN'-Fibroblasts£»
#'CD2','CD3D','CD3E','CD3G'-T£»'CD79A','CD79B'-B£»
#'CD14','LYZ'-Myeloid cells £»'CSF3R','S100A8','S100A9'- Neutrophils

new.cluster.ids<-c("Neotrophil","Cancer","Cancer","Cancer","Cancer",
                  "Myeloid","Cancer","Epithelial","Neotrophil",
                  "T_cell","Fibroblast","Fibroblast","Neotrophil",
                  "Endothelial","B_cell","B_cell","Mast_cell")
names(new.cluster.ids)<-levels(scRNA_harmony)
scRNA_harmony<-RenameIdents(scRNA_harmony,new.cluster.ids)
pdf('4-lung.pdf')
DimPlot(scRNA_harmony,reduction="umap",label=TRUE,pt.size=0.5)+NoLegend()
dev.off()
#name type
type=c("Neotrophil","Cancer","Cancer","Cancer","Cancer",
       "Myeloid","Cancer","Epithelial","Neotrophil",
       "T_cell","Fibroblast","Fibroblast","Neotrophil",
       "Endothelial","B_cell","B_cell","Mast_cell")
scRNA_harmony@meta.data$celltype=NA
View(scRNA_harmony@meta.data)
for (i in 1:17) {
  scRNA_harmony@meta.data$celltype[which(scRNA_harmony@meta.data$seurat_clusters%in%c(i-1))]<-type[i]
}
table(scRNA_harmony@meta.data$seurat_clusters,scRNA_harmony@meta.data$celltype)
pdf('5-lung.pdf')
DimPlot(scRNA_harmony,reduction="umap",group.by='celltype')
dev.off()


diffgene <- FindMarkers(object = scRNA_harmony, ident.1 = 'Epithelial', ident.2 = 'Cancer',group.by='celltype',min.pct = 0.25,logfc.threshold = 0.5,test.use = 'roc')
?FindMarkers
idx <- !grepl("^MT-", rownames(diffgene))
df_new <- diffgene[idx,]
save(df_new,file='df_new.RData')







