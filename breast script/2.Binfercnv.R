Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("/home/bsb/breast/scRNA/SCTB/")
set.seed(2023)
load('scRNA3.RData')
#load library
library(Seurat)
library(infercnv)
library(dplyr)
#re-cluster
scRNA2 <- SCTransform(scRNA2, vars.to.regress = "percent.mt", verbose = FALSE)
scRNA2 <- RunPCA(scRNA2, features = VariableFeatures(scRNA2))
pdf('10-bcsubsetelbowplot.pdf')
ElbowPlot(scRNA2, ndims=50, reduction="pca") 
dev.off()
pc.num=1:20
scRNA2 <- FindNeighbors(scRNA2, dims = pc.num)
library(clustree)
objB=scRNA2
objB<-FindClusters(objB,resolution=seq(0.1,1.2,by=0.1))
pdf('11-bcsubsetclustree.pdf',width = 9,height =9 )
clustree(objB)
dev.off()
scRNA2 <- FindClusters(scRNA2, resolution = 1.0)
scRNA2 = RunUMAP(scRNA2, dims = pc.num)
View(scRNA2@meta.data)
a<-DimPlot(scRNA2, reduction = "umap",group.by='celltype',label = TRUE) 
b<-DimPlot(scRNA2,reduction = 'umap',group.by = 'seurat_clusters',label = F)
pdf('12-Tandepi1.pdf')
a
dev.off()
pdf('12-Tandepi2.pdf')
b
dev.off()

#marker
Marker1= c('CDH1','EPCAM','ESR1','KRT18',"KRT19",
           'CD2','LCK','CD247','CD96','IL7R')
pdf('13-epidot.pdf',width=10,height = 7)
DotPlot(scRNA2,features=Marker1)
dev.off()
VlnPlot(scRNA2,features = Marker1,pt.size = 0,ncol = 10)
FeaturePlot(scRNA2,reduction='umap', features = Marker1 )


pos=read.table("./hg38_gencode_v27.txt")
pos=subset(pos,V2!='chrY' & V2!='chrM')
table(pos$V2)
#pos1=read.table('human.gene(1)(2).positions')
#table(pos1$V2)
#delete repeat gene
#pos1=distinct(pos,V7,.keep_all = TRUE)
#position reorder
#rownames(pos1)=pos1$V7
#pos2=select(pos1,V7,V2,V3,V4)
#View(pos2)
write.table(pos, 'geneLocate.txt', row.names=F, col.names=F, sep='\t')
#extract express-info
exprMatrix <- as.matrix(GetAssayData(scRNA2, slot='counts'))
cellAnnota <- subset(scRNA2@meta.data, select='seurat_clusters')
groupFiles='groupFiles.txt'
dim(exprMatrix)
write.table(cellAnnota,file =" groupFiles.txt",sep = '\t',col.names = F)
#create infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file=" groupFiles.txt",
                                    delim="\t",
                                    gene_order_file= "geneLocate.txt",
                                    ref_group_names=c('0','3','11'))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir=  'cnv1/' ,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   #  cluster_by_groups£ºidentify sourse,then cluster
                             hclust_method="ward.D2", plot_steps=F,
                             HMM=FALSE,
                             denoise = TRUE)
##identify malignant cell
grp=read.table("cnv1/infercnv.observation_groupings.txt",sep = "",header = T)
obs=read.table("cnv1/infercnv.observations.txt",header = T,check.names =F)
max(obs)
min(obs)
obs[obs>0.8 & obs<0.9]=2
obs[obs>=0.9 & obs<0.93]=1
obs[obs>=0.93 & obs<1.07]=0
obs[obs>=1.07 & obs<1.1]=1
obs[obs>=1.1 & obs<1.2]=2
scores=as.data.frame(colSums(obs))
scores$cluster=grp$Annotation.Group
colnames(scores)=c("score","cluster")
library(ggpubr)
pdf('14-ggbox.pdf')
ggboxplot(scores,"cluster","score",fill = "cluster")
dev.off()
#find relations between scores cluster and scRNA1 cluster
View(scRNA2@meta.data)
table(scRNA2@meta.data$seurat_clusters) 
table(scores$cluster) #4 and 3,corresponding 1 and 4

#?FindMarkers
cluster1 <- subset(scRNA2, idents = 1)
cluster4 <- subset(scRNA2, idents = 4)
diffgene <- FindMarkers(object = scRNA2, ident.1 = 1, ident.2 = 4,min.pct = 0.25,logfc.threshold = 0.5,test.use = 'roc')
#remove mt
idx <- !grepl("^MT-", rownames(diffgene))
df_new <- diffgene[idx,]
normal=subset(df_new,avg_log2FC>0)
cancer=subset(df_new,avg_log2FC<0)
save(df_new,file='df_new.RData')
save.image('./infercnv.RData')
save(cancer,file = 'cancer.RData')
save(normal,file="normal.RData")
#log>0 non-maglinant£¬log<0 maglinant
nor=AverageExpression(object = cluster1,features = 'CADM2')
nor
can=AverageExpression(object = cluster4,features = 'CADM2')
can


#another try
#low=subset(scores,score<1000)
#low=rownames(low)
#high=subset(scores,score>3000)
#high=rownames(high)
#scRNA2@meta.data$name=rownames(scRNA2@meta.data)
#View(scRNA2@meta.data)
#scRNA2@meta.data$diff=ifelse(scRNA2@meta.data$name %in% high, "high",
#                             ifelse( scRNA2@meta.data$name%in%low,'low','others') )
#table(scRNA2@meta.data$diff)
#marker_genes <- FindMarkers(object=scRNA2,ident.1 ="low",ident.2 ="high",by.group=diff,min.pct=0.1,logfc.threshold=0.25)
#Idents(scRNA2) = "diff"
#marker_genes <- FindMarkers(object = scRNA2, ident.1 = low, ident.2 = high, min.pct = 0.25,logfc.threshold = 0.5,test.use = 'roc')
#idx <- !grepl("^MT-", rownames(marker_genes))
#df_new <- marker_genes[idx,]
#save(df_new,file='df_new.RData')
