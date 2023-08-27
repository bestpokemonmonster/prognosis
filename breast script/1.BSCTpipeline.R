Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
#set work path
getwd()
setwd('/home/bsb/breast/scRNA/GSM5956094/')
set.seed(2023)
#load library
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
#load data
scRNA.counts=Read10X('./')
scRNA=CreateSeuratObject(scRNA.counts,min.cells = 3,project = 'BC.B',min.features = 500)
#mt
scRNA[['percent.mt']]<-PercentageFeatureSet(scRNA,pattern = '^MT-')
#red blood cell
HB.genes<-c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
#View(scRNA@meta.data)
#draw
pdf('../SCTB/1-quality.pdf')
col.num<-length(levels(scRNA@active.ident))
rainbow(col.num)
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #if not show point，set pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
dev.off()
###relation of different parameter
plot1=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
plot1
plot2
plot3
pdf('../SCTB/2-bcscplot1.pdf')
pearplot
dev.off()
#filter
scRNA1 <- subset(scRNA, subset = nFeature_RNA < 5000 & percent.mt < 15 & percent.HB < 1 & nCount_RNA < 40000)
scRNA
scRNA1
expB <- as.data.frame(scRNA1[["RNA"]]@counts)
save(expB,file = '../SCTB/expB.RData')
write.table(expB,'../SCTB/expB.txt',sep = '\t') 
#SCTransform
scRNA1 <- SCTransform(scRNA1, vars.to.regress = "percent.mt", verbose = FALSE)#排除线粒体影响归一化
##results are stored in：
scRNA1@assays$SCT
#or
scRNA1[["SCT"]]
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1)) 
plot4 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident") 
plot4
pdf('../SCTB/3-bcelbowplot.pdf')
ElbowPlot(scRNA1, ndims=50, reduction="pca") 
dev.off()
pc.num=1:33
#clustering
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num)
library(clustree)
objB=scRNA1
objB<-FindClusters(objB,resolution=seq(0.1,1.2,by=0.05))
pdf('../SCTB/4-bcclustree.pdf')
clustree(objB)
dev.off()
###resolution decided according data，while cells more then 20000，resolution usually more than 1.0
scRNA1 <- FindClusters(scRNA1, resolution = 0.2)
#TSNE
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
###label = TRUE means show annotation
pdf('../SCTB/5-bctsne.pdf')
DimPlot(scRNA1, reduction = "tsne",label = TRUE) 
dev.off()
#UMAP
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
pdf('../SCTB/5-bcumap.pdf')
DimPlot(scRNA1, reduction = "umap") 
dev.off()


save(scRNA1,file='../SCTB/scRNA1.RData')
rm(list=ls())
load('../SCTB/scRNA1.RData')
pc.num=1:33
#remove doublet
library(DoubletFinder)
#test best parameter
sweep.res.list <- paramSweep_v3(scRNA1, PCs = pc.num, sct = T)
#log normalization，sct = F（default）,set T if use SCTransform
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #see the best parameter point
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #get best pK value
##
DoubletRate = ncol(scRNA1)*8*1e-6 
#estimate proportion of homo-doublets
homotypic.prop <- modelHomotypic(scRNA1$seurat_clusters) 
#count doublets
nExp_poi <- round(DoubletRate*ncol(scRNA1)) 
#correct doublets 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
##identify doublets
scRNA1<- doubletFinder_v3(scRNA1, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
##Classification results in scRNA_harmony@meta.data
colnames(scRNA1@meta.data)
pdf('../SCTB/6-doublecell.pdf')
DimPlot(scRNA1, reduction = "umap", group.by = "DF.classifications_0.25_0.02_630")
dev.off()
save(scRNA1,file='../SCTB/scRNA2.RData')
View(scRNA1@meta.data)
table(scRNA1@meta.data$DF.classifications_0.25_0.02_630)
scRNA1=subset(scRNA1,DF.classifications_0.25_0.02_630=='Singlet')
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num)
scRNA1 <- FindClusters(scRNA1, resolution = 0.2)
#SingleR annotation
library(celldex)
refdata<-celldex::HumanPrimaryCellAtlasData()
save(refdata,file = '../SCTB/refdata.RData')
#load('../SCTB/refdata.RData')
#refdata@colData@rownames
#unique(refdata$label.main)
testdata <- GetAssayData(scRNA1, slot="data")
###start single R
library(SingleR)
?SingleR
predictions <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main)
table(predictions$labels)
predictions@listData[['scores']][1:5,1:5]
###make annotation file
celltype=data.frame(cell=rownames(scRNA1@meta.data),seurat=scRNA1@meta.data$seurat_clusters,
                    predict=predictions$labels)
rownames(celltype)=celltype$cell
celltype=celltype[,2:3]
sort(table(celltype[,1]))
table(celltype[,1:2])

#most cells in every cluster are named as cell type
lalala <- as.data.frame(table(celltype[,1:2]))
finalmap <- lalala %>% group_by(seurat) %>% top_n(n = 1, wt = Freq)#find most cells in every cluster
finalmap2 <-finalmap[order(finalmap$seurat),]$predict#find cell type from 0:n cluster
print(finalmap2)
testname <- scRNA1
new.cluster.ids <- as.character(finalmap2)
names(new.cluster.ids) <- levels(testname)
testname <- RenameIdents(testname, new.cluster.ids)
p1<- DimPlot(testname,label = F)
p2<-DimPlot(testname,group.by = 'seurat_clusters',label = T)
pdf('../SCTB/7-bcannotation.pdf',width=8,height = 6)
p1
dev.off()
pdf('../SCTB/7-bcannotation2.pdf',width=8,height = 6)
p2
dev.off()
#epithelial cell
Marker1= c('CDH1','EPCAM','ESR1','KRT18',"KRT19",
           'CD2','LCK','CD247','CD96','IL7R')
pdf('../SCTB/8-epithelialdot.pdf',width=10,height = 7)
DotPlot(testname,features=Marker1)
dev.off()

#extract epithelial and T cells，add to metadata
View(scRNA1@meta.data)
scRNA1@meta.data$celltype=celltype$predict
p3<-DimPlot(scRNA1,group.by = 'seurat_clusters',label = T)
p4<-DimPlot(scRNA1,group.by = 'celltype',label = T)
pdf('../SCTB/9-umap.pdf')
p3+p4
dev.off()
scRNA2 <- subset(scRNA1, subset= celltype=='Epithelial_cells' | celltype=='T_cells')
View(scRNA2@meta.data)
save(scRNA2,file = '../SCTB/scRNA3.RData')
save.image('../SCTB/BSCTpipeline.RData')

