

### Step1-3  Single-Cell Analyze：Monocyte, DC Cell
### 整理时间： 2022/7/22
### 作者： 庞建宇


## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)
gc()

# 设置工作路径
setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step1-3/")



library(Seurat)
load(file = "../Step1-1/res0.4_NSCLC.Integrate.Rdata")
MonDC  = NSCLC.Integrate[,NSCLC.Integrate@active.ident %in% c("Monocyte","Dendritic Cell")] #提取Monocyte和Dendritic Cell簇
rm(NSCLC.Integrate)
gc()


MonDC <- NormalizeData(MonDC)
MonDC <- FindVariableFeatures(MonDC, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(MonDC), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(MonDC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "MonDC-Sub-Top10-VarGene.png",width = 20,height = 10,path = "../../Fig/Step1-3/")


# 数据归一化 + 线形降维
all.genes <- rownames(MonDC)
MonDC <- ScaleData(MonDC, features = all.genes)
# 线性降维 对缩放的数据执行PCA.默认情况下，只使用前面确定的变量特性作为输入，但是如果想选择不同的子集，可以使用features参数来定义。
MonDC <- RunPCA(MonDC, features = VariableFeatures(object = MonDC))
# Examine and visualize PCA results a few different ways
# 查看PCA结果
print(MonDC[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(MonDC, dims = 1:2, reduction = "pca")
ggsave(filename = "MonDC-Sub-PCA.png",width = 16,height = 10,path = "../../Fig/Step1-3/")

DimPlot(MonDC, reduction = "pca", raster=FALSE)
ggsave(filename = "MonDC-Sub-PCA2.png",width = 16,height = 10,path = "../../Fig/Step1-3/") 

DimHeatmap(MonDC, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
ggsave(filename = "MonDC-Sub-PC1_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-3/")

DimHeatmap(MonDC, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
ggsave(filename = "MonDC-Sub-PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-3/")

save(MonDC,file = "res0.4_MonDC.Rdata")


load(file = "res0.4_MonDC.Rdata")
MonDC <- JackStraw(MonDC, num.replicate = 100) # 34m 60s
MonDC <- ScoreJackStraw(MonDC, dims = 1:20)
JackStrawPlot(MonDC, dims = 1:20)
ggsave(filename = "res0.4-MonDC-Sub-JackPlot.png",width = 16,height = 10,path = "../../Fig/Step1-3/")

ElbowPlot(MonDC)#肘部图 
# 综合以上方法，选择13个主成成分作为参数用于后续分析。
ggsave(filename = "res0.4-MonDC-Sub-ElbowPlot.png",width = 12,height = 10,path = "../../Fig/Step1-3/")


#细胞聚类 KNN算法
library(clustree)
MonDC <- FindNeighbors(MonDC, dims = 1:13)#dims = 1:13 即选取前13个主成分来分类细胞。
MonDC <- FindClusters(object = MonDC,
                      resolution = c(seq(0,1,by = 0.1)))
clustree(MonDC@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "MonDC-Sub-resolution(0-1).png",width = 20,height = 14,path = "../../Fig/Step1-3/")


#选取resolution = 0.3 作为后续分析参数
# Assign identity of clusters
Idents(object = MonDC) <- "RNA_snn_res.0.3"
MonDC@meta.data$seurat_clusters = MonDC@meta.data$RNA_snn_res.0.3
head(Idents(MonDC), 5)#查看前5个细胞的分类ID


# 非线性降维 UMAP/TSNE
# UMAP
MonDC <- RunUMAP(MonDC, dims = 1:13)
DimPlot(MonDC, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "MonDC-Sub-UMAP-label.png",width = 18,height = 12,path = "../../Fig/Step1-3/")


# 找每个簇的差异基因
library(dplyr)
MonDCMarker <- FindAllMarkers(MonDC, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <-MonDCMarker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(MonDC, features = top10$gene) + NoLegend()
ggsave(filename = "MonDC-Sub-Top10-MarkerGene.png",width = 16,height = 10,path = "../../Fig/Step1-3/")

save(top10,MonDCMarker,file = "res0.4_MonDC-Sub-Markers.Rdata")
save(MonDC, file = "res0.4_MonDC.Rdata")#保存Rdata，用于后续分析



load(file = "res0.4_MonDC-Sub-Markers.Rdata")
load(file = "res0.4_MonDC.Rdata")
# 自动化注释  SingleR
sce = MonDC
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.3
clusters=sce@meta.data$seurat_clusters

Blue.ref <- celldex::BlueprintEncodeData()
pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Mona.ref <- celldex::MonacoImmuneData()
pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

Nover.ref <- celldex::NovershternHematopoieticData()
pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")



cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    Blue=pred.Blue.ref$labels,
                    DICE=pred.DICE.ref$labels,
                    HPCA=pred.HPCA.ref$labels,
                    Mona=pred.Mona.ref$labels,
                    Nover=pred.Nover.ref$labels )

head(cellType)
sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']


pro='MonDC-Sub-SingleR_anno_res0.3'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 18,height = 14,path = "../../Fig/Step1-3/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 18,height = 14,path = "../../Fig/Step1-3/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 18,height = 14,path = "../../Fig/Step1-3/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 18,height = 14,path = "../../Fig/Step1-3/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 18,height = 14,path = "../../Fig/Step1-3/")

MonDC = sce

save(top10,MonDCMarker,cellType,file = "res0.4_MonDC-Sub-Markers.Rdata")
save(MonDC, file = "res0.4_MonDC.Rdata")#保存Rdata，用于后续分析



write.csv(cellType,file = "MonDC-Sub-cellType.csv")
write.csv(MonDCMarker, file =  "MonDC-Sub-MyeloidCellMarker.csv")
write.csv(top10, file = "MonDC-Sub-top10gene.csv")


#C0,C2,C5,C6,C7,C9,C11：Monocyte       Marker:"FBP1" , "MARCO" ,"RETN"
#C1,C12：Macrophage                    Marker:"FOLR2"   ,"PLTP" , "SEPP1" ,"CCL18"
#C3,C4：Dendritic Cell                 Marker:"CD48" , "FCN1"  , "IL1B" ,"LST1"
#C8：Granulocyte-Monocyte Progenitor   Marker:"KIAA0101" , "UBE2C"  , "PTTG1"
#10：Plasmacytoid Dendritic Cell       Marker:"GZMB" , "JCHAIN"  , "ITM2C","C12orf75"
#C13：Granulosa cell                   Marker:"BIRC3" , "CCL17"  , "CCL19"
#C14：Neutrophil                       Marker:"PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R"

load( file = "res0.4_MonDC.Rdata")#保存Rdata，用于后续分析
library(dplyr)
MonMarker <- FindMarkers(MonDC, ident.1 = c(0,2,5,6,7,9,11), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
MonMarker_top10 <-MonMarker  %>% top_n(n = 10, wt = avg_log2FC)

MacMarker <- FindMarkers(MonDC, ident.1 = c(1,12), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
MacMarker_top10 <-MacMarker  %>% top_n(n = 10, wt = avg_log2FC)

DCMarker <- FindMarkers(MonDC, ident.1 = c(3,4), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
DCMarker_top10 <-DCMarker  %>% top_n(n = 10, wt = avg_log2FC)

new.cluster.ids <- c("Monocyte", "Macrophage", "Monocyte", "Dendritic Cell", "Dendritic Cell", "Monocyte", "Monocyte", 
                     "Monocyte", "Granulocyte-Monocyte Progenitor", "Monocyte","Plasmacytoid Dendritic Cell","Monocyte",
                     "Macrophage","Granulosa cell","Neutrophil")

names(new.cluster.ids) <- levels(MonDC)
MonDC <- RenameIdents(MonDC, new.cluster.ids)
DimPlot(MonDC, reduction = "umap", label = TRUE, pt.size = 0.5,raster=FALSE)
ggsave(filename = "MonDC-Sub-UMAP_Cellmarker.png",width = 14,height = 10,path = "../../Fig/Step1-3/")

# VlnPlot(MonDC, features = "PLIN2")
# ggsave(filename = "Dendritic Cell_CNP.png",width = 18,height = 12,path = "../2*GSE_fig")
# 
# FeaturePlot(MonDC, features = "CCL17")
# ggsave(filename = "Dendritic Cell_CNP(1).png",width = 18,height = 12,path = "../2*GSE_fig")

genes_to_check = c("FBP1" , "MARCO" ,"RETN","FOLR2"   ,"PLTP" , "SEPP1" ,"CCL18","CD48" , "FCN1"  , "IL1B" ,"LST1",
                   "KIAA0101" , "UBE2C"  , "PTTG1","GZMB" , "JCHAIN"  , "ITM2C","C12orf75","BIRC3" , "CCL17"  , "CCL19",
                   "PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R")
DotPlot(MonDC,features = unique(genes_to_check)) + RotatedAxis()
ggsave(filename = "MonDC-Sub-All Marker.png",width = 18,height = 12,path = "../../Fig/Step1-3/")

save(MonDC, file = "res0.4_MonDC.Rdata")#保存Rdata，用于后续分析


# 2022.6.20  美化图片
load(file = 'res0.4_MonDC.Rdata')
library(Seurat)
library(ggsci)
library(ggplot2)
pal = pal_ucscgb(alpha = 0.7)(7)
table(MonDC@active.ident)
DimPlot(MonDC, reduction = "umap", label = T, label.box = T, cols = pal, raster = T) + NoLegend()
ggsave(filename = 'UMAP_MonDC_7type.pdf',width = 8,height = 6,path = '../../Fig/Step1-3/')

MonDC@meta.data[["Rename"]] = MonDC@active.ident
Idents(object = MonDC) <- "RNA_snn_res.0.3"
table(MonDC@active.ident)
pal2 = pal_ucscgb(alpha = 0.7)(15)
DimPlot(MonDC, reduction = "umap", label = T, label.box = T, cols = pal2, raster = T) + NoLegend()
ggsave(filename = 'UMAP_MonDC_res0.3.pdf',width = 8,height = 6,path = '../../Fig/Step1-3/')


new.cluster.ids <- c("Monocyte", "Macrophage", "Monocyte", "Dendritic Cell", "Dendritic Cell", "Monocyte", "Monocyte", 
                     "Monocyte", "GMP", "Monocyte","Plasmacytoid Dendritic Cell","Monocyte",
                     "Macrophage","Granulosa cell","Neutrophil")
names(new.cluster.ids) <- levels(MonDC)
MonDC <- RenameIdents(MonDC, new.cluster.ids)
table(MonDC@active.ident)
MonDC@meta.data[["Rename2"]] = MonDC@active.ident
DimPlot(MonDC, reduction = "umap", label = T, label.box = T, cols = pal2, raster = T) + NoLegend()
ggsave(filename = 'UMAP_MonDC_7type2.pdf',width = 8,height = 6,path = '../../Fig/Step1-3/')


Idents(object = MonDC) <- "Rename"
table(MonDC@active.ident)

save(MonDC,file = 'res0.4_MonDC.Rdata')

#C0,C2,C5,C6,C7,C9,C11：Monocyte       Marker:"FBP1" , "MARCO" ,"RETN"
#C1,C12：Macrophage                    Marker:"FOLR2"   ,"PLTP" , "SEPP1" ,"CCL18"
#C3,C4：Dendritic Cell                 Marker:"CD48" , "FCN1"  , "IL1B" ,"LST1"
#C8：Granulocyte-Monocyte Progenitor   Marker:"KIAA0101" , "UBE2C"  , "PTTG1"
#10：Plasmacytoid Dendritic Cell       Marker:"GZMB" , "JCHAIN"  , "ITM2C","C12orf75"
#C13：Granulosa cell                   Marker:"BIRC3" , "CCL17"  , "CCL19"
#C14：Neutrophil                       Marker:"PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R"
genes_to_check = c("FBP1" , "MARCO" ,"RETN","FOLR2"   ,"PLTP" , "SEPP1" ,"CCL18","CD48" , "FCN1"  ,
                   "KIAA0101" , "UBE2C"  , "PTTG1","GZMB" , "JCHAIN"  , "ITM2C","C12orf75","BIRC3" , "CCL17"  , "CCL19",
                   "PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R")
DotPlot(MonDC,features = unique(genes_to_check), cols = c("lightgrey", "orange")) + RotatedAxis()
ggsave(filename = "All Marker_MonDC.pdf",width = 12,height = 8,path = "../../Fig/Step1-3/")

# 'MARCO','SEPP1','CD48','KIAA0101','GZMB','BIRC3','PI3'
# genes_to_check = c('MARCO','SEPP1','CD48','KIAA0101','GZMB','BIRC3','PI3')
# FeaturePlot(MonDC, features = unique(genes_to_check), cols = c("lightgrey", "orange"),ncol = 4)
# ggsave(filename = "All Marker_MonDC2.png",width = 18,height = 9,path = "../2*GSE_fig/res0.4/Marker")



