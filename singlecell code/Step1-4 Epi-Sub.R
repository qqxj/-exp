

### Step1-4  Single-Cell Analyze：Epithelial/Cancer Cell
### 整理时间： 2022/7/22
### 作者： 庞建宇


## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)
gc()

setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step1-4/")



library(Seurat)
#  提取感兴趣的细胞簇进行亚聚类  # C6,C7,C8,C9,C12,C17,C23,C24,C26 : Epithelial/Cancer Cell 
load(file = "../Step1-1/res0.4_NSCLC.Integrate.Rdata")
Epi  = NSCLC.Integrate[,NSCLC.Integrate@active.ident %in% c("Epithelial/Cancer Cell")] #提取Epithelial/Cancer Cell簇
rm(NSCLC.Integrate)
gc()
Epi <- NormalizeData(Epi)
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(Epi), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(Epi)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "Epi-Sub-Top10-VarGene.png",width = 20,height = 10,path = "../../Fig/Step1-4/")


# 数据归一化 + 线形降维
all.genes <- rownames(Epi)
Epi <- ScaleData(Epi, features = all.genes)
# 线性降维 对缩放的数据执行PCA.默认情况下，只使用前面确定的变量特性作为输入，但是如果想选择不同的子集，可以使用features参数来定义。
Epi <- RunPCA(Epi, features = VariableFeatures(object = Epi))
# Examine and visualize PCA results a few different ways
# 查看PCA结果
print(Epi[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Epi, dims = 1:2, reduction = "pca")
ggsave(filename = "Epi-Sub-PCA.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

DimPlot(Epi, reduction = "pca", raster=FALSE)
ggsave(filename = "Epi-Sub-PCA2.png",width = 16,height = 10,path = "../../Fig/Step1-4/") 

DimHeatmap(Epi, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
ggsave(filename = "Epi-Sub-PC1_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

DimHeatmap(Epi, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
ggsave(filename = "Epi-Sub-PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

save(Epi,file = "res0.4_Epi.Rdata")



load(file = "res0.4_Epi.Rdata")

Epi <- JackStraw(Epi, num.replicate = 100) # 24m 18s
Epi <- ScoreJackStraw(Epi, dims = 1:20)
JackStrawPlot(Epi, dims = 1:20)
ggsave(filename = "res0.4-Epi-Sub-JackPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")
ElbowPlot(Epi)#肘部图 
# 综合以上方法，选择13个主成成分作为参数用于后续分析。
ggsave(filename = "res0.4-Epi-Sub-ElbowPlot.png",width = 12,height = 10,path = "../../Fig/Step1-4/")


#细胞聚类 KNN算法
library(clustree)
Epi <- FindNeighbors(Epi, dims = 1:13)#dims = 1:13 即选取前13个主成分来分类细胞。
Epi <- FindClusters(object = Epi,
                    resolution = c(seq(0,1,by = 0.1)))
clustree(Epi@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "Epi-Sub-resolution(0-1).png",width = 20,height = 14,path = "../../Fig/Step1-4/")


#选取resolution = 0.2 作为后续分析参数
# Assign identity of clusters
Idents(object = Epi) <- "RNA_snn_res.0.2"
Epi@meta.data$seurat_clusters = Epi@meta.data$RNA_snn_res.0.2
head(Idents(Epi), 5)#查看前5个细胞的分类ID


# 非线性降维 UMAP/TSNE
# UMAP
Epi <- RunUMAP(Epi, dims = 1:13)
DimPlot(Epi, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "Epi-Sub-UMAP-label.png",width = 18,height = 12,path = "../../Fig/Step1-4/")


# 找每个簇的差异基因
library(dplyr)
EpiMarker <- FindAllMarkers(Epi, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <-EpiMarker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Epi, features = top10$gene) + NoLegend()
ggsave(filename = "Epi-Sub-Top10-MarkerGene.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

save(top10,EpiMarker,file = "Epi-Sub-Markers.Rdata")
save(Epi, file = "res0.4_Epi.Rdata")#保存Rdata，用于后续分析



load(file = "Epi-Sub-Markers.Rdata")
load(file = "res0.4_Epi.Rdata")

# 自动化注释  SingleR
sce = Epi
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.2
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


pro='Epi-Sub-SingleR_anno_res0.2'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

Epi = sce

save(top10,EpiMarker,cellType,file = "Epi-Sub-Markers.Rdata")
save(Epi, file = "res0.4_Epi.Rdata")#保存Rdata，用于后续分析

write.csv(cellType,file = "Epi-Sub-cellType.csv")
write.csv(EpiMarker, file =  "Epi-Sub-Marker.csv")
write.csv(top10, file = "Epi-Sub-top10gene.csv")




# CopyKAT
load(file = "res0.4_Epi.Rdata")
scRNA <- Epi
counts <- as.matrix(scRNA@assays$RNA@counts)

##运行时间较长  Time difference of 2.45014 days
cnv <- copykat(rawmat=counts, ngene.chr=5, sam.name="NSCLC", n.cores=8)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
saveRDS(cnv, "cnv.rds")

CNV <- readRDS(file = "cnv.rds")
# NSCLC_copykat_prediction.txt是上一步生成的结果
mallignant <- read.delim("/NSCLC_copykat_prediction.txt")
mall <- mallignant
rownames(mall) <- mall[,1]
length(unique(mall$cell.names)) # 有重复
mall <- mall[!duplicated(mall$cell.names),] # 去除重复
rownames(mall) <- mall[,1]
names(mall)
b <- data.frame('copykat.pred' = mall[,-1])
rownames(b) <- rownames(mall)
table(mall$copykat.pred)

# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(scRNA, metadata = b)
table(scRNA@meta.data[["copykat.pred"]])
# p1 <- DimPlot(scRNA, group.by = "RNA_snn_res.0.2", label = T)
# # p2 <- DimPlot(scRNA, group.by = "copykat.pred") + scale_color_manual(values = c("red", "gray"))
# p2 <- DimPlot(scRNA, group.by = "copykat.pred")
# pc <- p1 + p2
# print(pc)
# ggsave("pred_mallignant.png", pc, width = 12, height = 8,path = "./Fig/Step1-4/")

Epi <- scRNA
p1 <- DimPlot(Epi, group.by = "RNA_snn_res.0.2", label = T)
p2 <-DimPlot(Epi, group.by = "copykat.pred")
pc <- p1 + p2
ggsave("pred_mallignant.png", pc, width = 14, height = 8,path = "../../Fig/Step1-4/")

save(Epi,file = "res0.4_Epi.Rdata")

table(Epi@meta.data[["copykat.pred"]])
CancerCell  = Epi[,Epi@meta.data$copykat.pred %in% "aneuploid"] #提取Cancer Cell簇
Normal_Epi  = Epi[,Epi@meta.data$copykat.pred %in% "diploid"] #提取Normal_Epi簇

save(CancerCell,file = "res0.4_CancerCell.Rdata")
save(Normal_Epi,file = "res0.4_Normal_Epi.Rdata")



load(file = "res0.4_Normal_Epi.Rdata")
library(ggplot2)

Normal_Epi <- NormalizeData(Normal_Epi)
Normal_Epi <- FindVariableFeatures(Normal_Epi, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(Normal_Epi), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(Normal_Epi)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "Normal_Epi-Sub-Top10-VarGene.png",width = 20,height = 10,path = "../../Fig/Step1-4/")


# 数据归一化 + 线形降维
all.genes <- rownames(Normal_Epi)
Normal_Epi <- ScaleData(Normal_Epi, features = all.genes)
# 线性降维 对缩放的数据执行PCA.默认情况下，只使用前面确定的变量特性作为输入，但是如果想选择不同的子集，可以使用features参数来定义。
Normal_Epi <- RunPCA(Normal_Epi, features = VariableFeatures(object = Normal_Epi))
# Examine and visualize PCA results a few different ways
# 查看PCA结果
print(Normal_Epi[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Normal_Epi, dims = 1:2, reduction = "pca")
ggsave(filename = "Normal_Epi-Sub-PCA.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

DimPlot(Normal_Epi, reduction = "pca", raster=FALSE)
ggsave(filename = "Normal_Epi-Sub-PCA2.png",width = 16,height = 10,path = "../../Fig/Step1-4/") 

DimHeatmap(Normal_Epi, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
ggsave(filename = "Normal_Epi-Sub-PC1_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

DimHeatmap(Normal_Epi, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
ggsave(filename = "Normal_Epi-Sub-PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

save(Normal_Epi,file = "res0.4_Normal_Epi.Rdata")



load(file = "res0.4_Normal_Epi.Rdata")
Normal_Epi <- JackStraw(Normal_Epi, num.replicate = 100) # 18m 14s
Normal_Epi <- ScoreJackStraw(Normal_Epi, dims = 1:20)
JackStrawPlot(Normal_Epi, dims = 1:20)
ggsave(filename = "res0.4-Normal_Epi-Sub-JackPlot.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

ElbowPlot(Normal_Epi)#肘部图 
# 综合以上方法，选择12个主成成分作为参数用于后续分析。
ggsave(filename = "res0.4-Normal_Epi-Sub-ElbowPlot.png",width = 12,height = 10,path = "../../Fig/Step1-4/")


#细胞聚类 KNN算法
library(clustree)
Normal_Epi <- FindNeighbors(Normal_Epi, dims = 1:12)#dims = 1:12 即选取前12个主成分来分类细胞。
Normal_Epi <- FindClusters(object = Normal_Epi,
                           resolution = c(seq(0,1,by = 0.1)))
clustree(Normal_Epi@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "Normal_Epi-Sub-resolution(0-1).png",width = 20,height = 14,path = "../../Fig/Step1-4/")


#选取resolution = 0.2 作为后续分析参数
# Assign identity of clusters
Idents(object = Normal_Epi) <- "RNA_snn_res.0.2"
Normal_Epi@meta.data$seurat_clusters = Normal_Epi@meta.data$RNA_snn_res.0.2
head(Idents(Normal_Epi), 5)#查看前5个细胞的分类ID


# 非线性降维 UMAP/TSNE
# UMAP
Normal_Epi <- RunUMAP(Normal_Epi, dims = 1:12)
DimPlot(Normal_Epi, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "Normal_Epi-Sub-UMAP-label.png",width = 18,height = 12,path = "../../Fig/Step1-4/")


# 找每个簇的差异基因
library(dplyr)
Normal_EpiMarker <- FindAllMarkers(Normal_Epi, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <-Normal_EpiMarker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Normal_Epi, features = top10$gene) + NoLegend()
ggsave(filename = "Normal_Epi-Sub-Top10-MarkerGene.png",width = 16,height = 10,path = "../../Fig/Step1-4/")

save(top10,Normal_EpiMarker,file = "res0.4_Normal_Epi-Sub-Markers.Rdata")
save(Normal_Epi,file = "res0.4_Normal_Epi.Rdata")#保存Rdata，用于后续分析



load(file = "res0.4_Normal_Epi-Sub-Markers.Rdata")
load(file = "res0.4_Normal_Epi.Rdata")

# 自动化注释  SingleR
sce = Normal_Epi
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.2
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


pro='Normal_Epi-Sub-SingleR_anno_res0.2'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 18,height = 14,path = "../../Fig/Step1-4/")

Normal_Epi = sce

save(top10,Normal_EpiMarker,cellType,file = "res0.4_Normal_Epi-Sub-Markers.Rdata")
save(Normal_Epi,file = "res0.4_Normal_Epi.Rdata")#保存Rdata，用于后续分析

write.csv(cellType,file = "Normal_Epi-Sub-cellType.csv")
write.csv(Normal_EpiMarker, file =  "Normal_Epi-Sub-MyeloidCellMarker.csv")
write.csv(top10, file = "Normal_Epi-Sub-top10gene.csv")



load(file = "res0.4_Normal_Epi.Rdata")
load(file = "res0.4_Normal_Epi-Sub-Markers.Rdata")
library(dplyr)
library(Seurat)
Idents(object = Normal_Epi) <- "RNA_snn_res.0.2"
for (i in 0:18) {
  dat <-  top10 %>% filter(cluster == i)
  assign(paste0('cluster',i),dat)
  
}
# C0    Marker:"KRT17","GPX2","NTRK2",C3     "ALDH3A1","NDUFA4L2",C0
c0gene = cluster0[,7]
genes_to_check = c("KRT17","GPX2","ALDH3A1","NDUFA4L2","NTRK2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C1    Marker:"SFTPC","SFTPA1","SFTPA2","PGCC","SFTPD",C4,C12       "PLA2G1B","PEBP4" C1
c1gene = cluster1[,7]
print(c1gene)
genes_to_check = c("SFTPC","SFTPA1","SFTPA2","SFTPD","PGC","PEBP4","PLA2G1B")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C2    Marker: "C15orf48",C7  "LCN2",C10  "CDKN2A",C2
c2gene = cluster2[,7]
print(c2gene)
genes_to_check = c( "C15orf48","CDKN2A","LCN2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C3    Marker: "AKR1C1","CSTA","AKR1C2",C0      "NTS","NIPSNAP2","NSD3","IGFL2-AS1",C3
c3gene = cluster3[,7]
print(c3gene)
genes_to_check = c("NTS","AKR1C1","NIPSNAP2","NSD3","IGFL2-AS1","CSTA","AKR1C2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C4    Marker: "TFPI2",C5     "GNG11",C4
c4gene = cluster4[,7]
print(c4gene)
genes_to_check = c("TFPI2","GNG11")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C5    Marker: "TFPI2",C4     "HOXB-AS3","NMU",C5
c5gene = cluster5[,7]
print(c5gene)
genes_to_check = c("HOXB-AS3","NMU","TFPI2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C6    Marker:"CST6","NUPR2","PAEP","TCN1",C6
c6gene = cluster6[,7]
print(c6gene)
genes_to_check = c("CST6","NUPR2","PAEP","TCN1")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C7    Marker: "ASS1","BIRC3","TNNC2",C7
c7gene = cluster7[,7]
print(c7gene)
genes_to_check = c("ASS1","BIRC3","TNNC2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C8    Marker: "RACK1","ATP5F1E",C0,C3,C8,C10      "XAF1",C8
c8gene = cluster8[,7]
print(c8gene)
genes_to_check = c("RACK1","ATP5F1E","XAF1")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C9    Marker: "CYP4B1",C16   "MYL9",C15  "CLIC3",C1   "AGER",C9
c9gene = cluster9[,7]
print(c9gene)
genes_to_check = c("AGER","CYP4B1","MYL9","CLIC3")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C10    Marker: "SELENOW","ATP5F1E","SLK",C0,C3,C8      "LCN2",C10
c10gene = cluster10[,7]
print(c10gene)
genes_to_check = c("SELENOW","ATP5F1E","SLK","LCN2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C11    Marker: "MUC5B",C18   "CHI3L1",C3     "HP","CTSE","LINC00922","LTF"C11
c11gene = cluster11[,7]
print(c11gene)
genes_to_check = c("HP","MUC5B","CHI3L1","CTSE","LINC00922","LTF"  )
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C12    Marker: "GMDS",C11     "ECM1","PLCG2""XIST","HSPA6",C12
c12gene = cluster12[,7]
print(c12gene)
genes_to_check = c("ECM1","PLCG2","XIST","HSPA6")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C13    Marker: "IGFBP5",C16     "DEFB1""FGB",C13
c13gene = cluster13[,7]
print(c13gene)
genes_to_check = c("DEFB1","IGFBP5","FGB")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C14    Marker: "CXCR4","CCL5","SRGN","CD52","TRAC","RGS1","CORO1A","CD3D" ,C14
c14gene = cluster14[,7]
print(c14gene)
genes_to_check = c("CXCR4","CCL5","SRGN","CD52","TRAC","RGS1","CORO1A","CD3D")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C15    Marker:"MDFI","MS4A15""TAC3",C15
c15gene = cluster15[,7]
print(c15gene)
genes_to_check = c("TAC3","MDFI","MS4A15")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C16    Marker:"C9orf24","C20orf85","C1orf194","FAM183A","RSPH1","TMEM190","C5orf49","C11orf88",C16
c16gene = cluster16[,7]
print(c16gene)
genes_to_check = c("C9orf24","C20orf85","C1orf194","FAM183A","RSPH1","TMEM190","C5orf49","C11orf88")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C17    Marker:"BPIFB2","TFF1","BPIFA1","ZG16B",C17
c17gene = cluster17[,7]
print(c17gene)
genes_to_check = c("BPIFB2","TFF1","BPIFA1","ZG16B")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()


# C18    Marker:"SNHG25","SLC26A2","ASRGL1","KLK11",C18
c18gene = cluster18[,7]
print(c18gene)
genes_to_check = c("SNHG25","SLC26A2","ASRGL1","KLK11" )
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()





#C0,C3,C4,C5,C7,C10,C11,C13：Basal Cell    Marker:"KRT17","GPX2",    "NTS","NIPSNAP2",    "TFPI2",   "NMU",   "ASS1", "BIRC3",     "LCN2","SLK","SELENOW",  "LTF"  , "DEFB1","FGB"
#C1：Pulmonary Alveolar Type II Cell       Marker:"PEBP4 ,"PLA2G1B" 
#C2：FOXN4+ Cell                           Marker:"CDKN2A" 
#C6：Luminal Epithelial Cell               Marker:"TCN1" 
#C8,C9：SLC16A7+ Cell                      Marker:"XAF1" ,"AGER"
#C12:Ionocyte Cell                         Marker:"PLCG2"
#C14:Langerhans Cell                       Marker:"RGS1","SRGN","CORO1A"
#C15:Unkonw
#C16:Ciliated Cell                         Marker:"C9orf24","C20orf85","C1orf194","FAM183A"
#C17,C18:Secretory Cell                    Marker:"ZG16B","SLC26A2"

new.cluster.ids <- c("Basal Cell", "Pulmonary Alveolar Type II Cell", "FOXN4+ Cell", "Basal Cell", "Basal Cell", "Basal Cell", "Luminal Epithelial Cell", 
                     "Basal Cell", "SLC16A7+ Cell", "SLC16A7+ Cell","Basal Cell","Basal Cell","Ionocyte Cell","Basal Cell","Langerhans Cell","Unkonw",
                     "Ciliated Cell","Secretory Cell","Secretory Cell")
library(ggplot2)
names(new.cluster.ids) <- levels(Normal_Epi)
Normal_Epi <- RenameIdents(Normal_Epi, new.cluster.ids)
DimPlot(Normal_Epi, reduction = "umap", label = TRUE, pt.size = 0.5,raster=FALSE)
ggsave(filename = "Normal_Epi-Sub-UMAP_Cellmarker.png",width = 14,height = 10,path = "../../Fig/Step1-4/")


genes_to_check = c("KRT17","GPX2","NIPSNAP2","TFPI2",
                   "PEBP4" ,"PLA2G1B","CDKN2A","TCN1","XAF1" ,"AGER","PLCG2","RGS1","SRGN","CORO1A","C9orf24","C20orf85","C1orf194","FAM183A","ZG16B","SLC26A2")
DotPlot(Normal_Epi,features = unique(genes_to_check)) + RotatedAxis()
ggsave(filename = "Normal_Epi-Sub-All Marker.png",width = 18,height = 12,path = "../../Fig/Step1-4/")

save(Normal_Epi,file = "res0.4_Normal_Epi.Rdata")#保存Rdata，用于后续分析


# library(dplyr)
# MonMarker <- FindMarkers(MonDC, ident.1 = c(0,2,5,6,7,9,11), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
# MonMarker_top10 <-MonMarker  %>% top_n(n = 10, wt = avg_log2FC)
# 
# MacMarker <- FindMarkers(MonDC, ident.1 = c(1,12), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
# MacMarker_top10 <-MacMarker  %>% top_n(n = 10, wt = avg_log2FC)
# 
# DCMarker <- FindMarkers(MonDC, ident.1 = c(3,4), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
# DCMarker_top10 <-DCMarker  %>% top_n(n = 10, wt = avg_log2FC)


# 2022.6.20  美化图片
library(Seurat)
library(ggsci)
library(ggplot2)
load(file = "res0.4_Epi.Rdata")
pal = pal_ucscgb(alpha = 0.7)(7)

Idents(object = Epi) <- "copykat.pred"
table(Epi@active.ident)
DimPlot(Epi, reduction = "umap", label = T, label.box = T, cols = c(pal[5],pal[4],pal[1]), raster = T)  + NoLegend()
ggsave(filename = 'UMAP_Epi_Cpaykat.pdf',width = 8,height = 6,path = '../../Fig/Step1-4/')


load(file = 'res0.4_Normal_Epi.Rdata')
pal2 = pal_ucscgb(alpha = 0.5)(10)
table(Normal_Epi@active.ident)
DimPlot(Normal_Epi, reduction = "umap", label = T, label.box = T, cols = pal2, raster = T)  + NoLegend()
ggsave(filename = 'UMAP_Normal_Epi_9type.pdf',width = 10,height = 6,path = '../../Fig/Step1-4/')

Normal_Epi@meta.data[["Rename"]] = Normal_Epi@active.ident
table(Normal_Epi@meta.data$Rename)

Idents(object = Normal_Epi) <- "RNA_snn_res.0.2"
table(Normal_Epi@active.ident)
pal3 = pal_ucscgb(alpha = 0.7)(20)
DimPlot(Normal_Epi, reduction = "umap", label = T, label.box = T, cols = pal3, raster = T)  + NoLegend()
ggsave(filename = 'UMAP_Normal_Epi_res0.2.pdf',width = 8,height = 6,path = '../../Fig/Step1-4/')


Idents(object = Normal_Epi) <- "Rename"
table(Normal_Epi@active.ident)

save(Normal_Epi,file = 'res0.4_Normal_Epi.Rdata')

#C0,C3,C4,C5,C7,C10,C11,C13：Basal Cell    Marker:"KRT17","GPX2",    "NTS","NIPSNAP2",    "TFPI2",   "NMU",   "ASS1", "BIRC3",     "LCN2","SLK","SELENOW",  "LTF"  , "DEFB1","FGB"
#C1：Pulmonary Alveolar Type II Cell       Marker:"PEBP4 ,"PLA2G1B" 
#C2：FOXN4+ Cell                           Marker:"CDKN2A" 
#C6：Luminal Epithelial Cell               Marker:"TCN1" 
#C8,C9：SLC16A7+ Cell                      Marker:"XAF1" ,"AGER"
#C12:Ionocyte Cell                         Marker:"PLCG2"
#C14:Langerhans Cell                       Marker:"RGS1","SRGN","CORO1A"
#C15:Unkonw
#C16:Ciliated Cell                         Marker:"C9orf24","C20orf85","C1orf194","FAM183A"
#C17,C18:Secretory Cell                    Marker:"ZG16B","SLC26A2"
genes_to_check = c("KRT17","GPX2","NIPSNAP2","TFPI2",
                   "PEBP4" ,"PLA2G1B","CDKN2A","TCN1","XAF1" ,"AGER","PLCG2","RGS1","SRGN",
                   "CORO1A","C9orf24","C20orf85","C1orf194","FAM183A","ZG16B","SLC26A2")
DotPlot(Normal_Epi,features = unique(genes_to_check), cols = c("lightgrey", "orange")) + RotatedAxis()
ggsave(filename = "All Marker_Normal_Epi.pdf",width = 12,height = 8,path = "../../Fig/Step1-4/")


# 'GPX2','TFPI2'   ,'PEBP4',"CDKN2A",'TCN1',  "XAF1",  'PLCG2','CORO1A','C9orf24','SLC26A2'
# genes_to_check = c( 'GPX2','TFPI2'   ,'PEBP4',"CDKN2A",'TCN1',  "XAF1",  'PLCG2','CORO1A','C9orf24','SLC26A2')
# FeaturePlot(Normal_Epi, features = unique(genes_to_check), cols = c("lightgrey", "orange"),ncol = 5)
# ggsave(filename = "All Marker_Normal_Epi2.png",width = 22,height = 9,path = "../2*GSE_fig/res0.4/Marker")

