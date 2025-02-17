
## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)
gc()

# 设置工作路径
setwd("/home/datahup/tyh/SCRNA/Rdata/")


library(Seurat)
load(file = "../Step1-1/res0.4_NSCLC.Integrate.Rdata")
#  提取感兴趣的细胞簇进行亚聚类  # C0,C1,C10,C15,C19,,C2,C16,C22  Plasma_cell+T Cell+B Cell
TBcell  = NSCLC.Integrate[,NSCLC.Integrate@meta.data$seurat_clusters %in% c(0,1,10,15,19,2,16,22)]#提取淋巴系免疫细胞簇
rm(NSCLC.Integrate)
gc()

TBcell <- NormalizeData(TBcell)
TBcell <- FindVariableFeatures(TBcell, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(TBcell), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(TBcell)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "TBcell-Sub-Top10-VarGene.png",width = 20,height = 10,path = "../../Fig/Step1-2/")


# 数据归一化 + 线形降维
all.genes <- rownames(TBcell)
TBcell <- ScaleData(TBcell, features = all.genes)
# 线性降维 对缩放的数据执行PCA.默认情况下，只使用前面确定的变量特性作为输入，但是如果想选择不同的子集，可以使用features参数来定义。
TBcell <- RunPCA(TBcell, features = VariableFeatures(object = TBcell))
# Examine and visualize PCA results a few different ways
# 查看PCA结果
print(TBcell[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(TBcell, dims = 1:2, reduction = "pca")
ggsave(filename = "TBcell-Sub-PCA.png",width = 16,height = 10,path = "../../Fig/Step1-2/")

DimPlot(TBcell, reduction = "pca", raster=FALSE)
ggsave(filename = "TBcell-Sub-PCA2.png",width = 16,height = 10,path = "../../Fig/Step1-2/") 

DimHeatmap(TBcell, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
ggsave(filename = "TBcell-Sub-PC1_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-2/")

DimHeatmap(TBcell, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
ggsave(filename = "TBcell-Sub-PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-2/")

save(TBcell,file = "res0.4_TBcell.Rdata")


load(file = "res0.4_TBcell.Rdata")
TBcell <- JackStraw(TBcell, num.replicate = 100) # 01h 39m 55s
TBcell <- ScoreJackStraw(TBcell, dims = 1:20)
JackStrawPlot(TBcell, dims = 1:20)
ggsave(filename = "res0.4_TBcell_JackPlot.png",width = 16,height = 10,path = "../../Fig/Step1-2/")
ElbowPlot(TBcell)#肘部图 
# 综合以上方法，选择15个主成成分作为参数用于后续分析。
ggsave(filename = "ElbowPlot.png",width = 12,height = 10,path = "../../Fig/Step1-2/")


#细胞聚类 KNN算法
library(clustree)
TBcell <- FindNeighbors(TBcell, dims = 1:15)#dims = 1:15 即选取前15个主成分来分类细胞。
TBcell <- FindClusters(object = TBcell,
                       resolution = c(seq(0,1,by = 0.1)))
clustree(TBcell@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "TBcell-Sub-resolution(0-1)(1).png",width = 20,height = 14,path = "../../Fig/Step1-2/")


#选取resolution = 0.4 作为后续分析参数
# Assign identity of clusters
Idents(object = TBcell) <- "RNA_snn_res.0.4"
TBcell@meta.data$seurat_clusters = TBcell@meta.data$RNA_snn_res.0.4
head(Idents(TBcell), 5)#查看前5个细胞的分类ID


# 非线性降维 UMAP
# UMAP
TBcell <- RunUMAP(TBcell, dims = 1:15)
DimPlot(TBcell, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "TBcell-Sub-UMAP-label-res0.4(1).png",width = 18,height = 12,path = "../../Fig/Step1-2/")



# 找每个簇的差异基因
library(dplyr)
TBcellMarker <- FindAllMarkers(TBcell, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <-TBcellMarker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(TBcell, features = top10$gene) + NoLegend()
ggsave(filename = "TBcell-Sub-Top10-MarkerGene.png",width = 16,height = 10,path = "../../Fig/Step1-2/")

save(top10,TBcellMarker,file = "res0.4_TBcell-Sub-Markers.Rdata")
save(TBcell, file = "res0.4_TBcell.Rdata")#保存Rdata，用于后续分析



load(file = "res0.4_TBcell.Rdata")
# 细胞注释
sce = TBcell
library(celldex)
library(SingleR)
library(Seurat)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.4
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


pro='TB cell-Sub-SingleR_anno_res0.4'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 18,height = 14,path = "../../Fig/Step1-2/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 18,height = 14,path = "../../Fig/Step1-2/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 18,height = 14,path = "../../Fig/Step1-2/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 18,height = 14,path = "../../Fig/Step1-2/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 18,height = 14,path = "../../Fig/Step1-2/")

TBcell = sce

save(top10,TBcellMarker,cellType,file = "res0.4_TBcell-Sub-Markers.Rdata")
save(TBcell, file = "res0.4_TBcell.Rdata")#保存Rdata，用于后续分析

write.csv(cellType,file = "TBcell-Sub-cellType_res0.4.csv")
write.csv(TBcellMarker, file =  "TBcell-Sub-markers_res0.4.csv")
write.csv(top10, file = "TBcell-Sub-top10gene_res0.4.csv")



load(file = "res0.4_TBcell.Rdata")

# VlnPlot(TBcell, features = "SPOCK2")
# ggsave(filename = "Dendritic Cell_CNP.png",width = 18,height = 12,path = "./Fig/Step1-2/")
# 
# FeaturePlot(TBcell, features = "IL7R")
# ggsave(filename = "Dendritic Cell_CNP(1).png",width = 18,height = 12,path = "./Fig/Step1-2/")

# C0,C9 : Natural killer (NK) Cell            Marker:"FGFBP2","GNLY","FCGR3A","XCL1","XCL2","KLRC1"
# C1,C5 : CD4+ T Memory Cell                  Marker:"NOSIP","CCR7"
# C2,C8 : CD8+ T Cell                         Marker:"GZMK","CD8A","CD8B"
# C3,C6,C13 : B Cell                          Marker:"MS4A1","CD79A","VPREB3"
# C4,C11 : Natural killer T (NKT) Cell        Marker:"LMNA","RGCC","TUBA1B","TUBB","HIST1H4C","TYMS"
# C7 : Regulatory T (Treg) cell               Marker:"TIGIT","BATF","CORO1B"
# C10,C12 : Plasma Cell                       Marker:"IGLC2","IGKC","JCHAIN","IGHG4","IGHG1"


new.cluster.ids <- c("NK Cell", "CD4+ T Memory Cell", "CD8+ T Cell", "B Cell", "NKT Cell", "CD4+ T Memory Cell", "B Cell", 
                     "Regulatory T (Treg) Cell", "CD8+ T Cell", "NK Cell","Plasma Cell","NKT Cell",
                     "Plasma Cell","B Cell")

names(new.cluster.ids) <- levels(TBcell)
TBcell <- RenameIdents(TBcell, new.cluster.ids)
DimPlot(TBcell, reduction = "umap", label = TRUE, pt.size = 0.5,raster=FALSE)
ggsave(filename = "TBcell-Sub-UMAP_Cellmarker.png",width = 14,height = 10,path = "../../Fig/Step1-2/")

save(TBcell, file = "res0.4_TBcell.Rdata")

genes_to_check = c("FGFBP2","GNLY","FCGR3A","NOSIP","CCR7","GZMK","CD8A","CD8B","MS4A1","CD79B","VPREB3",
                   "LMNA","RGCC","TUBA1B","TUBB","TIGIT","BATF","CORO1B","JCHAIN","IGHG4","IGHG1")
DotPlot(TBcell,features = unique(genes_to_check)) + RotatedAxis()
ggsave(filename = "TBcell-Sub-All Marker.png",width = 18,height = 12,path = "../../Fig/Step1-2/")



# load(file = "../../NSCLC/Rproject/NSCLC/2*GSE_project/2*GSE_project/res0.4_AllMarkergene.Rdata")
genes_to_check = c("CD74"  ,   "CD79A" ,   "HLA-DPA1", "HLA-DPB1", "HLA-DQA1" ,"HLA-DQB1" ,"HLA-DRA"  ,"MS4A1"   ,
                   "SMIM14" ,  "VPREB3","CD74"  ,   "CD79A"   , "CD79B"   , "CD83"   ,  "HLA-DPB1" ,"HLA-DQB1", "HLA-DRA" , "HLA-DRB1",
                   "IGHD"   ,  "TCL1A" )
DotPlot(TBcell,features = unique(genes_to_check)) + RotatedAxis()
ggsave(filename = "TBcell-Sub-All Marker.png",width = 18,height = 12,path = "../../Fig/Step1-2/")




# 2022.6.20  美化图片
load(file = 'res0.4_TBcell.Rdata')
library(ggsci)
library(Seurat)
library(ggplot2)
pal = pal_ucscgb()(7)
DimPlot(TBcell, reduction = "umap", label = T, label.box = T, cols = pal) + NoLegend()
ggsave(filename = 'UMAP_TBCell_7type.pdf',width = 8,height = 6,path = '../../Fig/Step1-2/')


TBcell@meta.data[["Rename"]] = TBcell@active.ident

Idents(object = TBcell) <- "RNA_snn_res.0.4"
table(TBcell@active.ident)
pal2 = pal_ucscgb()(14)
DimPlot(TBcell, reduction = "umap", label = T, label.box = T, cols = pal2) + NoLegend()
ggsave(filename = 'UMAP_TBCell_res0.4.pdf',width = 8,height = 6,path = '../../Fig/Step1-2/')

Idents(object = TBcell) <- "Rename"
table(TBcell@active.ident)

save(TBcell,file = 'res0.4_TBcell.Rdata')


# C0,C9 : Natural killer (NK) Cell            Marker:"FGFBP2","GNLY","FCGR3A","XCL1","XCL2","KLRC1"
# C1,C5 : CD4+ T Memory Cell                  Marker:"NOSIP","CCR7"
# C2,C8 : CD8+ T Cell                         Marker:"GZMK","CD8A","CD8B"
# C3,C6,C13 : B Cell                          Marker:"MS4A1","CD79A","VPREB3"
# C4,C11 : Natural killer T (NKT) Cell        Marker:"LMNA","RGCC","TUBA1B","TUBB","HIST1H4C","TYMS"
# C7 : Regulatory T (Treg) cell               Marker:"TIGIT","BATF","CORO1B"
# C10,C12 : Plasma Cell                       Marker:"IGLC2","IGKC","JCHAIN","IGHG4","IGHG1"
genes_to_check = c("FGFBP2","GNLY","FCGR3A","NOSIP","CCR7","GZMK","CD8A","CD8B","MS4A1","CD79B","VPREB3",
                   "LMNA","RGCC","TUBA1B","TUBB","TIGIT","BATF","CORO1B","JCHAIN","IGHG4","IGHG1")
DotPlot(TBcell,features = unique(genes_to_check), cols = c("lightgrey", "orange")) + RotatedAxis()
ggsave(filename = "All Marker_TBcell.pdf",width = 12,height = 8,path = "../../Fig/Step1-2/")


# 'GNLY','CCR7','CD8A','MS4A1','RGCC','TIGIT','IGHG1'
# genes_to_check = c('GNLY','CCR7','CD8A','MS4A1','RGCC','TIGIT','IGHG1')
# FeaturePlot(TBcell, features = unique(genes_to_check), cols = c("lightgrey", "orange"),ncol = 4)
# ggsave(filename = "All Marker_TBcell2.png",width = 18,height = 9,path = "./Fig/Step1-2/")

