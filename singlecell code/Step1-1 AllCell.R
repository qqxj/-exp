### Step1-1  Single-Cell Analyze：All Cell
### 整理时间： 2022/7/22
### 作者： 庞建宇
###------


## 清空工作环境
rm(list = ls())
options(stringsAsFactors = F)


# 设置工作路径
setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step1-1/")
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)

## Step1 整合


# 1 GSE131907
scRNA1 <- readRDS("/home/pjy/NSCLC/NSCLC/Data/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
NSCLC1 <- CreateSeuratObject(counts = scRNA1, project = "NSCLC1",min.cells = 10, min.features = 500)
NSCLC1

# 2 GSE148071
fs_NSCLC7=list.files(path = "/home/pjy/NSCLC/NSCLC/Data/GSE148071",pattern = 'exp.txt.gz')
fs_NSCLC7
setwd("/home/pjy/NSCLC/NSCLC/Data/GSE148071")
sceList_NSCLC7 <-  lapply(fs_NSCLC7, function(x){
  a = read.table(x , quote = '')
  scRNA7 = a
  p=str_split(x,'_',simplify = T)[,2]
  NSCLC7 <- CreateSeuratObject(counts = scRNA7, project = p, min.cells = 10, min.features = 500)
})
NSCLC7 = sceList_NSCLC7

for (i in 1:length(NSCLC7)){
  names(NSCLC7)[i] <- paste("NSCLC7_", i, sep = "")
}

NSCLC7 <- merge(NSCLC7[["NSCLC7_1"]], y = c(NSCLC7[["NSCLC7_2"]], NSCLC7[["NSCLC7_3"]], NSCLC7[["NSCLC7_4"]], NSCLC7[["NSCLC7_5"]],
                                            NSCLC7[["NSCLC7_6"]], NSCLC7[["NSCLC7_7"]], NSCLC7[["NSCLC7_8"]], NSCLC7[["NSCLC7_9"]], NSCLC7[["NSCLC7_10"]],
                                            NSCLC7[["NSCLC7_11"]], NSCLC7[["NSCLC7_12"]], NSCLC7[["NSCLC7_13"]], NSCLC7[["NSCLC7_14"]], NSCLC7[["NSCLC7_15"]],
                                            NSCLC7[["NSCLC7_16"]], NSCLC7[["NSCLC7_17"]], NSCLC7[["NSCLC7_18"]], NSCLC7[["NSCLC7_19"]], NSCLC7[["NSCLC7_20"]],
                                            NSCLC7[["NSCLC7_21"]], NSCLC7[["NSCLC7_22"]], NSCLC7[["NSCLC7_23"]], NSCLC7[["NSCLC7_24"]], NSCLC7[["NSCLC7_25"]],
                                            NSCLC7[["NSCLC7_26"]], NSCLC7[["NSCLC7_27"]], NSCLC7[["NSCLC7_28"]], NSCLC7[["NSCLC7_29"]], NSCLC7[["NSCLC7_30"]],
                                            NSCLC7[["NSCLC7_31"]], NSCLC7[["NSCLC7_32"]], NSCLC7[["NSCLC7_33"]], NSCLC7[["NSCLC7_34"]], NSCLC7[["NSCLC7_35"]],
                                            NSCLC7[["NSCLC7_36"]], NSCLC7[["NSCLC7_37"]], NSCLC7[["NSCLC7_38"]], NSCLC7[["NSCLC7_39"]], NSCLC7[["NSCLC7_40"]],
                                            NSCLC7[["NSCLC7_41"]], NSCLC7[["NSCLC7_42"]]), add.cell.ids = c("7_1", "7_2", "7_3", "7_4", "7_5", "7_6", "7_7", "7_8", "7_9", "7_10",
                                                                                                            "7_11", "7_12", "7_13", "7_14", "7_15", "7_16", "7_17", "7_18", "7_19", "7_20",
                                                                                                            "7_21", "7_22", "7_23", "7_24", "7_25", "7_26", "7_27", "7_28", "7_29", "7_30",
                                                                                                            "7_31", "7_32", "7_33", "7_34", "7_35", "7_36", "7_37", "7_38", "7_39", "7_40", 
                                                                                                            "7_41", "7_42"),project = "NSCLC7")

NSCLC7

NSCLC_list <- list(NSCLC1 = NSCLC1, NSCLC2 = NSCLC7)

setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step1-1/")
save(NSCLC_list,file = "NSCLC_list.Rdata")

for (i in 1:length(NSCLC_list)) {
  NSCLC_list[[i]] <- NormalizeData(NSCLC_list[[i]])
  NSCLC_list[[i]] <- FindVariableFeatures(NSCLC_list[[i]], selection.method = "vst", nfeatures = 1500)
}

features <- SelectIntegrationFeatures(object.list = NSCLC_list, nfeatures = 1500)

NSCLC.anchors <- FindIntegrationAnchors(object.list = NSCLC_list, anchor.features = features)
save(NSCLC.anchors,features,file="NSCLC_anchors.Rdata")

##利用锚点整合数据，运行时间较长
NSCLC.Integrate <- IntegrateData(anchorset = NSCLC.anchors)
dim(NSCLC.Integrate)
NSCLC.Integrate
save(NSCLC.Integrate,file= "NSCLC_Integrate.Rdata")

load(file = "NSCLC_Integrate.Rdata")

NSCLC.Integrate
#An object of class Seurat 
#34160 features across 286987 samples within 2 assays 
#Active assay: integrated (1500 features, 1500 variable features)
#1 other assay present: RNA


#  对所有细胞进行重分类
#  获取细胞类型
Idents(object = NSCLC.Integrate)
levels(NSCLC.Integrate)
table(Idents(NSCLC.Integrate))  # 获取每个细胞类型的细胞数目表格
# 其他的一些细胞类型的处理
# Stash cell identity classes
NSCLC.Integrate[["old.ident"]] <- Idents(object = NSCLC.Integrate)
# Set identity classes
Idents(object = NSCLC.Integrate) <- "NSCLC_cell"
Idents(object = NSCLC.Integrate, cells = 1:286987) <- "NSCLC_cell"
Idents(object = NSCLC.Integrate)
levels(NSCLC.Integrate)
table(Idents(NSCLC.Integrate))
#NSCLC.Integrate[["old.ident"]] <- Idents(object = NSCLC.Integrate)
NSCLC.Integrate[["orig.ident"]] <- Idents(object = NSCLC.Integrate)



# Step2 QC
NSCLC.Integrate[["percent.mt"]] <- PercentageFeatureSet(NSCLC.Integrate, assay = "RNA",pattern = "^MT-")
head(NSCLC.Integrate@meta.data)
# nFeature_RNA代表每个细胞测到的基因数目。
# nCount_RNA代表每个细胞测到所有基因的表达量之和。
# percent.mt代表测到的线粒体基因的比例。
VlnPlot(NSCLC.Integrate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "QCbefore.png",width = 20, height = 12, path = "../../Fig/Step1-1/")

plot1 <- FeatureScatter(NSCLC.Integrate, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(NSCLC.Integrate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)
plot1 + plot2
ggsave(filename = "QCbefore2.png",width = 20, height = 12, path = "../../Fig/Step1-1/")

#根据上图进行质控
NSCLC.Integrate <- subset(NSCLC.Integrate, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 
                          & nCount_RNA > 200 & nCount_RNA < 35000)
NSCLC.Integrate
#An object of class Seurat 
#34160 features across 202424 samples within 2 assays 
#Active assay: integrated (1500 features, 1500 variable features)
#1 other assay present: RNA

VlnPlot(NSCLC.Integrate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "QCafter.png",width = 20, height = 12, path = "../../Fig/Step1-1/")

plot1 <- FeatureScatter(NSCLC.Integrate, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(NSCLC.Integrate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)
plot1 + plot2
ggsave(filename = "QCafter2.png",width = 20, height = 12, path = "../../Fig/Step1-1/")



# Step 3 标准化数据+找高变基因
DefaultAssay(NSCLC.Integrate) <- "RNA"
NSCLC.Integrate <- NormalizeData(NSCLC.Integrate)
NSCLC.Integrate <- FindVariableFeatures(NSCLC.Integrate, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(NSCLC.Integrate), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(NSCLC.Integrate)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "Top10-VarGene.png",width = 20,height = 10,path = "../../Fig/Step1-1/")



# Step 4 数据归一化 + 线形降维
all.genes <- rownames(NSCLC.Integrate)
NSCLC.Integrate <- ScaleData(NSCLC.Integrate, features = all.genes)
# 线性降维 对缩放的数据执行PCA.默认情况下，只使用前面确定的变量特性作为输入，但是如果想选择不同的子集，可以使用features参数来定义。
NSCLC.Integrate <- RunPCA(NSCLC.Integrate, features = VariableFeatures(object = NSCLC.Integrate))
# Examine and visualize PCA results a few different ways
# 查看PCA结果
print(NSCLC.Integrate[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(NSCLC.Integrate, dims = 1:2, reduction = "pca")
ggsave(filename = "PCA.png",width = 16,height = 10,path = "../../Fig/Step1-1/")

DimPlot(NSCLC.Integrate, reduction = "pca", raster=FALSE)
ggsave(filename = "PCA2.png",width = 16,height = 10,path = "../../Fig/Step1-1/") 

DimHeatmap(NSCLC.Integrate, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
ggsave(filename = "PC1_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-1/")

DimHeatmap(NSCLC.Integrate, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
ggsave(filename = "PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step1-1/")

NSCLC.Integrate <- JackStraw(NSCLC.Integrate, num.replicate = 100) # 03h 18m 04s 
NSCLC.Integrate <- ScoreJackStraw(NSCLC.Integrate, dims = 1:20)
JackStrawPlot(NSCLC.Integrate, dims = 1:20)
ggsave(filename = "JackPlot.png",width = 16,height = 10,path = "../../Fig/Step1-1/")

ElbowPlot(NSCLC.Integrate)#肘部图 
# 综合以上方法，选择12个主成成分作为参数用于后续分析。
ggsave(filename = "ElbowPlot.png",width = 12,height = 10,path = "../../Fig/Step1-1/")



# Step 5 细胞聚类 KNN算法
library(clustree)
NSCLC.Integrate <- FindNeighbors(NSCLC.Integrate, dims = 1:12)#dims = 1:12 即选取前12个主成分来分类细胞。
NSCLC.Integrate <- FindClusters(object = NSCLC.Integrate,
                                resolution = c(seq(0,1,by = 0.1)))
clustree(NSCLC.Integrate@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "resolution(0-1).png",width = 20,height = 14,path = "../../Fig/Step1-1/")


Idents(object = NSCLC.Integrate) <- "RNA_snn_res.0.4"
NSCLC.Integrate@meta.data$seurat_clusters = NSCLC.Integrate@meta.data$RNA_snn_res.0.4
head(Idents(NSCLC.Integrate), 5)#查看前5个细胞的分类ID




# Step 6 UMAP/TSNE
# UMAP
NSCLC.Integrate <- RunUMAP(NSCLC.Integrate, dims = 1:12)
DimPlot(NSCLC.Integrate, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "UMAP-label-res0.4.png",width = 18,height = 12,path = "../../Fig/Step1-1/")

# TSNE
NSCLC.Integrate <- RunTSNE(NSCLC.Integrate, dims = 1:12)
DimPlot(NSCLC.Integrate, reduction = "tsne", label = TRUE,raster=FALSE)
ggsave(filename = "TSNE-label-res0.4.png",width = 18,height = 12,path = "../../Fig/Step1-1/")




##Step 7  找每个簇的差异基因
library(dplyr)
NSCLC.markers <- FindAllMarkers(NSCLC.Integrate, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <-NSCLC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(NSCLC.Integrate, features = top10$gene) + NoLegend()
ggsave(filename = "Top10-MarkerGene.png",width = 16,height = 10,path = "../../Fig/Step1-1/")



save(top10,NSCLC.markers,file = "res0.4_Markers.Rdata")
save(NSCLC.Integrate, file = "res0.4_NSCLC.Integrate.Rdata")#保存Rdata，用于后续分析


##Step 8 自动化注释  SingleR
load(file = "res0.4_Markers.Rdata")
load(file = "res0.4_NSCLC.Integrate.Rdata")
sce = NSCLC.Integrate
library(celldex)
library(SingleR)
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


pro='SingleR_anno_res0.4'
DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Blue', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Blue.png'),width = 18,height = 14,path = "../../Fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_DICE', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_DICE.png'),width = 18,height = 14,path = "../../Fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_HPCA', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_HPCA.png'),width = 18,height = 14,path = "../../Fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Mona', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Mona.png'),width = 18,height = 14,path = "../../Fig/Step1-1/")

DimPlot(sce,reduction = "umap",label=T,group.by = 'singleR_Nover', raster=FALSE)
ggplot2::ggsave(filename = paste0(pro,'_UMAP_Nover.png'),width = 18,height = 14,path = "../../Fig/Step1-1/")

NSCLC.Integrate = sce

save(NSCLC.Integrate, file = "res0.4_NSCLC.Integrate.Rdata")

# write.csv(cellType,file = "cellType_res0.4.csv")
# write.csv(NSCLC.markers, file =  "NSCLC_markers_res0.2.csv")
# write.csv(top10, file = "top10gene_res0.2.csv")


load(file = "res0.4_NSCLC.Integrate.Rdata")
VlnPlot(NSCLC.Integrate, features = "LST1")
# ggsave(filename = "NQO1.png",width = 18,height = 12,path = "../../Fig/Step1-1//res0.4")
FeaturePlot(NSCLC.Integrate, features = "COL3A1")
# ggsave(filename = "NQO1(1).png",width = 18,height = 12,path = "../../Fig/Step1-1//res0.4")


# C0,C1,C10,C15.C19 : T Cell                                 Marker:"TRBC1","CD3D","CD3E"
# C2 : B cell                                                Marker:"CD79A","MS4A1","CD79B"
# C3 : Monocyte                                              Marker:"LST1","FCN1"
# C4,C5,C13,C20 : Dendritic Cell                             Marker:C1QB,C1QC,C1QA   C13:"MMP12","PLIN2","IGSF6"
# C6,C7,C8,C9,C12,C17,C23,C24,C26 : Epithelial/Cancer cell   Marker:KRT19,KRT18,KRT8
# C11,C21  Fibroblasts                                       Marker:"DCN","FBLN1","LUM","COL1A1" ,"COL3A1"
# C14 : Mast Cell                                            Marker:"TPSB2","TPSAB1","CPA3","CTSG"
# C16,C22  Plasma_cell                                       Marker:"IGLC2","IGHG1","IGHG3"
# C18,C27: Endothelial cell                                  Marker:"GNG11", "RAMP2" ,"VWF"
# C25  Oligodendrocytes                                      Marker:"PLP1", "MBP", "GPM6B", "APLP1"

new.cluster.ids <- c("T Cell", "T Cell", "B Cell", "Monocyte", "Dendritic Cell", "Dendritic Cell", "Epithelial/Cancer Cell", 
                     "Epithelial/Cancer Cell", "Epithelial/Cancer Cell", "Epithelial/Cancer Cell","T Cell","Fibroblast",
                     "Epithelial/Cancer Cell","Dendritic Cell","Mast Cell","T Cell","Plasma Cell","Epithelial/Cancer Cell",
                     "Endothelial Cell","T Cell","Dendritic Cell","Fibroblast","Plasma Cell","Epithelial/Cancer Cell",
                     "Epithelial/Cancer Cell","Oligodendrocytes","Epithelial/Cancer Cell","Endothelial Cell")
names(new.cluster.ids) <- levels(NSCLC.Integrate)
NSCLC.Integrate <- RenameIdents(NSCLC.Integrate, new.cluster.ids)
DimPlot(NSCLC.Integrate, reduction = "umap", label = TRUE, pt.size = 0.5,raster=FALSE)
ggsave(filename = "UMAP_Cellmarker.png",width = 14,height = 10,path = "../../Fig/Step1-1/")

DimPlot(NSCLC.Integrate, reduction = "tsne", label = TRUE, pt.size = 0.5,raster=FALSE) 
ggsave(filename = "TSNE_Cellmarker.png",width = 14,height = 10,path = "../../Fig/Step1-1/")

save(NSCLC.Integrate, file = "res0.4_NSCLC.Integrate.Rdata")


#NSCLC.Integrate@active.ident = NSCLC.Integrate@active.ident
genes_to_check = c("TRBC1","CD3D","CD3E","CD79A","MS4A1","CD79B","LST1","FCN1","C1QB","C1QC","C1QA",
                   "KRT19","KRT18","KRT8","DCN","FBLN1","LUM","COL1A1" ,"COL3A1","TPSB2","TPSAB1",
                   "CPA3","CTSG","IGLC2","IGHG1","IGHG3","GNG11", "RAMP2" ,"VWF","PLP1", "MBP", "GPM6B", "APLP1")
DotPlot(NSCLC.Integrate,features = unique(genes_to_check)) + RotatedAxis()
ggsave(filename = "All Marker.png",width = 18,height = 12,path = "../../Fig/Step1-1/")



# 2022.6.20  美化图片
load(file = 'res0.4_NSCLC.Integrate.Rdata')
load(file = 'res0.4_Markers.Rdata')
library(ggsci)
library(Seurat)
library(ggplot2)
pal = pal_ucscgb(alpha = 0.8)(26)
pal2 = pal_lancet(alpha = 0.8)(3)
pal3 = pal_ucscgb(alpha = 0.8)(10)


class(NSCLC.Integrate@active.ident)
table(NSCLC.Integrate@active.ident)
Idents(object = NSCLC.Integrate) <- "RNA_snn_res.0.4"
table(NSCLC.Integrate@active.ident)

DimPlot(NSCLC.Integrate, reduction = "umap", label = T, label.box = T, cols = c(pal,pal2)) + NoLegend()
ggsave(filename = 'UMAP_res0.4.pdf',width = 8,height = 6,path = '../../Fig/Step1-1/')


new.cluster.ids <- c("T Cell", "T Cell", "B Cell", "Monocyte", "Dendritic Cell", "Dendritic Cell", "Epithelial/Cancer Cell", 
                     "Epithelial/Cancer Cell", "Epithelial/Cancer Cell", "Epithelial/Cancer Cell","T Cell","Fibroblast",
                     "Epithelial/Cancer Cell","Dendritic Cell","Mast Cell","T Cell","Plasma Cell","Epithelial/Cancer Cell",
                     "Endothelial Cell","T Cell","Dendritic Cell","Fibroblast","Plasma Cell","Epithelial/Cancer Cell",
                     "Epithelial/Cancer Cell","Oligodendrocytes","Epithelial/Cancer Cell","Endothelial Cell")
names(new.cluster.ids) <- levels(NSCLC.Integrate)
NSCLC.Integrate <- RenameIdents(NSCLC.Integrate, new.cluster.ids)
table(NSCLC.Integrate@active.ident)

NSCLC.Integrate@meta.data[["First"]] = NSCLC.Integrate@active.ident


DimPlot(NSCLC.Integrate, reduction = "umap", label = T, label.box = T,cols = pal3) + NoLegend()
ggsave(filename = 'UMAP_10type.pdf',width = 8,height = 6,path = '../../Fig/Step1-1/')




genes_to_check = c("TRBC1","CD3D","CD3E","CD79A","MS4A1","CD79B","LST1","FCN1","C1QB","C1QC","C1QA",
                   "KRT19","KRT18","KRT8","DCN","FBLN1","LUM","COL1A1" ,"COL3A1","TPSB2","TPSAB1",
                   "CPA3","CTSG","IGLC2","IGHG1","IGHG3","GNG11", "RAMP2" ,"VWF","PLP1", "MBP", "GPM6B", "APLP1")

DotPlot(NSCLC.Integrate , features = genes_to_check,cols = c("lightgrey", "orange")) + RotatedAxis()
ggsave(filename = "10type_marker.pdf",width = 12,height = 8,path = "../../Fig/Step1-1/")


# CD3D、MS4A1、IGHG1、TPSB2、FCN1、C1QC、DCN、RAMP2、KRT19、PLP1
# genes_to_check = c('CD3D','MS4A1','IGHG1','TPSB2','FCN1','C1QC','DCN','RAMP2','KRT19','PLP1')
# FeaturePlot(NSCLC.Integrate, features = unique(genes_to_check), cols = c("lightgrey", "orange"),ncol = 5)
# ggsave(filename = "All Marker3.png",width = 26,height = 9,path = "../2*GSE_fig/res0.4/Marker")


# 
# Idents(object = NSCLC.Integrate) <- "Rename"
# table(NSCLC.Integrate@active.ident)
# 
# save(NSCLC.Integrate, file = "res0.4_NSCLC.Integrate.Rdata")



genes_to_check = c('CD177','LUCAT1','CSRNP1','CXCR2','MS4A7','RETN')
DotPlot(NSCLC.Integrate , features = genes_to_check,cols = c("lightgrey", "orange")) + RotatedAxis()



# 堆叠小提琴图
# library(Seurat)
# library(ggplot2)
# modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
#   p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
#     xlab("") + ylab(feature) + ggtitle("") +
#     theme(legend.position = "none",
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.ticks.y = element_line(),
#           axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
#           plot.margin = plot.margin )
#   return(p)
# }
# 
# ## main function
# StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
#   plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
#   plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
#     theme(axis.text.x=element_text(), axis.ticks.x = element_line())
#   p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
#   return(p)
# }
# 
# #配色方案
# my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
#                 '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
#                 '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
#                 '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
#                 '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
#                 '#968175')
# pdf(file = './RePlot(2022.6.20)/10type_marker.pdf',width = 12,height = 10)
# StackedVlnPlot(NSCLC.Integrate, genes_to_check, pt.size=0, cols=my36colors)
# dev.off()

load(file = 'res0.4_NSCLC.Integrate.Rdata')
table(Idents(NSCLC.Integrate))
dat = subset(NSCLC.Integrate, Idents(NSCLC.Integrate) != 'Unkonw')

genes_to_check = c('MS4A1', 'CD79B', 'VPREB3', # B cell
                   'KRT17', 'GPX2', 'NIPSNAP2', 'TFPI2',# Basal Cell
                   "WFDC2", "SAA1", "SCGB2A2","C6orf15", # Cancer Cell
                   'NOSIP', 'CCR7', # CD4+ T Memory Cell
                   'GZMK', 'CD8A', 'CD8B', # CD8+ T Cell
                   'C9orf24', 'C20orf85', 'C1orf194', 'FAM183A', # Ciliated Cell
                   'CD48', 'FCN1', # Dendritic Cell
                   'GNG11', 'RAMP2', 'VWF', # Endothelial Cell
                   'DCN', 'FBLN1', 'LUM', 'COL1A1', 'COL3A1', # Fibroblast
                   'CDKN2A',# FOXN4+ Cell
                   'KIAA0101', 'UBE2C', 'PTTG1', # Granulocyte−Monocyte Progenitor
                   'BIRC3', 'CCL17', 'CCL19', # Granulosa cell
                   'PLCG2', # Ionocyte Cell
                   'RGS1', 'SRGN', 'CORO1A', # Langerhans Cell
                   'TCN1', # Luminal Epithelial Cell
                   'FOLR2', 'PLTP', 'SEPP1', 'CCL18', # Macrophage
                   'TPSB2', 'TPSAB1', 'CPA3', 'CTSG', # Mast Cell
                   'FBP1', 'MARCO', 'RETN', # Monocyte
                   'PI3', 'IL1R2', 'LRG1', 'SLC25A37', 'CSF3R', # Neutrophil
                   'FGFBP2', 'GNLY', 'FCGR3A', # NK Cell
                   'LMNA', 'RGCC', 'TUBA1B', 'TUBB', # NKT cell
                   'PLP1', 'MBP', 'GPM6B', 'APLP1', # Oligodendrocytes
                   'JCHAIN', 'IGHG4', 'IGHG1', # Plasma Cell
                   'GZMB', 'JCHAIN', 'ITM2C', 'C12orf75', # Plasmacytoid Dendritic Cell
                   'PEBP4', 'PLA2G1B',# Pulmonary Alveolar Type II Cell
                   'TIGIT', 'BATF', 'CORO1B',  # Regulatory T (Treg) Cell
                   'ZG16B', 'SLC26A2', # Secretory Cell
                   'XAF1', 'AGER' # SLC16A7+ Cell
)
library(scRNAtoolVis)
library(ggsci)
# pal = pal_ucscgb(alpha = 0.8)(26)
# pal2 = pal_lancet(alpha = 0.8)(3)
pdf(file = '../../Fig/Step1-1/Markergene.pdf',width = 12,height = 16)
AverageHeatmap(object = NSCLC.Integrate,
               markerGene = genes_to_check,
               # htCol = c(pal,pal2)
               )
dev.off()
