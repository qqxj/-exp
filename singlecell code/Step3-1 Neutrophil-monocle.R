

### Step3-1  Single-Cell Analyze：Neutrophils Monocle
### 整理时间： 2022/7/18
### 作者： 庞建宇


rm(list = ls())
options(stringsAsFactors = F)

setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step3-1/')
library(DDRTree)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(monocle)



load(file = "../Step1-3/res0.4_MonDC.Rdata")
#  提取感兴趣的细胞簇进行亚聚类  #   Plasma_cell+T Cell+B Cell
table(MonDC@active.ident)
Neutrophil  = MonDC[,MonDC@active.ident %in% "Neutrophil"]
rm(MonDC)
gc()


Neutrophil <- NormalizeData(Neutrophil)
Neutrophil <- FindVariableFeatures(Neutrophil, selection.method = "vst", nfeatures = 2000)
# 查看最高变的10个基因
top10 <- head(VariableFeatures(Neutrophil), 10)
# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(Neutrophil)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "Neutrophil-Sub-Top10-VarGene.png",width = 20,height = 10,path = "../../Fig/Step3-1/")


# 数据归一化 + 线形降维
all.genes <- rownames(Neutrophil)
Neutrophil <- ScaleData(Neutrophil, features = all.genes)
# 线性降维 对缩放的数据执行PCA.默认情况下，只使用前面确定的变量特性作为输入，但是如果想选择不同的子集，可以使用features参数来定义。
Neutrophil <- RunPCA(Neutrophil, features = VariableFeatures(object = Neutrophil))
# Examine and visualize PCA results a few different ways
# 查看PCA结果
print(Neutrophil[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Neutrophil, dims = 1:2, reduction = "pca")
ggsave(filename = "Neutrophil-Sub-PCA.png",width = 16,height = 10,path = "../../Fig/Step3-1/")

DimPlot(Neutrophil, reduction = "pca", raster=FALSE)
ggsave(filename = "Neutrophil-Sub-PCA2.png",width = 16,height = 10,path = "../../Fig/Step3-1/") 

DimHeatmap(Neutrophil, dims = 1, cells = 500, balanced = TRUE)#1个PC 500个细胞
ggsave(filename = "Neutrophil-Sub-PC1_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step3-1/")

DimHeatmap(Neutrophil, dims = 1:15, cells = 500, balanced = TRUE)#15个PC
ggsave(filename = "Neutrophil-Sub-PC15_HeatmapPlot.png",width = 16,height = 10,path = "../../Fig/Step3-1/")


save(Neutrophil,file = "../Step3-1/res0.4_Neutrophil.Rdata")
load(file = "../Step3-1/res0.4_Neutrophil.Rdata")


#细胞聚类 KNN算法
library(clustree)
Neutrophil <- FindNeighbors(Neutrophil, dims = 1:20)#dims = 1:20 即选取前20个主成分来分类细胞。
Neutrophil <- FindClusters(object = Neutrophil,
                           resolution = c(seq(0,1,by = 0.1)))
clustree(Neutrophil@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "Neutrophil-Sub-resolution(0-1).png",width = 20,height = 14,path = "../../Fig/Step3-1/")


#选取resolution = 0.7 作为后续分析参数
# Assign identity of clusters
Idents(object = Neutrophil) <- "RNA_snn_res.0.7"
Neutrophil@meta.data$seurat_clusters = Neutrophil@meta.data$RNA_snn_res.0.7
head(Idents(Neutrophil), 5)#查看前5个细胞的分类ID


# 非线性降维 UMAP/TSNE
# UMAP
Neutrophil <- RunUMAP(Neutrophil, dims = 1:20)
DimPlot(Neutrophil, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(filename = "Neutrophil-Sub-UMAP-label.png",width = 8,height = 6,path = "../../Fig/Step3-1/")



# monocle  
# Step1 Seurat -> monocle
data <- as(as.matrix(Neutrophil@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Neutrophil@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


# Step 2 估计size factor和离散度
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
# 过滤低质量的细胞
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))


# Step 3 细胞分类
# Clustering cells without marker genes  无监督方法
HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)

# 可视化 tSNE
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "CellType",
                   markers = c("MYF5", "ANPEP"))
HSMM <- clusterCells(HSMM, num_clusters = 10)
plot_cell_clusters(HSMM)


# Step 4 构建轨迹
# Step 4.1: 选择定义过程的基因 三种方法均为无监督

#使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(Neutrophil)
# diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
# HSMM <- setOrderingFilter(HSMM, diff.genes)
# plot_ordering_genes(HSMM)

##使用seurat选择的高变基因
# var.seurat <- VariableFeatures(Neutrophil)
# # var.seurat <- Neutrophil@assays[["RNA"]]@var.features
# HSMM <- setOrderingFilter(HSMM, var.seurat)
# plot_ordering_genes(HSMM)

##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


# Step 4.2 降维
HSMM <- reduceDimension(HSMM, max_components = 3,
                        method = 'DDRTree')

# Step 4.3 按照轨迹排序细胞
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")+
  facet_wrap(~seurat_clusters, nrow = 1)
ggsave(filename = 'seurat_clusters.png',width = 10,height = 5,path = '../../Fig/Step3-1/')


plot_cell_trajectory(HSMM, color_by = "State")

plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 1)
ggsave(filename = 'State.png',width = 10,height = 5,path = '../../Fig/Step3-1/')



# p1 = plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
# p2 = plot_cell_trajectory(HSMM, color_by = "State")
# p3 = plot_cell_trajectory(HSMM, color_by = "Pseudotime")
# 
# plotc <- p1|p2|p3
# plotc
# ggsave(filename = 'Clu+State+Tra.png',width = 10,height = 5,plot = plotc,path = '/home/pjy/NSCLC/NSCLC/Rproject/NSCLC/2*GSE_project/2*GSE_fig/res0.4/Neutrophil-monocle')


library(dplyr)
Marker <- unique(deg.cluster$gene)
Time_diff <- differentialGeneTest(HSMM[Marker,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# num_clusters为人为设置的聚类
p=plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
p
ggsave("Time_Marker.png", p, width = 10, height = 12,path = '../../Fig/Step3-1/')

# 提取每个簇的基因
p$tree_row

clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)




# 判断哪种状态对应于快速扩散
# Neutrophil   Marker:"PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R"
# blast_genes <- row.names(subset(fData(HSMM),
#                                 gene_short_name %in% c( "LINC01681", "SPATA22","CMC2","CCL4")))
# plot_genes_jitter(HSMM[blast_genes,],
#                   grouping = "State",
#                   min_expr = 0.1)
# 
# HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
#                                           num_cells_expressed >= 10))
# HSMM_filtered <- HSMM[HSMM_expressed_genes,]
# my_genes <- row.names(subset(fData(HSMM_filtered),
#                              gene_short_name %in% c("PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R")))
# cds_subset <- HSMM_filtered[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
# 
# plot_genes_in_pseudotime(cds_subset, color_by =  "State")
# 
# genes <- c("PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R")
# png(filename = "PI3+IL1R2+LRG1+SLC25A37+CSF3R_State_2.png",width = 1200,height = 800)
# p1 <- plot_genes_jitter(HSMM[genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(HSMM[genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(HSMM[genes,], color_by = "State")
# plotc <- p1|p2|p3
# dev.off()




# Step 5 差异分析
# Finding Genes that Change as a Function of Pseudotime
# to_be_tested <- row.names(subset(fData(HSMM),
#                                  gene_short_name %in% c("LINC01681", "SPATA22","CMC2","CCL4")))
# cds_subset <- HSMM[to_be_tested,]
# 
# diff_test_res <- differentialGeneTest(cds_subset,
#                                       fullModelFormulaStr = "~sm.ns(Pseudotime)")
# 
# diff_test_res[,c("gene_short_name", "pval", "qval")]
# plot_genes_in_pseudotime(cds_subset, color_by ="State")



#这里是把排序基因（disp.genes）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
library(dplyr)
Time_diff <- differentialGeneTest(HSMM[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# num_clusters为人为设置的聚类
p=plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)

ggsave("Time_heatmapAll.png", p, width = 10, height = 12,path = "../../Fig/Step3-1/")

# 提取每个簇的基因
p$tree_row
# Call:
#   hclust(d = d, method = method)
# 
# Cluster method   : ward.D2 
# Number of objects: 604 
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
# clustering
# 1   2   3 
# 181 261 162 
write.csv(clustering, "Time_clustering_all.csv", row.names = F)



# Finding Genes that Distinguish Cell Type or State
to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("UBC", "NCAM1", "ANPEP")))
cds_subset <- HSMM[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~State")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset,
                  grouping = "State",
                  color_by = "State",
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)


# Step 6 单细胞轨迹的“分支”分析
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "State")
BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)




## pData(HSMM)取出的是HSMM对象中HSMM@phenoData@data的内容
library(ggpubr)
library(Seurat)
df <- pData(HSMM) 
table(df$State)
df$NewState <- ifelse(df$State == '1', 'State1',
                  ifelse(df$State == '2','State2',
                      ifelse(df$State == '3', 'State3',
                         ifelse(df$State =='4','State2','State4'))))
table(df$NewState)

State <- data.frame(State = df$NewState, row.names = rownames(df), Cell = rownames(df))
State <- State[order(State$State),] 


table(State$State)
Neutrophil <- AddMetaData(Neutrophil, metadata = State)
Idents(object = Neutrophil) <- "State"
head(Idents(Neutrophil), 5)#查看前5个细胞的分类ID
table(Idents(Neutrophil))

Statemarker <- FindAllMarkers(Neutrophil)
table(Statemarker$cluster)
top50 <-Statemarker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(Neutrophil, features = top50$gene) + NoLegend()
ggsave(filename = 'State_Marker50.png',width = 8,height = 12,path = '../../Fig/Step3-1/')
save(Statemarker,file = 'Neutrophil_Statemarker.Rdata')

library(dplyr)
gene0 <- unique(Statemarker$gene)
Time_diff <- differentialGeneTest(HSMM[gene0,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# num_clusters为人为设置的聚类
p <- plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters = 4,show_rownames=T, return_heatmap=T)
ggsave(filename = 'State_Marker_pseudotime.png',plot = p,width = 8,height = 12,path = '../../Fig/Step3-1/')

# 提取每个簇的基因
# p$tree_row
# clusters <- cutree(p$tree_row, k = 4)
# clustering <- data.frame(clusters)
# clustering[,1] <- as.character(clustering[,1])
# colnames(clustering) <- "Gene_Clusters"
# table(clustering)

# clu1 <- subset(clustering,clustering$Gene_Clusters == 1)
# clu2 <- subset(clustering,clustering$Gene_Clusters == 2)
# clu3 <- subset(clustering,clustering$Gene_Clusters == 3)
# clu4 <- subset(clustering,clustering$Gene_Clusters == 4)
# clu5 <- subset(clustering,clustering$Gene_Clusters == 5)

# table(Statemarker$cluster)
# Statemarker1 <- subset(Statemarker,Statemarker$cluster== 1)
# Statemarker2 <- subset(Statemarker,Statemarker$cluster== 2)
# Statemarker3 <- subset(Statemarker,Statemarker$cluster== 3)
# Statemarker4 <- subset(Statemarker,Statemarker$cluster== 4)
# Statemarker5 <- subset(Statemarker,Statemarker$cluster== 5)


HSMM@phenoData@data$State <- df$NewState
p1 <- plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
p2 <- plot_cell_trajectory(HSMM, color_by = "State")
p3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plotc <- p1|p2|p3
plotc
ggsave(filename = 'Clu+State+Tra.png',width = 12,height = 5,plot = plotc,path = '../../Fig/Step3-1/')


plot_cell_trajectory(HSMM,color_by = "State")+facet_wrap(~State,nrow=1)
ggsave(filename = 'State.png',width = 10,height = 5,path = '../../Fig/Step3-1/')


my_genes <- row.names(subset(fData(HSMM),
                             gene_short_name %in% c('CD177', 'CSRNP1', 'CXCR2', 'LUCAT1', 'MS4A7', 'RETN')))
cds_subset <- HSMM[my_genes,]


plot_genes_in_pseudotime(cds_subset, color_by =  "State")
plot_genes_jitter(cds_subset, grouping = "State", color_by = "State")
plot_genes_violin(cds_subset, grouping = "State", color_by = "State")

save(Neutrophil,HSMM,file = "res0.4_Neutrophil.Rdata")
# kegg_list <- NULL
# go_list <- NULL
# for (i in 1:6) {
# 
#   dat <- paste0('State','_',i)
#   dat <- subset(ALLmarker,ALLmarker$cluster == i)
#   dat$SYMBOL <- dat$gene
# 
# ENTREZID1 <- bitr(unique(dat$SYMBOL), fromType = "SYMBOL",
#                 toType = c("ENTREZID"),
#                 OrgDb = org.Hs.eg.db)
# 
# kk.all <- enrichKEGG(gene = ENTREZID1$ENTREZID,
#                      organism = 'hsa',
#                      pvalueCutoff = 0.9,
#                      qvalueCutoff = 0.9) 
# kegg_list[[i]] <- kk.all@result
# 
# 
# titlename <- paste0('State',i)
# 
# dotplot(kk.all,showCategory=10, title = paste0(titlename,'_','Enrichment: KEGG pathway'))
# ggsave(filename = paste0(titlename,'_','Enrichment: KEGG1 pathway.png'),width = 12,height = 8 ,path = '/home/pjy/NSCLC/NSCLC/fig/Neutrophil/6State_kegg')
# 
# barplot(kk.all,showCategory=10, title ="Enrichment: KEGG pathway")
# ggsave(filename =paste0(titlename,'_','Enrichment: KEGG2 pathway.png'),width = 12,height = 8,path = '/home/pjy/NSCLC/NSCLC/fig/Neutrophil/6State_kegg')
# 
# 
# go.all <- enrichGO(gene           = ENTREZID1$ENTREZID, 
#                    OrgDb          = org.Hs.eg.db,
#                    ont            = 'all', 
#                    pAdjustMethod  = "BH",
#                    pvalueCutoff   = 0.9,
#                    qvalueCutoff = 0.9,
#                    readable       = TRUE)
# 
# go_list[[i]] <- go.all@result
# 
# dotplot(go.all, split="ONTOLOGY",showCategory = 5,
#         title =paste0(titlename,'_','Enrichment: GO pathway'))+ facet_grid(ONTOLOGY~., scale="free") 
# ggsave(filename = paste0(titlename,'_','Enrichment: GO pathway.png'),width = 12,height = 8,path = '/home/pjy/NSCLC/NSCLC/fig/Neutrophil/6State_go')
# 
# }


# 通路富集

setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step3-1/')
load(file = 'Neutrophil_Statemarker.Rdata')

library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
GO_GSEA_list <- list()
KEGG_GSEA_list <- list()

for (i in 1:4) {
  
name <- paste0('State',i)
dat <- subset(Statemarker,Statemarker$cluster == name)
dat$SYMBOL <- dat$gene

ENTREZID1 <- bitr(unique(dat$SYMBOL), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)

dat <- subset(dat,dat$SYMBOL %in% ENTREZID1$SYMBOL)
dat <- merge(dat,ENTREZID1)
geneList=dat$avg_log2FC
names(geneList)=dat$ENTREZID
geneList=sort(geneList,decreasing = T)

#GO GSEA
GO <- gseGO(
  geneList, #gene_fc
  ont = "ALL",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,#人类注释基因
  keyType = "ENTREZID",
  pvalueCutoff = 0.9,
  pAdjustMethod = "BH")#p值校正方法

sortGO<-GO[order(GO$enrichmentScore, decreasing = T),]#按照enrichment score从高到低排序
head(sortGO)
dim(sortGO)
# write.table(sortGO,paste0(paste0('State',i),'_','gsea_sortGO.txt')) #保存结果
GO_GSEA_list[[i]] <- sortGO
GO_GSEA_list[[i+4]] <- GO

#KEGG GSEA
KEGG <- gseKEGG(
  geneList,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.9,
  pAdjustMethod = "BH")


sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]#按照enrichment score从高到低排序
head(sortKEGG)
dim(sortKEGG)
# write.table(sortKEGG,paste0(paste0('State',i),'_','gsea_sortKEGG.txt')) #保存结果
KEGG_GSEA_list[[i]] <- sortKEGG
KEGG_GSEA_list[[i+4]] <- KEGG

}

save(KEGG_GSEA_list,GO_GSEA_list,file = 'State_GSEA.Rdata')
load(file = 'State_GSEA.Rdata')

# 查看富集情况-KEGG
library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggsci)
pal = pal_ucscgb()(20)
State1_KEGG <- as.data.frame(KEGG_GSEA_list[1]) %>% filter(pvalue < 0.05)
write.csv(State1_KEGG,file = 'State1_KEGG.csv')
paths <- State1_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[5]],paths, pvalue_table = TRUE,title = 'State1 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:9])
ggsave(filename = 'State1 KEGG Enrichment.png',width = 12,height = 8,path = '../../Fig/Step3-1/')


State2_KEGG <- as.data.frame(KEGG_GSEA_list[2]) %>% filter(pvalue < 0.05)
write.csv(State2_KEGG,file = 'State2_KEGG.csv')
paths <- State2_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:2])
# paths <- c("hsa05171", "hsa03010")#选取你需要展示的通路ID
# gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment')
ggsave(filename = 'State2 KEGG Enrichment.png',width = 12,height = 8,path = '../../Fig/Step3-1/')


State3_KEGG <- as.data.frame(KEGG_GSEA_list[3]) %>% filter(pvalue < 0.05)
write.csv(State3_KEGG,file = 'State3_KEGG.csv')
paths <- State3_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:8])
# paths <- c("hsa05140", "hsa04613","hsa03010","hsa05171")#选取你需要展示的通路ID
# gseaplot2(KEGG_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 KEGG Enrichment')
ggsave(filename = 'State3 KEGG Enrichment.png',width = 12,height = 8,path = '../../Fig/Step3-1/')


State4_KEGG <- as.data.frame(KEGG_GSEA_list[4]) %>% filter(pvalue < 0.05)
write.csv(State4_KEGG,file = 'State4_KEGG.csv')
paths <- State4_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:6])
# paths <- c("hsa05323", "hsa05171","hsa04144","hsa03010")#选取你需要展示的通路ID
# gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment')
ggsave(filename = 'State4 KEGG Enrichment.png',width = 12,height = 8,path = '../../Fig/Step3-1/')


# 多个向量取交集
KEGG_inter <- Reduce(intersect,list(State1_KEGG$Description,State2_KEGG$Description,State3_KEGG$Description,State4_KEGG$Description))
KEGG_inter
# [1] "Coronavirus disease - COVID-19"

venn_list_KEGG <- list(State1 = State1_KEGG$Description,
                       State2 = State2_KEGG$Description,
                       State3 = State3_KEGG$Description,
                       State4 = State4_KEGG$Description)
# Venn
library(VennDiagram)
venn.diagram(venn_list_KEGG, filename = 'KEGG_inter.png', imagetype = 'png', 
             fill = c('red', 'blue','green','Purple'), alpha = 0.50, 
             cat.col = c('red', 'blue','green','Purple'), cat.cex = 0.9, cat.fontfamily = 'serif',
             col = c('red', 'blue','green','Purple'), cex = 1.5, fontfamily = 'serif')
inter_kegg <- get.venn.partitions(venn_list_KEGG) # 提取交集基因

# UpSet
library(UpSetR)
data2<-fromList(venn_list_KEGG)
g1<- upset(data2,nset = 4,order.by = "freq",
           queries = list(list(query = intersects, params = list("State1","State2","State3","State4"),
                               active = T, color = 'red'))) # 要先加载UpSetR
g1


# 查看富集情况-GO
State1_GO <- as.data.frame(GO_GSEA_list[1]) %>% filter(pvalue < 0.05)
State1_GO <- State1_GO[order(State1_GO$ONTOLOGY),]
table(State1_GO$ONTOLOGY)
write.csv(State1_GO,file = 'State1_GO.csv')

State1_GO_BP <- subset(State1_GO,State1_GO$ONTOLOGY == 'BP')
paths <- c("GO:0006396", "GO:2000116","GO:0002526","GO:0003012")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[5]],paths, pvalue_table = TRUE,title = 'State1 GO BP Enrichment')
ggsave(filename = 'State1 GO BP Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State1_GO_CC <- subset(State1_GO,State1_GO$ONTOLOGY == 'CC')
paths <- c("GO:0016604", "GO:0062023","GO:0031012","GO:0030312")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[5]],paths, pvalue_table = TRUE,title = 'State1 GO CC Enrichment')
ggsave(filename = 'State1 GO CC Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State1_GO_MF <- subset(State1_GO,State1_GO$ONTOLOGY == 'MF')
paths <- c("GO:0005215", "GO:0048306","GO:0005125","GO:0005509")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[5]],paths, pvalue_table = TRUE,title = 'State1 GO MF Enrichment')
ggsave(filename = 'State1 GO MF Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')



State2_GO <- as.data.frame(GO_GSEA_list[2]) %>% filter(pvalue < 0.05)
State2_GO <- State2_GO[order(State2_GO$ONTOLOGY),]
table(State2_GO$ONTOLOGY)
write.csv(State2_GO,file = 'State2_GO.csv')

State2_GO_BP <- subset(State2_GO,State2_GO$ONTOLOGY == 'BP')
paths <- c("GO:0006935", "GO:0006629","GO:0002181","GO:0006412")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 GO BP Enrichment')
ggsave(filename = 'State2 GO BP Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State2_GO_CC <- subset(State2_GO,State2_GO$ONTOLOGY == 'CC')
paths <- c("GO:0031226", "GO:0005887","GO:0044391","GO:1990904")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 GO CC Enrichment')
ggsave(filename = 'State2 GO CC Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State2_GO_MF <- subset(State2_GO,State2_GO$ONTOLOGY == 'MF')
paths <- c("GO:0005102", "GO:0003735","GO:0005198","GO:0003723")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 GO MF Enrichment')
ggsave(filename = 'State2 GO MF Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')



State3_GO <- as.data.frame(GO_GSEA_list[3]) %>% filter(pvalue < 0.05)
State3_GO <- State3_GO[order(State3_GO$ONTOLOGY),]
table(State3_GO$ONTOLOGY)
write.csv(State3_GO,file = 'State3_GO.csv')

State3_GO_BP <- subset(State3_GO,State3_GO$ONTOLOGY == 'BP')
paths <- c("GO:0050729", "GO:0009612","GO:0002181","GO:0048732")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 GO BP Enrichment')
ggsave(filename = 'State3 GO BP Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State3_GO_CC <- subset(State3_GO,State3_GO$ONTOLOGY == 'CC')
paths <- c("GO:0015629", "GO:0005856","GO:0022625","GO:0005840")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 GO CC Enrichment')
ggsave(filename = 'State3 GO CC Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State3_GO_MF <- subset(State3_GO,State3_GO$ONTOLOGY == 'MF')
paths <- c("GO:0015631", "GO:0003924","GO:0003735","GO:0005198")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 GO MF Enrichment')
ggsave(filename = 'State3 GO MF Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')



State4_GO <- as.data.frame(GO_GSEA_list[4]) %>% filter(pvalue < 0.05)
State4_GO <- State4_GO[order(State4_GO$ONTOLOGY),]
table(State4_GO$ONTOLOGY)
write.csv(State4_GO,file = 'State4_GO.csv')

State4_GO_BP <- subset(State4_GO,State4_GO$ONTOLOGY == 'BP')
paths <- c("GO:0070372", "GO:0070371","GO:0009615","GO:0006397")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 GO BP Enrichment')
ggsave(filename = 'State4 GO BP Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State4_GO_CC <- subset(State4_GO,State4_GO$ONTOLOGY == 'CC')
paths <- c("GO:0005788", "GO:0012507","GO:0030667","GO:0005887")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 GO CC Enrichment')
ggsave(filename = 'State4 GO CC Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')

State4_GO_MF <- subset(State4_GO,State4_GO$ONTOLOGY == 'MF')
paths <- c("GO:0005125", "GO:0005126","GO:0004888","GO:0060089")#选取你需要展示的通路ID
gseaplot2(GO_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 GO MF Enrichment')
ggsave(filename = 'State4 GO MF Enrichment.png',width = 10,height = 8,path = '../../Fig/Step3-1/')


# 多个向量取交集
# GO_BP_inter <- Reduce(intersect,list(State1_KEGG$Description,State2_KEGG$Description,State3_KEGG$Description,State4_KEGG$Description))
# KEGG_inter
# # [1] "Coronavirus disease - COVID-19"


# Venn
# BP
venn_list_GO_BP <- list(State1 = State1_GO_BP$Description,
                        State2 = State2_GO_BP$Description,
                        State3 = State3_GO_BP$Description,
                        State4 = State4_GO_BP$Description)
library(VennDiagram)
venn.diagram(venn_list_GO_BP, filename = 'GO_BP_inter.png', imagetype = 'png', 
             fill = c('red', 'blue','green','Purple'), alpha = 0.50, 
             cat.col = c('red', 'blue','green','Purple'), cat.cex = 0.9, cat.fontfamily = 'serif',
             col = c('red', 'blue','green','Purple'), cex = 1.5, fontfamily = 'serif')
inter_go_bp <- get.venn.partitions(venn_list_KEGG) # 提取交集基因

# UpSet
library(UpSetR)
data2<-fromList(venn_list_GO_BP)
upset(data2,nset = 4,order.by = "freq",
           queries = list(list(query = intersects, params = list("State1","State2","State3","State4"),
                               active = T, color = 'red'))) # 要先加载UpSetR

# CC
venn_list_GO_CC <- list(State1 = State1_GO_CC$Description,
                        State2 = State2_GO_CC$Description,
                        State3 = State3_GO_CC$Description,
                        State4 = State4_GO_CC$Description)
venn.diagram(venn_list_GO_CC, filename = 'GO_CC_inter.png', imagetype = 'png', 
             fill = c('red', 'blue','green','Purple'), alpha = 0.50, 
             cat.col = c('red', 'blue','green','Purple'), cat.cex = 0.9, cat.fontfamily = 'serif',
             col = c('red', 'blue','green','Purple'), cex = 1.5, fontfamily = 'serif')
inter_go_cc <- get.venn.partitions(venn_list_KEGG) # 提取交集基因

# UpSet
library(UpSetR)
data2<-fromList(venn_list_GO_CC)
upset(data2,nset = 4,order.by = "freq",
      queries = list(list(query = intersects, params = list("State1","State2","State3","State4"),
                          active = T, color = 'red'))) # 要先加载UpSetR


# MF
venn_list_GO_MF <- list(State1 = State1_GO_MF$Description,
                        State2 = State2_GO_MF$Description,
                        State3 = State3_GO_MF$Description,
                        State4 = State4_GO_MF$Description)
venn.diagram(venn_list_GO_MF, filename = 'GO_MF_inter.png', imagetype = 'png', 
             fill = c('red', 'blue','green','Purple'), alpha = 0.50, 
             cat.col = c('red', 'blue','green','Purple'), cat.cex = 0.9, cat.fontfamily = 'serif',
             col = c('red', 'blue','green','Purple'), cex = 1.5, fontfamily = 'serif')
inter_go_mf <- get.venn.partitions(venn_list_KEGG) # 提取交集基因

# UpSet
library(UpSetR)
data2<-fromList(venn_list_GO_MF)
upset(data2,nset = 4,order.by = "freq",
      queries = list(list(query = intersects, params = list("State1","State2","State3","State4"),
                          active = T, color = 'red'))) # 要先加载UpSetR

# 5维以上Venn图
# library(venn)
# #作图
# # png('venn_7.png', width = 1500, height = 1500, res = 200, units = 'px')
# venn(venn_list_KEGG, zcolor = 'style')
# # dev.off()



# 2022.6.8 
load(file = 'res0.4_Neutrophil.Rdata') 
load(file = 'Neutrophil_Statemarker.Rdata')
# p2 <- plot_cell_trajectory(HSMM, color_by = "State")
# p3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
# p2|p3
# ggsave(filename = 'State_Pseudotime.png',width = 10,height = 6,path = '../2*GSE_fig/res0.4/Neutrophil-monocle/')
# 
# plot_cell_trajectory(HSMM, color_by = "State") + facet_wrap(~State,nrow=1)
# ggsave(filename = 'State_Pseudotime.png',width = 10,height = 6,path = '../2*GSE_fig/res0.4/Neutrophil-monocle/')
# 
# library(dplyr)
# library(Seurat)
# top50 <-Statemarker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# DoHeatmap(Neutrophil, features = top50$gene) + NoLegend()
# ggsave(filename = 'State_Marker50.png',width = 10,height = 16,path = '/home/pjy/NSCLC/NSCLC/Rproject/NSCLC/2*GSE_project/2*GSE_fig/res0.4/Neutrophil-monocle')
# 


# 2022.6.20 重新作图
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step3-1/')
load(file = 'res0.4_Neutrophil.Rdata')
load(file = 'Neutrophil_Statemarker.Rdata')

library(ggsci)
p1 = plot_cell_trajectory(HSMM, 
                          color_by = "Pseudotime", # color_by 选择分组类型
                          show_branch_points = F,  # show_branch_points 是否展示分枝节点
                          show_tree = T)  +        # show_tree 是否展示连线
  scale_color_steps(low = 'yellow', high = 'purple',)
p1
ggsave(filename = 'Monocle_Pseudotime.pdf',plot = p1,width = 8,height = 6,path = '../../Fig/Step3-1/')

pal = pal_ucscgb(alpha = 0.8)(10)
p2 = plot_cell_trajectory(HSMM, color_by = "State",show_branch_points = F,show_tree = T) + 
  scale_color_discrete(names('')) +  # 隐藏图例标题
  scale_color_manual(values = pal[c(1,6,3,4)]) 

p2
ggsave(filename = 'Monocle_Rename.pdf',plot = p2,width = 8,height = 6,path = '../../Fig/Step3-1/')

p4 = p2 + facet_wrap(~State, nrow = 1) # ~State 为指定的分组
p4
ggsave(filename = 'Monocle_State.pdf',plot = p4,width = 12,height = 6,path = '../../Fig/Step3-1/')


library(dplyr)
library(Seurat)
top50 <-Statemarker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
pdf(file = '../../Fig/Step3-1//State_Marker50.pdf',width = 10,height = 18)
DoHeatmap(Neutrophil, features = top50$gene) 
dev.off()



load(file = 'State_GSEA.Rdata')
# 查看富集情况-KEGG
library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggsci)
pal = pal_ucscgb()(20)
State1_KEGG <- as.data.frame(KEGG_GSEA_list[1]) %>% filter(pvalue < 0.05)
paths <- State1_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[5]],paths, pvalue_table = TRUE,title = 'State1 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:9])
ggsave(filename = 'State1 KEGG Enrichment.pdf',width = 12,height = 8,path = '../../Fig/Step3-1/')


State2_KEGG <- as.data.frame(KEGG_GSEA_list[2]) %>% filter(pvalue < 0.05)
paths <- State2_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:2])
# paths <- c("hsa05171", "hsa03010")#选取你需要展示的通路ID
# gseaplot2(KEGG_GSEA_list[[6]],paths, pvalue_table = TRUE,title = 'State2 KEGG Enrichment')
ggsave(filename = 'State2 KEGG Enrichment.pdf',width = 12,height = 8,path = '../../Fig/Step3-1/')


State3_KEGG <- as.data.frame(KEGG_GSEA_list[3]) %>% filter(pvalue < 0.05)
paths <- State3_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[7]],paths, pvalue_table = TRUE,title = 'State3 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:8])
ggsave(filename = 'State3 KEGG Enrichment.pdf',width = 12,height = 8,path = '../../Fig/Step3-1/')


State4_KEGG <- as.data.frame(KEGG_GSEA_list[4]) %>% filter(pvalue < 0.05)
paths <- State4_KEGG$ID#选取你需要展示的通路ID
gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment',
          rel_heights = c(2, 0.3,0.3),subplots = c(1:3), color = pal[1:6])
# paths <- c("hsa05323", "hsa05171","hsa04144","hsa03010")#选取你需要展示的通路ID
# gseaplot2(KEGG_GSEA_list[[8]],paths, pvalue_table = TRUE,title = 'State4 KEGG Enrichment')
ggsave(filename = 'State4 KEGG Enrichment.pdf',width = 12,height = 8,path = '../../Fig/Step3-1/')


# 查看富集情况-GO
s1 = GO_GSEA_list[[1]]
s1 = s1 %>% filter(pvalue < 0.01)
s1a = s1[,c(1:3,5,7)]
s1a$State = rep('State1', length(rownames(s1a)))


s2 = GO_GSEA_list[[2]]
s2 = s2 %>% filter(pvalue < 0.01)
s2a = s2[,c(1:3,5,7)]
s2a$State = rep('State2', length(rownames(s2a)))


s3 = GO_GSEA_list[[3]]
s3 = s3 %>% filter(pvalue < 0.01)
s3a = s3[,c(1:3,5,7)]
s3a$State = rep('State3', length(rownames(s3a)))


s4 = GO_GSEA_list[[4]]
s4 = s4 %>% filter(pvalue < 0.01)
s4a = s4[,c(1:3,5,7)]
s4a$State = rep('State4', length(rownames(s4a)))


go_state = rbind(s1a,s2a,s3a,s4a)
go_bp = subset(go_state,go_state$ONTOLOGY == 'BP')
go_cc = subset(go_state,go_state$ONTOLOGY == 'CC')
go_mf = subset(go_state,go_state$ONTOLOGY == 'MF')


# BP
go_bp1 = go_bp %>% filter(abs(enrichmentScore) > 0.6) # BP途径过多 再次筛选
bp_s1 = subset(go_bp1,go_bp1$State == 'State1')
bp_s1 = bp_s1[,c(3,4)]
names(bp_s1)[2] = 'State1 ES'

bp_s2 = subset(go_bp1,go_bp1$State == 'State2')
bp_s2 = bp_s2[,c(3,4)]
names(bp_s2)[2] = 'State2 ES'

bp_s3 = subset(go_bp1,go_bp1$State == 'State3')
bp_s3 = bp_s3[,c(3,4)]
names(bp_s3)[2] = 'State3 ES'

bp_s4 = subset(go_bp1,go_bp1$State == 'State4')
bp_s4 = bp_s4[,c(3,4)]
names(bp_s4)[2] = 'State4 ES'


bp_union = full_join(bp_s1,full_join(bp_s2,full_join(bp_s3,bp_s4)))
rownames(bp_union) = bp_union$Description
bp_union = bp_union[,-1]
bp_union[is.na(bp_union)] = 0
str(bp_union)
bp_union = round(bp_union,digits = 3)


library(pheatmap)
p = pheatmap(bp_union,cluster_rows=F,cluster_cols=F,
             # display_numbers=T,number_format="%s.3f",
             border="white",
             fontsize_number=2,
             fontsize_col = 10,
             fontsize_row = 7,
             angle_col = 0,
             main = 'GSEA Enrichment Score in GO BP Terms',
             color = colorRampPalette(colors = c("purple","white","red"))(100),
             filename= '../../Fig/Step3-1/State1-4 GO BP union GSEA.pdf', width=8.5, height=6)



# CC
cc_s1 = subset(go_cc,go_cc$State == 'State1')
cc_s1 = cc_s1[,c(3,4)]
names(cc_s1)[2] = 'State1 ES'
# rownames(bp_s1) = bp_s1$Description
# bp_s1 = bp_s1[,-1]

cc_s2 = subset(go_cc,go_cc$State == 'State2')
cc_s2 = cc_s2[,c(3,4)]
names(cc_s2)[2] = 'State2 ES'
# rownames(bp_s2) = bp_s2$Description
# bp_s2 = bp_s2[,-1]

cc_s3 = subset(go_cc,go_cc$State == 'State3')
cc_s3 = cc_s3[,c(3,4)]
names(cc_s3)[2] = 'State3 ES'

cc_s4 = subset(go_cc,go_cc$State == 'State4')
cc_s4 = cc_s4[,c(3,4)]
names(cc_s4)[2] = 'State4 ES'


cc_union = full_join(cc_s1,full_join(cc_s2,full_join(cc_s3,cc_s4)))
rownames(cc_union) = cc_union$Description
cc_union = cc_union[,-1]
cc_union[is.na(cc_union)] = 0
str(cc_union)
cc_union = round(cc_union,digits = 3)


library(pheatmap)
p = pheatmap(cc_union,cluster_rows=F,cluster_cols=F,
             # display_numbers=T,number_format="%s.3f",
             border="white",
             fontsize_number=2,
             fontsize_col = 10,
             fontsize_row = 8,
             angle_col = 0,
             main = 'GSEA Enrichment Score in GO CC Terms',
             color = colorRampPalette(colors = c("purple","white","red"))(100),
             filename= '../../Fig/Step3-1/State1-4 GO CC union GSEA.pdf', width=8, height=6)


# MF
mf_s1 = subset(go_mf,go_mf$State == 'State1')
mf_s1 = mf_s1[,c(3,4)]
names(mf_s1)[2] = 'State1 ES'
# rownames(bp_s1) = bp_s1$Description
# bp_s1 = bp_s1[,-1]

mf_s2 = subset(go_mf,go_mf$State == 'State2')
mf_s2 = mf_s2[,c(3,4)]
names(mf_s2)[2] = 'State2 ES'
# rownames(bp_s2) = bp_s2$Description
# bp_s2 = bp_s2[,-1]

mf_s3 = subset(go_mf,go_mf$State == 'State3')
mf_s3 = mf_s3[,c(3,4)]
names(mf_s3)[2] = 'State3 ES'

mf_s4 = subset(go_mf,go_mf$State == 'State4')
mf_s4 = mf_s4[,c(3,4)]
names(mf_s4)[2] = 'State4 ES'


mf_union = full_join(mf_s1,full_join(mf_s2,full_join(mf_s3,mf_s4)))
rownames(mf_union) = mf_union$Description
mf_union = mf_union[,-1]
mf_union[is.na(mf_union)] = 0
str(mf_union)
mf_union = round(mf_union,digits = 3)


library(pheatmap)
p = pheatmap(mf_union,cluster_rows=F,cluster_cols=F,
             # display_numbers=T,number_format="%s.3f",
             border="white",
             fontsize_number=2,
             fontsize_col = 10,
             fontsize_row = 8,
             angle_col = 0,
             main = 'GSEA Enrichment Score in GO MF Terms',
             color = colorRampPalette(colors = c("purple","white","red"))(100),
             filename= '../../Fig/Step3-1/State1-4 GO MF union GSEA.pdf', width=8, height=6)


