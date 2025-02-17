

### Step6-1  探索预后基因功能
### 整理时间： 2022/7/22
### 作者： 庞建宇
#setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step6-1/')
setwd("/home/datahup/syj/Single.Cell/Rproject/save/step-6-1")
rm(list = ls())
options(stringsAsFactors = F)


library(dplyr)
library(DESeq2)
load(file = "/home/datahup/syj/Single.Cell/Rproject/save/step-2-0/TCGA_Step1 output.Rdata")
load(file = "/home/datahup/syj/Single.Cell/Rproject/save/step 4-1/Train_Result.Rdata") 
gene = rownames(gene_multi)
NSCLCcount[1:6,1:6]
dat = as.data.frame(t(NSCLCcount))
dat[1:6,1:6]

dat0 = subset(dat,substr(rownames(dat),14,15) < 10 )#从每个行名中提取第14到第15个字符
df = data.frame(gene = dat0[,gene], row.names = rownames(dat0)) 


#Step1 高低分组做差异分析
for (i in 1:length(gene)) {
  
df = df[order(df[,i]),]
a = quantile(df[,i],c(0.3,0.7))
df$group = ifelse(df[,i] <= a[1],'Low',
                       ifelse(df[,i] >= a[2],'High', 'no'))
table(df$group)  
df1 = subset(df,df$group == 'High' | df$group == 'Low')
#构建dds矩阵
dat1 = as.data.frame(t(dat0[rownames(df1),]))
range(dat1)
dat1 <- 2^dat1 -1 #转Count
dat1 <- round(dat1,digits = 0) # 取整
dat1[dat1 < 0 ] = 0
# table(is.na(dat1))
# range(dat1)
# which(dat1 ==22378253076 ,arr.ind = TRUE)# 如果你想知道一个值的行和列matrix或者data.frame，可以考虑使用arr.ind=TRUE
# # row col
# # LINC00676 37819 531
# dat1 <- dat1[-37819,]
countData <- dat1
countData[1:6,1:6]
# table(is.na(countData))
Group <- data.frame(group = df1$group, row.names = rownames(df1))
condition <- factor(Group$group)
table(condition)
head(condition)
# table(is.na(condition))

# 差异分析
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
head(dds)
dim(dds)
dds <- DESeq(dds) 
resultsNames(dds)
res <- results(dds)
summary(res)
DEG_INFO <- as.data.frame(res)
# 提取差异分析结果
table(res$padj<0.05) #取P值小于0.05的结果
res <- res[order(res$padj),]
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
rownames(resdata) <- resdata$Row.names
resdata <- resdata[,-1]
resdata <- na.omit(resdata)
# 确定差异表达倍数
# logFC_cutoff <- with(resdata,mean(abs(resdata$log2FoldChange)) + 2*sd(abs(resdata$log2FoldChange)) )
# # 取前两位小数
# logFC_cutoff <- round(logFC_cutoff, 2)

diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)# 所有差异基因
diff_gene_deseq2_df <- resdata[diff_gene_deseq2,1:6]
diff_gene_deseq2_df$state <- ifelse(diff_gene_deseq2_df$log2FoldChange > 2, 'Up','Down')
table(diff_gene_deseq2_df$state)

write.csv(diff_gene_deseq2_df,file = paste0(gene[i],'_DEG.csv'))
# Up_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange > logFC_cutoff)) 
# Up_gene_deseq2 <- row.names(Up_gene_deseq2) # 上调差异基因
# 
# Down_gene_deseq2 <-  subset(res,padj < 0.05 & (log2FoldChange <  -logFC_cutoff))
# Down_gene_deseq2 <- row.names(Down_gene_deseq2) # 下调差异基因


# 火山图
library(ggplot2)
for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                          'padj' = res$padj,
                          'State' = rep('No', length(res$log2FoldChange)))
up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 2), which(for_volcano$padj < 0.05))
down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -2), which(for_volcano$padj < 0.05))
for_volcano[up_sig_indices,'State'] <- 'Up'
for_volcano[down_sig_indices,'State'] <- 'Down'
for_volcano$State <- as.factor(for_volcano$State)
for_volcano$padj <- -log10(for_volcano$padj)

this_tile <- paste0('Cutoff for logFC is 2',
                    '\nThe number of Up gene is ',nrow(diff_gene_deseq2_df[diff_gene_deseq2_df$state =='Up',]) ,
                    '\nThe number of Down gene is ',nrow(diff_gene_deseq2_df[diff_gene_deseq2_df$state =='Down',]))


p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = State))+
  geom_point(size = I(0.7))+
  scale_color_manual(values = c('No'='black', 'Up' = 'red', 'Down' = 'blue'))+
  geom_vline(xintercept = c(2, -2), lty=2, size=I(0.4), colour = 'grey11')+
  geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2, size=I(0.1),colour = 'grey11')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'),
        panel.grid = element_blank())+
  labs(x='log2FoldChange', y = '-log10Pvalue')+
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5)) 

p

ggsave(filename = paste0(gene[i],'_vol.png'),plot = p,width = 8,height = 6,path = "../../Fig/Step6-1/")

}



# 相关性分析
# for (i in 1:length(gene)) {
#   
# exprSet <- dat0
# y <- as.numeric(exprSet[,gene[i]])
# colnames <- colnames(exprSet)
# cor_data_df <- data.frame(colnames)
# 
#   for (j in 1:length(colnames)){
#   test <- cor.test(as.numeric(exprSet[,j]),y,type="spearman")
#   cor_data_df[j,2] <- test$estimate
#   cor_data_df[j,3] <- test$p.value
#   }
# 
# names(cor_data_df) <- c("symbol","correlation","pvalue")
# head(cor_data_df)
# cor_data_df <- na.omit(cor_data_df)
# 
# library(dplyr)
# library(tidyr)
# # 有显著差异的基因集
# cor_data_sig <- cor_data_df %>%
#   filter(pvalue < 0.05) %>%
#   arrange(desc(correlation))
# # %>%dplyr::slice(1:500)
# write.csv(cor_data_sig,file = paste0(gene[i],'_cor.csv'))
# 
# }



#  KEGG-GSVA
kk_list = list()
kegg_gsva_list = list()
go_list = list()
go_gsva_list = list()

# 加载GO通路数据库
{
source("getGoTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
save(GO_DATA, file = "GO_DATA.RData")
findGO <- function(pattern, method = "key"){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method == "key"){
    pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
  } else if(method == "gene"){
    pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways) = "pathway"
  
  if(length(pathways) == 0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
} # 用于寻找 GO ID
getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
} # 获取 GO geneSet
load("GO_DATA.RData") # 载入数据 GO_DATA
}


# gene = gene[-5] # CSRNP1 没有通路交集，需要单独run
for (i in 1:length(gene)) {

  print(paste0('Now is ',gene[i]))
  
# 分组
df = df[order(df[,i]),]
a = quantile(df[,i],c(0.3,0.7))
df[,paste0(gene[i],'_group')] = ifelse(df[,i] <= a[1],'Low',
                            ifelse(df[,i] >= a[2],'High', 'no'))
table(df[,paste0(gene[1],'_group')])  
df1 = subset(df,df[,paste0(gene[i],'_group')] == 'High' | df[,paste0(gene[i],'_group')] == 'Low')


# deg
df2 = read.csv(file = paste0(gene[i],'_DEG.csv'))
rownames(df2) <- df2$X
df2 <- df2[,-1]
#corgene
# df3 = read.csv(file = paste0(gene[i],'_cor.csv'))
# df3 <- df3[,-1] %>% filter(abs(correlation) >= 0.2) # 挑选绝对值相关性大于0.2的基因


## KEGG 富集
##  ID转换
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# degene 
deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db)
# KEGG
kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                     organism = 'hsa',
                     pvalueCutoff = 0.05)

#corgene
# cor_entr <- bitr(df3$symbol, fromType = "SYMBOL",
#                  toType = c("ENTREZID"),
#                  OrgDb = org.Hs.eg.db)
# # KEGG
# kk_cor <- enrichKEGG(gene = cor_entr$ENTREZID,
#                      organism = 'hsa',
#                      pvalueCutoff = 0.05)
# 
# kk_list[[gene[i]]] <- c(deg = kk_deg,cor = kk_cor)

kk_list[[gene[i]]] <- kk_deg


# 提取通路里的基因集
library(KEGGREST) 
# listDatabases()  

# library(VennDiagram)
# venn_list <- list(kegg_deg_pathway = kk_deg@result[["ID"]], kegg_cor_pathway = kk_cor@result[["ID"]])
# venn.diagram(venn_list, filename = paste0('../Figer/KO&Overexp/pathway_inter/',gene[i],'_kegginter.png'), imagetype = 'png', 
#              fill = c('red', 'blue'), alpha = 0.50, 
#              cat.col = c('red', 'blue'), cat.cex = 0.9, cat.fontfamily = 'serif',
#              col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif')

# inter <- get.venn.partitions(venn_list) # 提取交集基因
# keggpathway <- inter$..values..[1]
# keggpathway <- intersect(kk_deg@result[["ID"]],kk_cor@result[["ID"]])
keggpathway <- kk_deg@result[["ID"]]

kegg_genelist <- list()
   for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))

    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称

    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva

    }

# library(msigdbr)
library(GSVA)
# library(pheatmap)

meta <- data.frame(group = df1[,paste0(gene[i],'_group')], row.names = rownames(df1))
df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))

# m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #选取物种人类
# msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
kegg_gsva <- gsva(df3, kegg_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
# pheatmap(kegg_gsva, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/pathway_inter/',gene[i],'_kegg_gsva_heatmap.png'), width=14, height=16)#绘制热图

kegg_gsva_list[[gene[i]]] <- kegg_gsva


# GO
go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                   OrgDb          = org.Hs.eg.db,
                   ont            = 'ALL', 
                   pAdjustMethod  = "BH",
                   pvalueCutoff   = 0.05, 
                   # qvalueCutoff   = 0.2, 
                   readable       = TRUE)

# go_cor <- enrichGO(gene           = cor_entr$ENTREZID, 
#                    OrgDb          = org.Hs.eg.db,
#                    ont            = 'ALL', 
#                    pAdjustMethod  = "BH",
#                    pvalueCutoff   = 0.05, 
#                    # qvalueCutoff   = 0.2, 
#                    readable       = TRUE)

# go_list[[gene[i]]] <- c(deg = go_deg,cor = go_cor)
go_list[[gene[i]]] <- go_deg

# venn_list2 <- list(go_deg_pathway = go_deg@result[["ID"]], go_cor_pathway = go_cor@result[["ID"]])
# venn.diagram(venn_list2, filename = paste0('../Figer/KO&Overexp/pathway_inter/',gene[i],'_gointer.png'), imagetype = 'png', 
#              fill = c('red', 'blue'), alpha = 0.50, 
#              cat.col = c('red', 'blue'), cat.cex = 0.9, cat.fontfamily = 'serif',
#              col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif')

# inter <- get.venn.partitions(venn_list) # 提取交集基因
# keggpathway <- inter$..values..[1]
# gopathway <- intersect(go_deg@result[["ID"]],go_cor@result[["ID"]])
gopathway <- go_deg@result[["ID"]]

go_genelist <- list()
# gopathway <- go_cor@result[["ID"]]


# 批量获取通路基因集
   for (x in gopathway) {
  go_genelist <- getGO(gopathway)
   }

go_gsva <- gsva(df3, go_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
# go_gsva2 <- gsva(df3, go_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva

# pheatmap(go_gsva, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/pathway_inter/',gene[i],'_go_gsva_heatmap.png'), width=14, height=16)#绘制热图


go_gsva_list[[gene[i]]] <- go_gsva

print(paste0(gene[i],' is over'))

}
# save(df1,kk_list,kegg_gsva_list,go_list,go_gsva_list,file = 'Step9_Result.Rdata')

save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,file = 'Step6-1 Result.Rdata')



# load(file = 'Step9_Result.Rdata')
# gene = rownames(gene_multi)
# kegg_deg_gsva = list()
# go_deg_gsva = list()
# # 分别做差异基因与相关基因的gsva
# for (i in 1:length(gene)) {
#   
#   #degene
#   df2 = read.csv(file = paste0(gene[i],'_DEG.csv'))
#   rownames(df2) <- df2$X
#   df2 <- df2[,-1]
#   kegg_deg_list <- kk_list[[gene[i]]][["deg"]]@result[["ID"]]
#   
#   # 提取通路基因
#   # KEGG
#   library(KEGGREST) 
#   kegg_genelist <- list()
#   for (j in kegg_deg_list) {
#     gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
#     #获取通路中gene信息 
#     # gs[[1]]$GENE 
#     #查找所有基因 
#     # print(paste0('Now is ',i))
#     
#     genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
#     genelist <- genes[1:length(genes)%%3 ==2] 
#     gs[[1]]$NAME # 通路名称
#     
#     kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
#     
#   }
#   
#   # library(msigdbr)
#   library(GSVA)
#   library(pheatmap)
#   
#   meta <- data.frame(group = df1[,paste0(gene[i],'_group')], row.names = rownames(df1))
#   df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
#   # m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #选取物种人类
#   # msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
#   kegg_gsva <- gsva(df3, kegg_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
#   pheatmap(kegg_gsva, show_rownames=T, show_colnames=F, annotation_col=meta,
#            cluster_cols = F, # 列不聚类
#            # treeheight_row = 0,#不展示聚类树
#            fontsize_row=5, filename= paste0('../Figer/KO&Overexp/deg/',gene[i],'_kegg_gsva_heatmap.png'))#绘制热图
#   
#   kegg_deg_gsva[[gene[i]]] <- kegg_gsva
#   
#   
#   # GO
#   gopathway <- go_list[[gene[i]]][["deg"]]@result[["ID"]]
#   go_genelist <- list()
#   # gopathway <- go_cor@result[["ID"]]
#   # 批量获取通路基因集
#   for (x in gopathway) {
#     go_genelist <- getGO(gopathway)
#   }
#   
#   go_gsva <- gsva(df3, go_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
#   # go_gsva2 <- gsva(df3, go_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
#   
#   pheatmap(go_gsva, show_rownames=T, show_colnames=F, annotation_col=meta,
#            cluster_cols = F, # 列不聚类
#            # treeheight_row = 0,#不展示聚类树
#            fontsize_row=5, filename= paste0('../Figer/KO&Overexp/deg/',gene[i],'_go_gsva_heatmap.png'))#绘制热图
#   
#   
#   go_deg_gsva[[gene[i]]] <- go_gsva
#   
#   
#   #corgene
#   library(clusterProfiler)
#   library(org.Hs.eg.db)
#   library(ggplot2)
#   library(stringr)
#   p1 = dotplot(kk_list[[gene[i]]][["cor"]],
#                showCategory=15, title = paste0('Enrichment: ',
#                                                gene[i],' Correlation KEGG pathway'))+ scale_color_continuous(low='purple', high='yellow') + scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) 
#   ggsave(filename = paste0(gene[i],' Correlation KEGG pathway.png'),plot = p1,width = 10,height = 8,path = '../Figer/KO&Overexp/cor/')
#   
#   
#   p2 = barplot(kk_list[[gene[i]]][["cor"]],
#                showCategory=15, title = paste0('Enrichment: ',
#                                                gene[i],' Correlation KEGG pathway')) + 
#     # scale_color_continuous(low='purple', high='yellow') 
#       scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) 
#   ggsave(filename = paste0(gene[i],' Correlation KEGG pathway2.png'),plot = p2,width = 10,height = 8,path = '../Figer/KO&Overexp/cor/')
#   
#   
#   p3 = dotplot(go_list[[i]][["cor"]], split="ONTOLOGY",showCategory = 6,
#           title = paste0('Enrichment: ',gene[i],' Correlation GO pathway'))+ facet_grid(ONTOLOGY~., scale="free")
#   
#   p4 = p3 + scale_color_continuous(low='purple', high='yellow')+ scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) 
#   
#   ggsave(filename = paste0(gene[i],' Correlation GO pathway.png'),plot = p4,width = 10,height = 6,path = '../Figer/KO&Overexp/cor/')
#   
#   
#   
#   }
# save(df1,kk_list,kegg_gsva_list,go_list,go_gsva_list,kegg_deg_gsva,go_deg_gsva,file = 'Step9_Result.Rdata')


# 2022.6.10 
# load(file = 'Step6-1 Result.Rdata')
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ggplot2)
# library(stringr)
# # a = kk_list[["MS4A7"]]@result 
# # a = a[order(a$Count,decreasing = T),]
# 
# # 重新批量画图
# genes = c('MS4A7','CD177','CSRNP1','CXCR2','LUCAT1','RETN')
# for (i in 1:length(genes)) {
# # KEGG
# p1 = dotplot(kk_list[[genes[i]]],
#              showCategory=15, title = paste0('Enrichment: ' ,genes[i] ,' KEGG pathway'))+ scale_color_continuous(low='purple', high='yellow') + scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) 
# p1
# ggsave(filename = paste0(genes[i] ,' KEGG pathway.png'),plot = p1,width = 10,height = 8,path = '../Figer/KO&Overexp/2022.6.10/')
# 
# 
# p2 = barplot(kk_list[[genes[i]]],
#              showCategory=15, title = paste0('Enrichment: ' ,genes[i] ,' KEGG pathway')) +
# scale_color_continuous(low='purple', high='yellow')+
# scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) 
# p2
# ggsave(filename = paste0(genes[i] ,' KEGG pathway2.png'),plot = p2,width = 10,height = 8,path = '../Figer/KO&Overexp/2022.6.10/')
# 
# 
# # GO
# p3 = dotplot(go_list[[genes[i]]], split="ONTOLOGY",showCategory = 6,
#              title = paste0('Enrichment: ',genes[i],' GO pathway'))+ facet_grid(ONTOLOGY~., scale="free")
# 
# p4 = p3 + scale_color_continuous(low='purple', high='yellow')+ scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) 
# p4
# ggsave(filename = paste0(genes[i],' GO pathway.png'),plot = p4,width = 10,height = 6,path = '../Figer/KO&Overexp/2022.6.10/')
# 
# library(GOplot)
# b = go_list[[genes[i]]]@result
# term = b[,c(1:3,7,9)]
# names(term) = c('Category','ID','Term','adj_pval','Genes')
# term$Genes <- gsub("/", ",", term$Genes)  # 用逗号替换/
# # DEG
# FC = read.csv(file = paste0(genes[i],'_DEG.csv'))
# rownames(FC) = FC[,1]
# FC1 = FC[,c(1,3)]
# names(FC1) = c('ID','logFC')
# # 合并
# circ = circle_dat(term,FC1)
# #指定通路
# d3 = b[order(b$Count,decreasing = T),]
# d4 = d3[1:5,3]#这个是指定你要做那几个GO的弦图
# #最终画图对象
# chord = chord_dat(circ,FC1,d4)
# head(chord)
# #画图
# png(filename = paste0("../Figer/KO&Overexp/2022.6.10/",genes[i],'_GO_chord.png'),width = 1200,height = 1200)#准备好一块画布。
# GOChord(chord,space = 0.02,gene.order = 'logFC',gene.space = 0.25,gene.size = 5)#进行绘图。
# dev.off()#然后必须关闭并保存这块画布。
# 
# 
# 
# # GSVA
# df1 = subset(df,df[,paste0(genes[i],'_group')] != 'no')
# table(df1[,paste0(genes[i],'_group')])
# meta <- data.frame(group = df1[,paste0(genes[i],'_group')], row.names = rownames(df1))
# # GO GSVA
# df2 = go_gsva_list[[genes[i]]]
# df3 = go_list[[genes[i]]]@result
# df3 = df3[order(df3$Count,decreasing = T),]
# df3 = df3[1:20,3]
# # Top20
# library(pheatmap)
# df4 = df2[df3,]
# pheatmap(df4, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/2022.6.10/',genes[i],'_Top20_go_gsva_heatmap.png'), width=12, height=8)#绘制热图
# # All
# pheatmap(df2, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/2022.6.10/',genes[i],'_go_gsva_heatmap.png'), width=12, height=8)#绘制热图
# 
# 
# # KEGG GSVA
# df5 = kegg_gsva_list[[genes[i]]]
# df6 = kk_list[[genes[i]]]@result
# df6 = df6[order(df6$Count,decreasing = T),]
# df7 = df6[1:20,2]
# 
# # a = gsub(rownames(df5),'- Homo sapiens (human)','')
# # Maf <- gsub(rownames(df5), pattern = " - Homo sapiens (human)", replacement = '')
# rownames(df5) <- str_sub(rownames(df5), start = 1, end = -24)
# 
# # Top20
# library(pheatmap)
# df8 = subset(df5,rownames(df5) %in% df7)
# pheatmap(df8, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/2022.6.10/',genes[i],'_Top20_kegg_gsva_heatmap.png'), width=12, height=8)#绘制热图
# # All
# pheatmap(df5, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/2022.6.10/',genes[i],'_kegg_gsva_heatmap.png'), width=12, height=8)#绘制热图
# 
# }



# 基于所有基因GSVA
# load(file = "../Step2-0/TCGA_Step1 output.Rdata")
# load(file = "../Step4-1/Train_Result.Rdata") 
# 
# gene = rownames(gene_multi)
# NSCLCcount[1:6,1:6]
# dat = as.data.frame(t(NSCLCcount))
# dat[1:6,1:6]
# dat0 = subset(dat,substr(rownames(dat),14,15) < 10 )
# df = data.frame(gene = dat0[,gene], row.names = rownames(dat0)) 
# 
# library(KEGGREST) 
# source("./Code/getGoTerm.R")
# # GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
# # save(GO_DATA, file = "GO_DATA.RData")
# findGO <- function(pattern, method = "key"){
#   
#   if(!exists("GO_DATA"))
#     load("GO_DATA.RData")
#   if(method == "key"){
#     pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
#   } else if(method == "gene"){
#     pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
#   }
#   
#   colnames(pathways) = "pathway"
#   
#   if(length(pathways) == 0){
#     cat("No results!\n")
#   } else{
#     return(pathways)
#   }
# } # 用于寻找 GO ID
# getGO <- function(ID){
#   
#   if(!exists("GO_DATA"))
#     load("GO_DATA.RData")
#   allNAME = names(GO_DATA$PATHID2EXTID)
#   if(ID %in% allNAME){
#     geneSet = GO_DATA$PATHID2EXTID[ID]
#     names(geneSet) = GO_DATA$PATHID2NAME[ID]
#     return(geneSet)     
#   } else{
#     cat("No results!\n")
#   }
# } # 获取 GO geneSet
# load("GO_DATA.RData") # 载入数据 GO_DATA
# 
# 
# all_kegg_gsva_list = list()
# all_go_gsva_list = list()
# for (i in 1:length(gene)) {
#   
# keggpathway <- kk_list[[gene[i]]]@result$ID
# 
# kegg_genelist <- list()
# for (j in keggpathway) {
#   gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
#   #获取通路中gene信息 
#   # gs[[1]]$GENE 
#   #查找所有基因 
#   # print(paste0('Now is ',i))
#   
#   genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
#   genelist <- genes[1:length(genes)%%3 ==2] 
#   gs[[1]]$NAME # 通路名称
#   
#   kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
#   
# }
# 
# # library(msigdbr)
# library(GSVA)
# library(pheatmap)
# 
# meta <- data.frame(group = df[,paste0(gene[i],'_group')], row.names = rownames(df))
# table(meta$group)
# meta = subset(meta,meta$group != 'no')
# df3 <- as.matrix(t(dat0[rownames(meta),]))
# 
# # m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #选取物种人类
# # msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
# kegg_gsva <- gsva(df3, kegg_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
# pheatmap(kegg_gsva, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = T, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/2022.6.10.allgene/',gene[i],'_kegg_gsva_heatmap.png'), width=12, height=8)#绘制热图
# 
# all_kegg_gsva_list[[gene[i]]] <- kegg_gsva
# 
# #GO
# gopathway <- go_list[[gene[i]]]@result$ID
# 
# go_genelist <- list()
# # gopathway <- go_cor@result[["ID"]]
# # 批量获取通路基因集
# for (x in gopathway) {
#   go_genelist <- getGO(gopathway)
# }
# 
# go_gsva <- gsva(df3, go_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
# # go_gsva2 <- gsva(df3, go_genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
# 
# pheatmap(go_gsva, show_rownames=T, show_colnames=F, annotation_col=meta,
#          cluster_cols = T, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          fontsize_row=5, filename= paste0('../Figer/KO&Overexp/2022.6.10.allgene/',gene[i],'_go_gsva_heatmap.png'), width=12, height=8)#绘制热图
# 
# 
# all_go_gsva_list[[gene[i]]] <- go_gsva
# }
# 
# save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,all_kegg_gsva_list,all_go_gsva_list,file = 'Step6-1 Result.Rdata')




# load(file = 'Step6-1 Result.Rdata')
# load("../Step4-1/Train_Result.Rdata")
# # GSVA
# genes = rownames(gene_multi)
# for (i in 1:length(genes)) {
#   
# df1 = subset(df,df[,paste0(genes[i],'_group')] != 'no')
# table(df1[,paste0(genes[i],'_group')])
# meta <- data.frame(group = df1[,paste0(genes[i],'_group')], row.names = rownames(df1))
# # GO GSVA
# df2 = go_gsva_list[[genes[i]]]
# df3 = go_list[[genes[i]]]@result
# df3 = df3[order(df3$Count,decreasing = T),]
# df3 = df3[1:20,3]
# # Top20
# library(pheatmap)
# df4 = df2[df3,]
# ann_colors = list(
#   group = c(Low = "purple", High = "yellow")
# )
# pheatmap(df4, show_rownames=T, show_colnames=F, annotation_col=meta,annotation_colors=ann_colors,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          main = paste0(genes[i],' Top20 GO Terms'),
#          color = colorRampPalette(colors = c("purple","white","yellow"))(100),legend_breaks = c(-0.5:0.5),
#          fontsize_row=8, filename= paste0('../Figer/KO&Overexp/2022.6.10/',genes[i],'_Top20_go_gsva_heatmap.png'), width=10, height=6)#绘制热图
# 
# # KEGG GSVA
# df5 = kegg_gsva_list[[genes[i]]]
# df6 = kk_list[[genes[i]]]@result
# df6 = df6[order(df6$Count,decreasing = T),]
# df7 = df6[1:20,2]
# 
# # a = gsub(rownames(df5),'- Homo sapiens (human)','')
# # Maf <- gsub(rownames(df5), pattern = " - Homo sapiens (human)", replacement = '')
# rownames(df5) <- str_sub(rownames(df5), start = 1, end = -24)
# 
# # Top20
# library(pheatmap)
# df8 = subset(df5,rownames(df5) %in% df7)
# ann_colors = list(
#   group = c(Low = "purple", High = "yellow")
# )
# pheatmap(df8, show_rownames=T, show_colnames=F, annotation_col=meta,annotation_colors=ann_colors,
#          cluster_cols = F, # 列不聚类
#          # treeheight_row = 0,#不展示聚类树
#          main = paste0(genes[i],' Top20 KEGG Terms'),
#          color = colorRampPalette(colors = c("purple","white","yellow"))(100),legend_breaks = c(-0.5:0.5),
#          fontsize_row=8, filename= paste0('../Figer/KO&Overexp/2022.6.10/',genes[i],'_Top20_kegg_gsva_heatmap.png'), width=10, height=6)#绘制热图
# }



# 6.28 GSVA得分差异分析 找到关键通路
load(file = 'Step6-1 Result.Rdata')
load("../Step4-1/Train_Result.Rdata")
gene = rownames(gene_multi)

KEGG_logFC_list = list()
GO_logFC_list = list()

for (i in gene) {
  
kk_re = kk_list[[i]]@result
kk_re = kk_re %>% filter(pvalue < 0.05)
kk_gsva_df = as.data.frame(all_kegg_gsva_list[[i]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = kk_gsva_df[kk_re$Description,]

go_re = go_list[[i]]@result
go_re = go_re %>% filter(pvalue < 0.05)
go_gsva_df = as.data.frame(all_go_gsva_list[[i]])
go_gsva_df = go_gsva_df[go_re$Description,]


group = data.frame(row.names = rownames(df), group = df[,paste0(i,'_group')], Sample = rownames(df))
group = group[group$group != 'no',]

# gsva矩阵差异分析
library(limma)
library(stringr)
# 分组矩阵
Group = data.frame(row.names = rownames(group), group = group$group)
table(Group$group)
Group = as.matrix(Group)
class(Group)

design = model.matrix(~0+factor(Group))
colnames(design) = levels(factor(Group))
rownames(design) = rownames(Group)
head(design)

# 比较矩阵
contrast.matrix<-makeContrasts("High-Low",
                               levels = design)
contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
# Contrasts
# Levels High-Low
# High        1
# Low        -1

# 差异分析函数
deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

# KEGG 差异分析
KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
KEGG_logFC_list[[i]] = KEGG_logFC

# GO 差异分析
GO_logFC = deg(go_gsva_df,design,contrast.matrix)
GO_logFC_list[[i]] = GO_logFC


}

save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,all_kegg_gsva_list,all_go_gsva_list,
     KEGG_logFC_list,GO_logFC_list,file = 'Step6-1 Result.Rdata')


# 作图
load(file = 'Step6-1 Result.Rdata')
library(dplyr)
library(tidyr)
library(tibble)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(stringr)
# MS4A7
# KEGG
KEGG_logFC = KEGG_logFC_list[["MS4A7"]]
kk_gsva_df = as.data.frame(all_kegg_gsva_list[["MS4A7"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("MS4A7",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in MS4A7 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'MS4A7 kegg gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


# GO
go_re = go_list[["MS4A7"]]@result
go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[["MS4A7"]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.65)

go_gsva_df = as.data.frame(all_go_gsva_list[["MS4A7"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in MS4A7 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'MS4A7 go gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')




# RETN
# KEGG
KEGG_logFC = KEGG_logFC_list[["RETN"]]
kk_gsva_df = as.data.frame(all_kegg_gsva_list[["RETN"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("RETN",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in RETN Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'RETN kegg gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


# GO
go_re = go_list[["RETN"]]@result
go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[["RETN"]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.1)

go_gsva_df = as.data.frame(all_go_gsva_list[["RETN"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in RETN Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'RETN go gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')



# CXCR2
# KEGG
KEGG_logFC = KEGG_logFC_list[["CXCR2"]]
kk_gsva_df = as.data.frame(all_kegg_gsva_list[["CXCR2"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("CXCR2",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in CXCR2 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'CXCR2 kegg gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


# GO
go_re = go_list[["CXCR2"]]@result
go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[["CXCR2"]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.25)

go_gsva_df = as.data.frame(all_go_gsva_list[["CXCR2"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in CXCR2 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'CXCR2 go gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')




# CD177
# KEGG
KEGG_logFC = KEGG_logFC_list[["CD177"]]
kk_gsva_df = as.data.frame(all_kegg_gsva_list[["CD177"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("CD177",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in CD177 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'CD177 kegg gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


# GO
go_re = go_list[["CD177"]]@result
go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[["CD177"]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.12)

go_gsva_df = as.data.frame(all_go_gsva_list[["CD177"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in CD177 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'CD177 go gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')



# CSRNP1
# KEGG
KEGG_logFC = KEGG_logFC_list[["CSRNP1"]]
kk_gsva_df = as.data.frame(all_kegg_gsva_list[["CSRNP1"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("CSRNP1",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in CSRNP1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'CSRNP1 kegg gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


# GO
go_re = go_list[["CSRNP1"]]@result
go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[["CSRNP1"]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05)

go_gsva_df = as.data.frame(all_go_gsva_list[["CSRNP1"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in CSRNP1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'CSRNP1 go gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')



# LUCAT1
# KEGG
KEGG_logFC = KEGG_logFC_list[["LUCAT1"]]
kk_gsva_df = as.data.frame(all_kegg_gsva_list[["LUCAT1"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("LUCAT1",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in LUCAT1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'LUCAT1 kegg gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


# GO
go_re = go_list[["LUCAT1"]]@result
go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[["LUCAT1"]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.1)
which(rownames(GO_logFC1) == 'oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen')
GO_logFC1 = GO_logFC1[-15,]

go_gsva_df = as.data.frame(all_go_gsva_list[["LUCAT1"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in LUCAT1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
ggsave(filename = 'LUCAT1 go gsva.pdf',width = 10,height = 8,path = '../../Fig/Step6-1/')


