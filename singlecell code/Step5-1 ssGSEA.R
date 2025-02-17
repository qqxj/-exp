

### Step5-1  探索预后模型的免疫浸润
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step5-1/')
rm(list = ls())
options(stringsAsFactors = F)


# 获取28种免疫细胞基因集
library(dplyr)
library(tidyverse)
geneSet <- read.csv("CellReports.txt",header = F,sep = "\t",) # 用EXCEL打开删除NA列
class(geneSet)
geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

Imu = as.list(c(l[["Central memory CD8 T cell"]],l[["Effector memeory CD8 T cell"]],l[["Activated CD4 T cell"]],
                l[["Central memory CD4 T cell"]],l[["Effector memeory CD4 T cell"]],l[["T follicular helper cell"]],
                l[["Gamma delta T cell"]],l[["Type 1 T helper cell"]],l[["Type 17 T helper cell"]],l[["Type 2 T helper cell"]],
                l[["Regulatory T cell"]],l[["Activated B cell"]],l[["Immature  B cell"]],l[["Memory B cell"]],l[["Natural killer cell"]],
                l[["CD56bright natural killer cell"]],l[["CD56dim natural killer cell"]],l[["Myeloid derived suppressor cell"]],
                l[["Natural killer T cell"]],l[["Activated dendritic cell"]],l[["Plasmacytoid dendritic cell"]],l[["Immature dendritic cell"]],
                l[["Macrophage"]],l[["Eosinophil"]],l[["Mast cell"]],l[["Monocyte"]],l[["Neutrophil"]]))

Imu = c(l[["Central memory CD8 T cell"]],l[["Effector memeory CD8 T cell"]],l[["Activated CD4 T cell"]])



# 除以4个内参基因的表达矩阵
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
refGene <- c("ACTB","GAPDH","TFRC","TUBB")
refExpr <- as.data.frame(t(NSCLCcount[refGene,]))
Expr <- NSCLCcount[refGene,]
Expr <- as.data.frame(t(NSCLCcount[-refGene,]))
Basal_CT <- as.data.frame(t(NSCLCcount))
gene <- colnames(Basal_CT)
# Basal_CT[,gene[1]]
for (i in gene) {
  Basal_CT[,i] = ((Basal_CT[,i] - Basal_CT$ACTB)+(Basal_CT[,i]-Basal_CT$GAPDH)+(Basal_CT[,i]-Basal_CT$TFRC)+(Basal_CT[,i]-Basal_CT$TUBB))/4
}
refGene <- c("ACTB","GAPDH","TFRC","TUBB")
for (i in 1:length(refGene)) {
  p = which(colnames(Basal_CT) == refGene[i])
  print(paste0(refGene[i],' in ',p))
}
# [1] "ACTB in 20765"
# [1] "GAPDH in 33795"
# [1] "TFRC in 12195"
# [1] "TUBB in 18467"
df <- as.data.frame(t(Basal_CT[,c(-20765,-33795,-12195,-18467)]))
df <- scale(df)
df <- as.data.frame(df)
# 导出数据后，需要去除第一行第一列的空值，并改为.txt格式
write.csv(df,file = "TCGA_scale.csv")


# 运行ssGSEA
# BiocManager::install("GSVA")
# BiocManager::install("limma")
library(GSVA)
library(limma)
# exprSet <- as.matrix(read.table(file = "TCGA_scale.txt"))#基因的表达量需要是矩阵，行为基因，列为样本
# exprSet[1:6,1:6]
# #开始进行ssGSEA
# ssgsea<- gsva(exprSet, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# #基因集需要是list为对象。默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。当输入表达式值是整数计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"
# Immune = as.data.frame(t(ssgsea))
# Immune$id <- rownames(Immune)
# save(Immune,file = 'Step7 scale_ssGSEA.Rdata')


# load(file = 'Step1 output.Rdata')
exprSet <- as.matrix(NSCLCcount)#基因的表达量需要是矩阵，行为基因，列为样本
exprSet[1:6,1:6]
#开始进行ssGSEA
ssgsea<- gsva(exprSet, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#基因集需要是list为对象。默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。当输入表达式值是整数计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"
Immune2 = as.data.frame(t(ssgsea))
Immune2$id <- rownames(Immune2)
save(Immune2,file = 'FPKM_ssGSEA.Rdata')



# 加载训练集
load(file = "../Step4-1/Train_Result.Rdata") 
df <- read.csv(file = 'TCGA_scale.csv')
rownames(df) <- df$X
df <- df[,-1]
df[1:6,1:6]
exprSet <- as.data.frame(t(df[rownames(gene_multi),]))
dim(exprSet)
# exprSet <- as.data.frame(t(df[rownames(gene_multi),rownames(sv2)]))
# dim(exprSet)
exprSet$RiskScore <- exprSet[,rownames(gene_multi)[1]]*gene_multi[1,1]+exprSet[,rownames(gene_multi)[2]]*gene_multi[2,1]+exprSet[,rownames(gene_multi)[3]]*gene_multi[3,1]+exprSet[,rownames(gene_multi)[4]]*gene_multi[4,1]+exprSet[,rownames(gene_multi)[5]]*gene_multi[5,1]+exprSet[,rownames(gene_multi)[6]]*gene_multi[6,1]
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore) , "Low","High")
table(exprSet$RiskGroup)
exprSet$id <- rownames(exprSet)
table(substr(rownames(exprSet),14,15))
exprSet$Type <- ifelse(substr(rownames(exprSet),14,15) == 11, 'Normal','Tumor')
table(exprSet$Type)
exprSet <- subset(exprSet,exprSet$Type == 'Tumor')


ImmSet <- merge(Immune2,exprSet,by = "id")
rownames(ImmSet) <- ImmSet$id
ImmSet2 <- ImmSet[,2:29]
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
dim(ImmSet2)
dim(ImmSet3)
ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == ImmSet$id ,ImmSet$RiskGroup,0)
table(ImmSet3$RiskGroup)
table(ImmSet$RiskGroup)
save(ImmSet,ImmSet3,Immune2,file = 'ssGSEA_output.Rdata')


# 小提琴图
library(ggplot2)
library(tidyverse)
library(ggpubr)
# plot
colnames(ImmSet3)
ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  # 小提琴图层
  geom_violin(position = position_dodge(0.9),alpha = 1.2,
              width = 1.8,trim = T,
              color = NA) +
  # 箱线图图层
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  # outlier.shape = 21,
  # outlier.shape = NA
  # 主题调整
  theme_bw(base_size = 16) +
  labs(x = "Immune Cell Type", y = 'ssGSEA Estimating Score') +
  theme(axis.text.x = element_text(angle = 65,hjust = 1,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  # 颜色设置
  # scale_fill_manual(values = c('Low'='#398AB9','High'='red'),
  #                   name = '') +
  scale_fill_manual(values = c("#FF0000",'#00CCFF'))+ 
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup,label = ..p.signif..),method = "wilcox.test")#kruskal.test  p.signif
  # 添加显著性标记
  # stat_compare_means(aes(group=RiskGroup),
  #                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
  #                                    
  #                                     symbols = c("***", "**", "*", "NS")),label = "p.signif",
  #                    label.y = 1.7,size = 4) + ylim(0.3,1.7)
ggsave(filename = '../../Fig/Step5-1/ssGSEA.pdf',width = 16,height = 10)



# 相关性分析
library(ggplot2)
library(ggstatsplot)
library(cowplot)

Immcell <- colnames(ImmSet)[2:29]
# 检查是否符合正态分布
for (i in 1:length(Immcell)) {
  x=i+1
a <- shapiro.test(ImmSet[,x])
ifelse(a[["p.value"]] > 0.05,
       print(paste0(Immcell[i],' is Yes')),
       print(paste0(Immcell[i],' is No')))
}


Immcell <- colnames(ImmSet)[2:29]
plotlist <- list()
library(rlang)
for (i in 1:length(Immcell)) {
  print(paste0('Now is ',Immcell[i]))
  xlab1 <- sym(Immcell[i])
  p1 <- ggscatterstats(data = ImmSet,
                       y = RiskScore,
                       x = !!xlab1,
                       # centrality.para = "mean",
                       # margins = "both",
                       bf.message = FALSE,#去除贝叶斯相关的统计值
                       type = "nonparamatric",#选择非参数检验
                       xfill = "#CC79A7",
                       yfill = "#009E73",
                       # marginal.type = "density", # 其中marginal.type可选 histograms，boxplots，density，violin，densigram (density + histogram)
                       title = paste0("Relationship between RiskScore and ",Immcell[i]))
  # print(p1)
  # ggsave(filename = paste0(Immcell[i],'.png'), plot = p1, width = 8, height = 6, path = '../Figer/Relationship')
  plotlist[[i]] = p1
}

plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],
          plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],nrow = 2)

save(ImmSet,ImmSet3,Immune2,plotlist,file = 'ssGSEA_output.Rdata')
# p1 <- ggscatterstats(data = ImmSet,
#                      y = RiskScore,
#                      x = 'Activated CD8 T cell',
#                      # centrality.para = "mean",
#                      # margins = "both",
#                      bf.message = FALSE,#去除贝叶斯相关的统计值
#                      type = "nonparamatric",#选择非参数检验
#                      xfill = "#CC79A7",
#                      yfill = "#009E73",
#                      # marginal.type = "density", # 其中marginal.type可选 histograms，boxplots，density，violin，densigram (density + histogram)
#                      title = "Relationship between RiskScore and Activated CD8 T cell")
# p1


#6.28 
load(file = 'ssGSEA_output.Rdata')
# 小提琴图
library(ggplot2)
library(tidyverse)
library(ggpubr)
# plot
colnames(ImmSet3)
ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  # 小提琴图层
  geom_violin(position = position_dodge(0.9),alpha = 1.2,
              width = 1.8,trim = T,
              color = NA) +
  # 箱线图图层
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  # outlier.shape = 21,
  # outlier.shape = NA
  # 主题调整
  theme_bw(base_size = 16) +
  labs(x = "Immune Cell Type", y = 'ssGSEA Estimating Score') +
  theme(axis.text.x = element_text(angle = 65,hjust = 1,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  # 颜色设置
  # scale_fill_manual(values = c('Low'='#398AB9','High'='red'),
  #                   name = '') +
  scale_fill_manual(values = c("#FF0000",'#00CCFF'))+ 
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup,label = ..p.signif..),method = "wilcox.test")#kruskal.test  p.signif
# 添加显著性标记
# stat_compare_means(aes(group=RiskGroup),
#                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
#                                    
#                                     symbols = c("***", "**", "*", "NS")),label = "p.signif",
#                    label.y = 1.7,size = 4) + ylim(0.3,1.7)
ggsave(filename = '../../Fig/Step5-1/ssGSEA.pdf',width = 16,height = 10)


# 计算相关性
Cor_imm = ImmSet[,c(2:29,36)]
Cor_imm = Cor_imm[,c(29,1:28)]
Cor_imm1 = cor(Cor_imm, method = "spearman")
head(Cor_imm1)
# 计算显著性
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(Cor_imm1)
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_risk.pdf',width = 14,height = 12)
corrplot(Cor_imm1, type="upper",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between Risk Scores and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
         #p.mat = p.mat, sig.level = 0.05 , insig = "blank")  #添加显著性

dev.off()


# 6gene immune
load(file = 'ssGSEA_output.Rdata')
cor_df = ImmSet[,1:35]
rownames(cor_df) = cor_df$id
cor_df = cor_df[,-1]
cor_df = cor_df[,c(29:34,1:28)]
# 计算相关性
cor_imm = cor(cor_df, method = "spearman")
head(cor_imm)
# 计算显著性
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(cor_imm)

# MS4A7
cor_imm1 = cor_imm[-c(2:6),-c(2:6)]
p.mat1 = p.mat[-c(2:6),-c(2:6)]
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_MS4A7.pdf',width = 14,height = 12)
corrplot(cor_imm1, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between MS4A7 Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
        # p.mat = p.mat1, sig.level = 0.05 , insig = "blank")  #添加显著性
dev.off()


# RETN
cor_imm2 = cor_imm[-c(1,3:6),-c(1,3:6)]
p.mat2 = p.mat[-c(1,3:6),-c(1,3:6)]
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_RETN.pdf',width = 14,height = 12)
corrplot(cor_imm2, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between RETN Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
         #p.mat = p.mat2, sig.level = 0.05 , insig = "blank")  #添加显著性
dev.off()


# CXCR2
cor_imm3 = cor_imm[-c(1,2,4:6),-c(1,2,4:6)]
p.mat3 = p.mat[-c(1,2,4:6),-c(1,2,4:6)]
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_CXCR2.pdf',width = 14,height = 12)
corrplot(cor_imm3, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between CXCR2 Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
         # p.mat = p.mat3, sig.level = 0.05 , insig = "blank")  #添加显著性
dev.off()


# CD177
cor_imm4 = cor_imm[-c(1:3,5:6),-c(1:3,5:6)]
p.mat4 = p.mat[-c(1:3,5:6),-c(1:3,5:6)]
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_CD177.pdf',width = 14,height = 12)
corrplot(cor_imm4, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between CD177 Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F,)
         #p.mat = p.mat4, sig.level = 0.05 , insig = "blank")  #添加显著性
dev.off()


# CSRNP1
cor_imm5 = cor_imm[-c(1:4,6),-c(1:4,6)]
p.mat5 = p.mat[-c(1:4,6),-c(1:4,6)]
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_CSRNP1.pdf',width = 14,height = 12)
corrplot(cor_imm5, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between CSRNP1 Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
        #  p.mat = p.mat5, sig.level = 0.05 , insig = "blank")  #添加显著性
dev.off()

# LUCAT1
cor_imm6 = cor_imm[-c(1:5),-c(1:5)]
p.mat6 = p.mat[-c(1:5),-c(1:5)]
library(corrplot)
pdf(file = '../../Fig/Step5-1/28Imm_LUCAT1.pdf',width = 14,height = 12)
corrplot(cor_imm6, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between LUCAT1 Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
        # p.mat = p.mat6, sig.level = 0.05 , insig = "blank")  #添加显著性

dev.off()

