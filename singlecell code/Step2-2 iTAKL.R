

### Step2-2  Single-Cell Analyze：iTAKL
### 整理时间： 2022/7/22
### 作者： 庞建宇

rm(list =ls())
options(stringsAsFactors = F)

setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step2-2/")

library(Seurat)
library(ggplot2)
library(monocle)

## 更新所有细胞信息
load(file = "../Step1-2/res0.4_TBcell.Rdata")
TBcell@meta.data$New <- TBcell@active.ident
pd1 <- new('AnnotatedDataFrame', data = TBcell@meta.data)
TBCell <- as.data.frame(pd1@data)
rm(TBcell)
gc()
df1 <- TBCell[,c(1:4,23)]
table(df1$New)



load(file = "../Step1-3/res0.4_MonDC.Rdata")
MonDC@meta.data$New <- MonDC@active.ident
pd2 <- new('AnnotatedDataFrame', data = MonDC@meta.data)
MonDCset <- as.data.frame(pd2@data)
rm(MonDC)
gc()
df2 <- MonDCset[,c(1:4,23)]
table(df2$New)



load(file = "../Step1-4/res0.4_CancerCell.Rdata")
pd3 <- new('AnnotatedDataFrame', data = CancerCell@meta.data)
CancerCellSet <- as.data.frame(pd3@data)
CancerCellSet$New <- rep('Cancer Cell',length(rownames(CancerCellSet)))
rm(CancerCell)
gc()
df3 <- CancerCellSet[,c(1:4,25)]
table(df3$New)



load(file = "../Step1-4/res0.4_Normal_Epi.Rdata")
Normal_Epi@meta.data$New <- Normal_Epi@active.ident
pd4 <- new('AnnotatedDataFrame', data = Normal_Epi@meta.data)
NormalEpiSet <- as.data.frame(pd4@data)
rm(Normal_Epi)
gc()
df4 <- NormalEpiSet[,c(1:4,25)]
table(df4$New) # Unkonw 773


AllSet <- rbind(df1,df2,df3,df4)


load(file = '../Step1-1/res0.4_NSCLC.Integrate.Rdata')
# 提取Cell注释信息
NSCLC.Integrate@meta.data$New <- NSCLC.Integrate@active.ident
pd <- new('AnnotatedDataFrame', data = NSCLC.Integrate@meta.data)
NSCLCinter <- as.data.frame(pd@data)

table(NSCLCinter$New)

NSCLCinter2 <- NSCLCinter[,c(1:4,23)]
colnames(NSCLCinter2)[5] <- 'old'


library(dplyr)
NSCLCinter2$id <- rownames(NSCLCinter2)
AllSet$id <- rownames(AllSet)

AllSet2 <- full_join(AllSet,NSCLCinter2,by = 'id')
rownames(AllSet2) <- AllSet2$id

str(AllSet2)
table(AllSet2$New)
table(is.na(AllSet2$New))

AllSet3 <- AllSet2[,c(5:6,11)]
table(is.na(AllSet3))

AllSet3[is.na(AllSet3)] <- 'Unkonw' # 替换NA

str(AllSet3)
AllSet3$New <- as.character(AllSet3$New)
table(AllSet3$New)
table(AllSet$New)
AllSet3$old <- as.character(AllSet3$old)

AllSet3$Rename <- ifelse(AllSet3$New == 'Unkonw', AllSet3$old,  AllSet3$New)
table(AllSet3$Rename)
table(AllSet$New)

AllSet3$Rename <- ifelse(AllSet3$Rename == 'Epithelial/Cancer Cell', 'Unkonw',  AllSet3$Rename)
table(AllSet3$Rename)

AllSet3$Rename <- as.factor(AllSet3$Rename)
table(AllSet3$Rename)
str(AllSet3$Rename)


# 添加新的细胞信息到Seurat对象中
NSCLC.Integrate <- AddMetaData(NSCLC.Integrate, metadata = AllSet3)
Idents(object = NSCLC.Integrate) <- "Rename"
head(Idents(NSCLC.Integrate), 5)#查看前5个细胞的分类ID
table(Idents(NSCLC.Integrate))

DimPlot(NSCLC.Integrate, reduction = 'umap',label = F)

save(NSCLC.Integrate,file = '../Step1-1/res0.4_NSCLC.Integrate.Rdata')



#重新分配病人分组
# 1 GSE131907
library(dplyr)
library(stringr)

scRNA1 <- readRDS("/home/pjy/NSCLC/NSCLC/Data/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
scRNA1[1:6,1:6]


# substring 可以只设置first参数，last参数缺省时则默认为1000000，指字符串的最大长度。


phe1 <- read.csv(file = '/home/pjy/NSCLC/NSCLC/Data/GSE131907/phe.csv')
df <- data.frame(row.names = colnames(scRNA1), patient = ifelse(substring(colnames(scRNA1),18) %in% phe1$Samples ,phe1$Patient.id, phe1$Patient.id),id=colnames(scRNA1))
table(df$patient)
unique(df$patient)
# df$id = rownames(df)
rm(scRNA1)
gc()


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


d = list()
for (i in 1:42) {
  b = as.data.frame(sceList_NSCLC7[[i]]@assays[["RNA"]]@counts@Dimnames[2])
  c = data.frame(row.names = b[,1], patient = rep(paste0('P',i),sceList_NSCLC7[[i]]@assays[["RNA"]]@counts@Dim[2]), id = b[,1])
  d[[i]] = c
}

df2 <- do.call(rbind,d) # 多个数据集合并
df3 <- rbind(df,df2)


# 提取细胞信息
e = as.data.frame(NSCLC.Integrate@assays[["RNA"]]@counts@Dimnames[[2]])
# 更改细胞名称与原始细胞名称一致
f = as.data.frame(e[188767:length(e[,1]),])
o = as.data.frame(f[1:6483,])
o$name <- substring(o[,1],5)
names(o)[1] <- 'oldid'

o1 = as.data.frame(f[6484:length(f[,1]),])
o1$name <- substring(o1[,1],6)
names(o1)[1] <- 'oldid'
f1 = rbind(o,o1)

f2 = as.data.frame(e[1:188766,])
names(f2) <- 'oldid'
f2$name <- f2$oldid

f3 <- rbind(f1,f2)

# 添加病人信息
df4 <- df3[f3$name,]
df4$name <- rownames(df4)
df5 <- merge(df4,f3,by = 'name')
rownames(df5) <- df5$oldid
df6 <- df5[,2:4]


# 添加病人信息到Seurat对象中
NSCLC.Integrate <- AddMetaData(NSCLC.Integrate, metadata = df6)
# Idents(object = NSCLC.Integrate) <- "Rename"
# head(Idents(NSCLC.Integrate), 5)#查看前5个细胞的分类ID
# table(Idents(NSCLC.Integrate))

DimPlot(NSCLC.Integrate, reduction = 'umap',group.by = 'patient',label = F,raster=T)
ggsave(filename = 'patient.pdf',width = 14,height = 8,path = '../../Fig/Step1-1/')

DimPlot(NSCLC.Integrate, reduction = 'umap',label = F)

save(NSCLC.Integrate,file = '../Step1-1/res0.4_NSCLC.Integrate.Rdata')




load(file = '../Step1-1/res0.4_NSCLC.Integrate.Rdata')
table(Idents(NSCLC.Integrate))
table(NSCLC.Integrate@meta.data$Rename)

# #27 cell
# datSet  = NSCLC.Integrate[,!(NSCLC.Integrate@meta.data$Rename %in% c('Oligodendrocytes','Unkonw'))]
# table(Idents(datSet))
# table(datSet@meta.data$Rename)
# 
# save(datSet,file = './Rdata/Step2-2/res0.4_NSCLC.Integrate_27cell.Rdata')

#28 cell
datSet  = NSCLC.Integrate[,!(NSCLC.Integrate@meta.data$Rename %in% 'Unkonw')]
table(Idents(datSet))
table(datSet@meta.data$Rename)

save(datSet,file = '../Step2-2/res0.4_NSCLC.Integrate_28cell.Rdata')

rm(NSCLC.Integrate)
gc()


# load(file = 'res0.4_NSCLC.Integrate_27cell.Rdata')
# Marker <- FindAllMarkers(datSet, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5) # 运行时间较长
# Marker_top10 <-Marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# save(Marker,Marker_top10,file = 'res0.4_NSCLC.Integrate_27cell_Marker.Rdata')
# 
# 
# load(file = 'res0.4_NSCLC.Integrate_27cell_Marker.Rdata')
# load(file = "Genedata.Rdata")
# Marker_data <-  Genedata[Marker_top10$gene,]
# Marker_data$symble <- rownames(Marker_data)
# Marker_data <- Marker_data[,c(28,1:27)]
# rownames(Marker_data) <- NULL
# write.table(Marker_data,file="top10MarGenedata.txt",sep = "\t")# EXCEL打开，删除第一列


# CD177 CSRNP1 CXCR2 LUCAT1 MS4A7 RETN
# install.packages("remotes")
# remotes::install_github("lyc-1995/MySeuratWrappers")
# library(MySeuratWrappers)

# plot_markers <- sort(c('CD177', 'CSRNP1', 'CXCR2', 'LUCAT1', 'MS4A7', 'RETN'))
# VlnPlot(datSet, features = plot_markers,stacked=T,pt.size=0.5)
# 
# 
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
# 
# #配色方案
# my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
#                 '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
#                 '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
#                 '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
#                 '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
#                 '#968175')
# StackedVlnPlot(datSet, plot_markers, pt.size=0, cols=my36colors) 



# load(file = './Rdata/Step2-2/res0.4_NSCLC.Integrate_28cell.Rdata')
# plot_markers <- sort(c('CD177', 'CSRNP1', 'CXCR2', 'LUCAT1', 'MS4A7', 'RETN'))
# library(reshape2)
# library(ggsci)
# library(dplyr)
# vln.df=as.data.frame(datSet[["RNA"]]@data[plot_markers,])
# vln.df[1:6,1:6]
# 
# vln.df$gene=rownames(vln.df)
# vln.df=melt(vln.df ,id="gene") # 转化为包含gene cell exp 三列的数据框
# colnames(vln.df)[c(2,3)]= c("CB", "exp")
# 
# head(vln.df)
# vln.df[1:10,]
# 
# anno=data.frame(Cell = datSet@active.ident)
# anno$CB <- rownames(anno)
# 
# vln.df=inner_join(vln.df ,anno,by="CB")
# vln.df$gene=factor(vln.df$gene,levels = plot_markers)#为了控制画图的基因顺序
# 
# vln.df %>% ggplot(aes(Cell,exp)) + geom_violin(aes(fill=gene),scale = "width") + 
#   facet_grid(vln.df$gene~. ,scales = "free_y")+
#   # scale_fill_brewer(palette = "set3",direction = 1) +
#   scale_x_discrete("") + scale_y_continuous("")+
#  theme_bw()+
#  labs(title='scRNA')+ # 添加标题，x轴，y轴内容
#  theme(
#    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1, face = "bold"),
#    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#    legend.position = "none", plot.title=element_text(hjust=0.5,face = "bold.italic",size = 14))+
#   scale_fill_lancet()
# 
# 




# iTAKL

library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
# load(file = 'res0.4_NSCLC.Integrate_27cell.Rdata')

load(file = '../Step2-2/res0.4_NSCLC.Integrate_28cell.Rdata')
# sdata <- readRDS(file = "~/Desktop/Seurat.rds")
# sdata <- readRDS(file = "iTALK/Seurat.rds")
# iTALK 要求的矩阵: 行为细胞，列为基因
# datSet@assays$RNA@counts
# iTalk_data2 <- as.matrix(datSet@assays$RNA@data)
iTalk_data2 <- datSet@assays$RNA@data
iTalk_data2 <- as.data.frame(t(iTalk_data2)) # 报错，稀疏矩阵太大，需要转换
# Error in asMethod(object) : 
#   Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102

# save(iTalk_data2,file = '27cell_countdata.Rdata')
# rm(list = ls())
# gc()
# # .rs.restartR()
# 
# 
# load(file = '27cell_countdata.Rdata')
# dim(iTalk_data2)
# # [1]  32660 197932
# str(iTalk_data2)
# options(max.print=1000000) 
# # 将稀疏矩阵数据重新录入matrix
# as_matrix <- function(mat){
#   
#   tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
#   
#   row_pos <- mat@i+1
#   col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
#   val <- mat@x
#   
#   for (i in seq_along(val)){
#     tmp[row_pos[i],col_pos[i]] <- val[i]
#   }
#   
#   row.names(tmp) <- mat@Dimnames[[1]]
#   colnames(tmp) <- mat@Dimnames[[2]]
#   return(tmp)
# } 
# 
# iTalk_data0 <- as_matrix(iTalk_data2)

# save(iTalk_data0,iTalk_data,file = '27cell_countdata.Rdata')


# write.table(iTalk_data0, 'cellphonedb_count.txt', sep='\t', quote=F)



# 28cell
dim(iTalk_data2)
# [1] 32660 198603
str(iTalk_data2)
options(max.print=1000000) 
# 将稀疏矩阵数据重新录入matrix
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
} 
iTalk_data0 <- as_matrix(iTalk_data2)

# save(iTalk_data0,iTalk_data2,file = '27cell_countdata.Rdata')
# write.table(iTalk_data0, 'cellphonedb_28cell_count.txt', sep='\t', quote=F)

meta_data <- cbind(rownames(datSet@meta.data), datSet@meta.data[,'Rename', drop=F])  
meta_data <- as.matrix(meta_data)
table(is.na(meta_data))
# meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
# write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
# write.table(meta_data, 'cellphonedb_28cell.txt', sep='\t', quote=F, row.names=F)

table(datSet@meta.data$patient)
meta_data0 <- cbind(rownames(datSet@meta.data), datSet@meta.data[,'patient', drop=F])  
# meta_data0 <- as.matrix(meta_data0)

library(dplyr)
meta_data1 <- inner_join(meta_data0,meta_data,by = 'rownames(datSet@meta.data)')
meta_data1[1:6,]
rownames(meta_data1) <- meta_data1[,1]
meta_data1 <- meta_data1[,-1]

# save(meta_data1,file = './Rdata/Step2-2/iTAKL.Rdata')
save(meta_data1,iTalk_data2,file = 'iTAKL_28cell.Rdata')
# iTalk_data2 <- read.delim(file = 'cellphonedb_count.txt')
iTalk_data0 <- iTalk_data2
iTalk_data0[1:6,1:6]
iTalk_data1 <- t(iTalk_data0)
iTalk_data1[1:6,1:6]

rm(iTalk_data0)
gc()
# save(iTalk_data0,iTalk_data,file = '27cell_countdata.Rdata')
# save(iTalk_data1,meta_data1,file = 'iTAKL.Rdata')
# 
# save(iTalk_data1,meta_data1,file = 'iTAKL_28cell.Rdata')

# load(file = 'iTAKL.Rdata')
# load(file = '27cell_countdata.Rdata')

iTalk_data2 <- as.data.frame(iTalk_data1)
iTalk_data2 <- cbind(iTalk_data2,meta_data1)
iTalk_data2[1:6,32660:32662]
# iTalk_data2 <- iTalk_data2[,-32661]

names(iTalk_data2)[32661:32662] <- c('compare_group','cell_type')
iTalk_data2[1:6,32660:32662]

# iTalk_data3 <- as.matrix(iTalk_data2)
# save(iTalk_data2,file = 'iTAKL2.Rdata')
save(iTalk_data2,file = 'iTAKL_28cell.Rdata')

rm(meta_data1)
gc()




# # iTALK 要求包含cell_type列
# table(datSet@active.ident)
# iTalk_data1$cell_type <- datSet@active.ident
# 
# 
# 
# # iTALK 要求包含compare_group列（多样本），表示每个细胞的生物学分组/样本
# iTalk_data1$compare_group <- datSet@meta.data$Group

unique(iTalk_data2$cell_type)
unique(iTalk_data2$compare_group)


# 配体-受体概览
# 通过所有细胞的高表达基因分析其中包含的配体-受体
library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggsci)
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step2-2/')

# load(file = 'iTAKL2.Rdata')
load(file = 'iTAKL_28cell.Rdata')

mycolors <- pal_igv("default")(51)

# my10colors <- my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282')

highly_exprs_genes <- rawParse(iTalk_data2, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data2$cell_type)
cell_col <- structure(mycolors[1:length(cell_types)], names=cell_types)

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}

iTalk_res2 <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:20,] 

# save(iTalk_res,iTalk_res2,file = 'iTAKL_result.Rdata')

save(iTalk_data2,iTalk_res,iTalk_res2,highly_exprs_genes,file = 'iTAKL_28cell_result.Rdata')



# 结果可视化
# load(file = 'iTAKL_result.Rdata')

load(file = 'iTAKL_28cell_result.Rdata')
table(iTalk_res$comm_type)

iTalk_res2 <- iTalk_res
name <- unique(iTalk_res$cell_from)
table(iTalk_res$cell_from)
table(iTalk_res$cell_to)

# ID太长 影响视觉效果 做简单转换
for (i in 1:length(name)) {
  iTalk_res2$cell_from[iTalk_res2$cell_from == name[i]] = paste0('C',i)
  iTalk_res2$cell_to[iTalk_res2$cell_to == name[i]] = paste0('C',i)
}  

table(iTalk_res2$cell_from)
table(iTalk_res2$cell_to)
unique(iTalk_res2$cell_from)
rename = data.frame(beforename = name, aftername = unique(iTalk_res2$cell_from))


df <- iTalk_res %>% filter(cell_from_mean_exprs >= 0.1 & cell_to_mean_exprs >= 0.1) %>% sort(cell_from_mean_exprs,decreasing = F)
table(df$comm_type)
unique(df$cell_from)
rownames(df) = NULL
write.csv(df,file = 'All_0.1_LR.csv')





# checkpoint
checkpoint <- subset(iTalk_res2,iTalk_res2$comm_type == 'checkpoint')

checkpoint2 <- subset(df,df$comm_type == 'checkpoint')
checkpoint2 <- checkpoint2[order(checkpoint2$cell_from_mean_exprs*checkpoint2$cell_to_mean_exprs,decreasing=T),][1:30,]
rownames(checkpoint2) <- NULL


# 如果有超过20组配体-受体结果，取前30进行展示
res0 <- checkpoint[order(checkpoint$cell_from_mean_exprs*checkpoint$cell_to_mean_exprs,decreasing=T),][1:30,]
table(checkpoint2$cell_from_mean_exprs %in% res0$cell_from_mean_exprs)
# TRUE 
# 30 


mycolors0 <- pal_igv()(28)
cell_types <- unique(iTalk_res$cell_from)
cell_col <- structure(mycolors0[1:length(cell_types)], names=cell_types)
NetView(checkpoint,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)


# mycolors2 <- pal_npg("nrc")(10)
mycolors <- pal_ucscgb()(26)
mycolors2 <- pal_d3('category20')(20)

genes <- unique(c(unique(res0$ligand),unique(res0$receptor)))
genes_col <- structure(mycolors2[1:length(genes)], names=genes)

cell = unique(c(unique(res0$cell_from),unique(res0$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

arr_col = c(rep("#1F77B4FF",2),"#FF7F0EFF",rep("#1F77B4FF",2),rep("#2CA02CFF",4),'#1F77B4FF',
            rep("#2CA02CFF",2),"#1F77B4FF","#2CA02CFF","#1F77B4FF","#2CA02CFF", "#D62728FF",
            rep("#2CA02CFF",2), rep("#1F77B4FF",2), rep("#2CA02CFF",3), "#1F77B4FF","#2CA02CFF",
            "#1F77B4FF", "#D62728FF", "#FF7F0EFF", "#1F77B4FF")

LRPlot(res0,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res0$cell_from_mean_exprs,
       link.arr.width=res0$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res0,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)




# cytokine
cytokine <- subset(iTalk_res2,iTalk_res2$comm_type == 'cytokine')

cytokine2 <- subset(df,df$comm_type == 'cytokine')
cytokine2 <- cytokine2[order(cytokine2$cell_from_mean_exprs*cytokine2$cell_to_mean_exprs,decreasing=T),][1:30,]
rownames(cytokine2) <- NULL


# 如果有超过20组配体-受体结果，取前30进行展示
res1 <- cytokine[order(cytokine$cell_from_mean_exprs*cytokine$cell_to_mean_exprs,decreasing=T),][1:30,]
table(cytokine2$cell_from_mean_exprs %in% res1$cell_from_mean_exprs)
# TRUE 
# 30 


mycolors0 <- pal_igv()(28)
cell_types <- unique(iTalk_res$cell_from)
cell_col <- structure(mycolors0[1:length(cell_types)], names=cell_types)
NetView(cytokine,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)


genes <- unique(c(unique(res1$ligand),unique(res1$receptor)))
genes_col <- structure(mycolors2[1:length(genes)], names=genes)

cell = unique(c(unique(res1$cell_from),unique(res1$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

arr_col = c(rep("#1F77B4FF",7),"#FF7F0EFF","#1F77B4FF","#2CA02CFF",
            rep("#1F77B4FF",3),"#D62728FF",rep("#1F77B4FF",4),"#9467BDFF",rep("#1F77B4FF",11))

LRPlot(res1,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res1,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)




# growth
growth <- subset(iTalk_res2,iTalk_res2$comm_type == 'growth factor')

growth2 <- subset(df,df$comm_type == 'growth factor')
growth2 <- growth2[order(growth2$cell_from_mean_exprs*growth2$cell_to_mean_exprs,decreasing=T),][1:30,]
rownames(growth2) <- NULL


mycolors0 <- pal_igv()(28)
cell_types <- unique(iTalk_res$cell_from)
cell_col <- structure(mycolors0[1:length(cell_types)], names=cell_types)
NetView(growth,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)


# 如果有超过20组配体-受体结果，取前30进行展示
res2 <- growth[order(growth$cell_from_mean_exprs*growth$cell_to_mean_exprs,decreasing=T),][1:30,]
table(growth2$cell_from_mean_exprs %in% res2$cell_from_mean_exprs)
# TRUE 
# 30 


genes <- unique(c(unique(res2$ligand),unique(res2$receptor)))
genes_col <- structure(mycolors2[1:length(genes)], names=genes)

cell = unique(c(unique(res2$cell_from),unique(res2$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

arr_col = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#FF7F0EFF","#1F77B4FF","#D62728FF","#1F77B4FF","#FF7F0EFF",
            "#D62728FF","#2CA02CFF","#2CA02CFF","#D62728FF","#1F77B4FF","#FF7F0EFF","#1F77B4FF","#1F77B4FF",
            "#2CA02CFF","#1F77B4FF","#1F77B4FF","#FF7F0EFF","#2CA02CFF","#FF7F0EFF","#FF7F0EFF","#FF7F0EFF",
            "#D62728FF","#1F77B4FF","#D62728FF","#FF7F0EFF","#FF7F0EFF","#D62728FF")

LRPlot(res2,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res2$cell_from_mean_exprs,
       link.arr.width=res2$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res2,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)




# other
other <- subset(iTalk_res2,iTalk_res2$comm_type == 'other')

other2 <- subset(df,df$comm_type == 'other')
other2 <- other2[order(other2$cell_from_mean_exprs*other2$cell_to_mean_exprs,decreasing=T),][1:30,]
rownames(other2) <- NULL


mycolors0 <- pal_igv()(28)
cell_types <- unique(iTalk_res$cell_from)
cell_col <- structure(mycolors0[1:length(cell_types)], names=cell_types)
NetView(other,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)



# 如果有超过20组配体-受体结果，取前30进行展示
res3 <- other[order(other$cell_from_mean_exprs*other$cell_to_mean_exprs,decreasing=T),][1:30,]
table(other2$cell_from_mean_exprs %in% res3$cell_from_mean_exprs)
# TRUE 
# 30 


genes <- unique(c(unique(res3$ligand),unique(res3$receptor)))
genes_col <- structure(mycolors2[c(5,2)], names=genes)

cell = unique(c(unique(res3$cell_from),unique(res3$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

arr_col = c(rep("#9467BDFF",30))

LRPlot(res3,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res3$cell_from_mean_exprs,
       link.arr.width=res3$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res3,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)


top30 = rbind(checkpoint2,cytokine2,growth2,other2)
write.csv(top30,file = 'top30LR.csv')

save(df,rename,iTalk_res,iTalk_res2,highly_exprs_genes,top30,file = 'iTAKL_28cell_result.Rdata')
write.csv(rename,file = 'rename.csv')


load(file = 'iTAKL_28cell_result.Rdata')

df1 <- iTalk_res2 %>% filter(cell_from_mean_exprs >= 0.1 & cell_to_mean_exprs >= 0.1) %>% sort(cell_from_mean_exprs,decreasing = F)
table(df1$comm_type)
unique(df1$cell_from)
rownames(df1) = NULL

table(df1$comm_type)

# mycolors0 <- pal_igv()(28)
mycolors <- pal_ucscgb()(26)
mycolors2 <- pal_d3('category20')(20)
mycolors3 <- pal_igv("default")(51)

# cancer checkpoint
checkpoint3 <- subset(df1,df1$comm_type == 'checkpoint')
table(checkpoint3$cell_from == 'C13'| checkpoint3$cell_to == 'C13')

cancer_checkpoint = subset(checkpoint3,checkpoint3$cell_from == 'C13' | checkpoint3$cell_to == 'C13')
res4 <- cancer_checkpoint[order(cancer_checkpoint$cell_from_mean_exprs*cancer_checkpoint$cell_to_mean_exprs,decreasing=T),][1:30,]

cell = unique(c(unique(res4$cell_from),unique(res4$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

genes <- unique(c(unique(res4$ligand),unique(res4$receptor)))
genes_col <- structure(mycolors3[1:length(genes)], names=genes)

# type <- unique(c(unique(res4$comm_type),unique(res4$comm_type)))
arr_col = c("#5050FFFF",rep("#CE3D32FF",5),rep("#5050FFFF",2),rep("#749B58FF",3),"#F0E685FF","#5050FFFF",rep("#F0E685FF",2),
            rep("#749B58FF",3),rep("#F0E685FF",2),"#466983FF","#BA6338FF","#F0E685FF","#466983FF","#BA6338FF","#5050FFFF",
            rep("#749B58FF",2),"#466983FF","#BA6338FF")


LRPlot(res4,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res4$cell_from_mean_exprs,
       link.arr.width=res4$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res4,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)



# cancer cytokine
checkpoint4 <- subset(df1,df1$comm_type == 'cytokine')
table(checkpoint4$cell_from == 'C13' | checkpoint4$cell_to == 'C13')

cancer_cytokine = subset(checkpoint4,checkpoint4$cell_from == 'C13' | checkpoint4$cell_to == 'C13')
res5 <- cancer_cytokine[order(cancer_cytokine$cell_from_mean_exprs*cancer_cytokine$cell_to_mean_exprs,decreasing=T),][1:30,]

cell = unique(c(unique(res5$cell_from),unique(res5$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

genes <- unique(c(unique(res5$ligand),unique(res5$receptor)))
genes_col <- structure(mycolors3[1:length(genes)], names=genes)

arr_col = c(rep("#5050FFFF",6), "#CE3D32FF", rep("#5050FFFF",2), "#CE3D32FF", "#749B58FF", "#F0E685FF","#466983FF",
            "#BA6338FF","#5DB1DDFF","#802268FF","#6BD76BFF","#5050FFFF","#F0E685FF","#5050FFFF","#749B58FF","#5050FFFF",
            "#802268FF","#5050FFFF","#BA6338FF","#5DB1DDFF","#F0E685FF","#6BD76BFF","#F0E685FF","#749B58FF")


LRPlot(res5,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res5$cell_from_mean_exprs,
       link.arr.width=res5$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res5,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)



# cancer growth factor
checkpoint5 <- subset(df1,df1$comm_type == 'growth factor')
table(checkpoint5$cell_from == 'C13' | checkpoint5$cell_to == 'C13')

cancer_growth = subset(checkpoint5,checkpoint5$cell_from == 'C13' | checkpoint5$cell_to == 'C13')
res6 <- cancer_growth[order(cancer_growth$cell_from_mean_exprs*cancer_growth$cell_to_mean_exprs,decreasing=T),][1:30,]

cell = unique(c(unique(res6$cell_from),unique(res6$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

genes <- unique(c(unique(res6$ligand),unique(res6$receptor)))
genes_col <- structure(mycolors3[1:length(genes)], names=genes)

arr_col = c(rep("#5050FFFF",4), rep("#CE3D32FF",2), '#5050FFFF', "#CE3D32FF", "#5050FFFF", rep("#CE3D32FF",3),
            rep("#5050FFFF",3), "#CE3D32FF", "#5050FFFF", rep("#CE3D32FF",12),"#5050FFFF")


LRPlot(res6,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res6$cell_from_mean_exprs,
       link.arr.width=res6$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res6,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)




# cancer other 
checkpoint6 <- subset(df1,df1$comm_type == 'other')
table(checkpoint6$cell_from == 'C13' | checkpoint6$cell_to == 'C13')

cancer_other = subset(checkpoint6,checkpoint6$cell_from == 'C13' | checkpoint6$cell_to == 'C13')
res7 <- cancer_other[order(cancer_other$cell_from_mean_exprs*cancer_other$cell_to_mean_exprs,decreasing=T),][1:30,]

cell = unique(c(unique(res7$cell_from),unique(res7$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

genes <- unique(c(unique(res7$ligand),unique(res7$receptor)))
genes_col <- structure(mycolors3[1:length(genes)], names=genes)

arr_col = c("#5050FFFF",rep("#CE3D32FF",2),"#5050FFFF",rep("#CE3D32FF",3), "#749B58FF","#5050FFFF",rep("#F0E685FF",2),
            rep("#466983FF",2), "#F0E685FF","#466983FF","#F0E685FF","#CE3D32FF","#F0E685FF","#466983FF",rep("#F0E685FF",2),
            "#5050FFFF","#F0E685FF","#5050FFFF","#BA6338FF",rep("#466983FF",3),"#BA6338FF","#F0E685FF")


LRPlot(res7,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=res7$cell_from_mean_exprs,
       link.arr.width=res7$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)

NetView(res7,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)



cancer_top30 = rbind(res4,res5,res6,res7)
rownames(cancer_top30) = NULL
write.csv(cancer_top30,file = 'cancer_top30.csv')

cancer_all = subset(df,df$cell_from == 'Cancer Cell' | df$cell_to == 'Cancer Cell')
rownames(cancer_all) = NULL
table(cancer_all$comm_type)
write.csv(cancer_all,file = 'cancer_all.csv')


# 2022.6.8  确定最终结果图
library(ggsci)
library(iTALK)
# Excel中选择关键的受体-配体对
LR = read.csv(file = 'LR.csv')
LR = LR[,-1]
mycolors = pal_ucscgb()(23)
cell = unique(c(unique(LR$cell_from),unique(LR$cell_to)))
cell_col <- structure(mycolors[1:length(cell)], names=cell)

mycolors2 = pal_d3("category20")(20)
genes <- unique(c(unique(LR$ligand),unique(LR$receptor)))
genes_col <- structure(mycolors2[1:length(genes)], names=genes)
genes_col

mycolors3 = pal_lancet()(9)
arr_col = c(rep(mycolors3[1],13),rep(mycolors3[2],18),rep(mycolors3[3],18))

pdf(file = 'iTAKL.pdf',width = 8,height = 6)
LRPlot(LR,
       datatype='mean count',
       cell_col = cell_col,
       gene_col = genes_col,
       link.arr.col = arr_col,
       # link.arr.lty = "solid", # dashed 虚线 solid实线
       link.arr.lwd=LR$cell_from_mean_exprs,
       link.arr.width=LR$cell_to_mean_exprs,
       # link.arr.type = 'triangle', # triangle 三角形  ellipse 椭圆
)
dev.off()


