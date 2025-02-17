
##----
### Step4-3  展示预后基因的表达量、生存图
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step4-4/')
rm(list = ls())
options(stringsAsFactors = F)


# 带显著性的箱图
# 导入所需的包
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(cowplot)

# 导入数据
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
load(file = "../Step2-0/TCGA_Step1 PheData.Rdata")
load("../Step4-1/Train_Result.Rdata")

Gene <- rownames(gene_multi)
# 检查分组情况
rownames(dat) <- dat$Id
a <- dat[dat$type == "Normal",]
table(substr(rownames(a), 14,15))
expr1 <- as.data.frame(t(NSCLCcount[Gene,rownames(sv1)]))
# expr1 <- as.data.frame(t(NSCLCcount[Gene,]))
table(substr(rownames(expr1), 14,15))
expr1$Group <- ifelse(substr(rownames(expr1), 14,15) == 11 ,"Normal","Tumor")
table(expr1$Group)
Gene


# 箱图
#----------------------
p1 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=MS4A7,colour = Group ), #箱线图
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=MS4A7,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=MS4A7), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='MS4A7',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.y = element_text(size = 12, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold.italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_blank(),#删除X轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p1

p2 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=RETN,colour = Group ), #箱线图
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=RETN,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=RETN), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='RETN',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.y = element_blank(),#删除y轴标签
        axis.title.x = element_blank(),#删除X轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p2

p3 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=CXCR2,colour = Group ), #箱线图
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=CXCR2,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=CXCR2), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='CXCR2',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.y = element_blank(),#删除y轴标签
        axis.title.x = element_blank(),#删除X轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p3

p4 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=CD177,colour = Group ), #箱线图
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=CD177,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=CD177), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='CD177',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.y = element_text(size = 12, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold.italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_blank(),#删除X轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p4

p5 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=CSRNP1,colour = Group ), #箱线图
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=CSRNP1,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=CSRNP1), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='CSRNP1',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.x = element_blank(),#删除X轴标签
        axis.title.y = element_blank(),#删除y轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p5

p6 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=LUCAT1,colour = Group ), #箱线图
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=LUCAT1,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=LUCAT1), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='LUCAT1',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.y = element_blank(),#删除y轴标签
        axis.title.x = element_blank(),#删除X轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p6

plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)


# 分半小提琴图
library(devtools)
# devtools::install_github("psyteachr/introdataviz")
library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)

gene <- rownames(gene_multi)
gene1 <- as.data.frame(t(NSCLCcount[gene[1],]))
names(gene1) <- "Expression"
gene1$Gene <- rep(gene[1],length(rownames(gene1)))
gene1$Type <- ifelse(substr(rownames(gene1), 14,15) == 11 ,"Normal","Tumor")


gene2 <- as.data.frame(t(NSCLCcount[gene[2],]))
names(gene2) <- "Expression"
gene2$Gene <- rep(gene[2],length(rownames(gene2)))
gene2$Type <- ifelse(substr(rownames(gene2), 14,15) == 11 ,"Normal","Tumor")


gene3 <- as.data.frame(t(NSCLCcount[gene[3],]))
names(gene3) <- "Expression"
gene3$Gene <- rep(gene[3],length(rownames(gene3)))
gene3$Type <- ifelse(substr(rownames(gene3), 14,15) == 11 ,"Normal","Tumor")


gene4 <- as.data.frame(t(NSCLCcount[gene[4],]))
names(gene4) <- "Expression"
gene4$Gene <- rep(gene[4],length(rownames(gene4)))
gene4$Type <- ifelse(substr(rownames(gene4), 14,15) == 11 ,"Normal","Tumor")


gene5 <- as.data.frame(t(NSCLCcount[gene[5],]))
names(gene5) <- "Expression"
gene5$Gene <- rep(gene[5],length(rownames(gene5)))
gene5$Type <- ifelse(substr(rownames(gene5), 14,15) == 11 ,"Normal","Tumor")


gene6 <- as.data.frame(t(NSCLCcount[gene[6],]))
names(gene6) <- "Expression"
gene6$Gene <- rep(gene[6],length(rownames(gene6)))
gene6$Type <- ifelse(substr(rownames(gene6), 14,15) == 11 ,"Normal","Tumor")

df <- rbind(gene1,gene2,gene3,gene4,gene5,gene6)
head(df)
colnames(df)
df$ID <- rownames(df)
rownames(df) <- NULL

# pal = pal_lancet(palette = 'lanonc')(9)
# pal = pal_npg(palette = 'nrc')(9)
mypalette = pal_ucscgb()(10)
ggplot(df,aes(x = Gene, y = Expression, fill = Type)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
        legend.position = 'top') + # 图例的位置
  # scale_fill_brewer(palette = 'Set1') +
  # scale_fill_jco(name = '') +
  labs(y="Expression(log2(Count+1))",x=NULL)+ # 添加标题，x轴，y轴内容
  # scale_color_lancet(name = '')+
  scale_fill_manual(values = mypalette[c(6,1)]) + 
  # scale_fill_manual(values = c(pal[1],pal[2]),name = '') +
  # ylim(1,8) + #限定y轴范围
  # 添加显著性标记
  stat_compare_means(aes(group=Type),
                     # symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "NS")),
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "NS")),
                     label.y = 18,size = 5)

ggsave(filename = '6gene.pdf', width = 8,height = 6, path = '../../Fig/Step4-4/')



# KM生存曲线
load("../Step4-1/Train_Result.Rdata")
load(file = "../Step4-2/Verify_TCGA283.Rdata")
Gene <- rownames(gene_multi)
Gene
# [1] "MS4A7"  "RETN"   "CXCR2"  "CD177"  "CSRNP1" "LUCAT1"
dat1 = data.frame(MS4A7 = sv2$MS4A7, RETN = sv2$RETN, CXCR2 = sv2$CXCR2, CD177 = sv2$CD177, CSRNP1 = sv2$CSRNP1,
                  LUCAT1 = sv2$LUCAT1, OS = sv2$OS, OS.time = sv2$OS.time, row.names = rownames(sv2), id = rownames(sv2))
dat2 = data.frame(MS4A7 = data_4$MS4A7, RETN = data_4$RETN, CXCR2 = data_4$CXCR2, CD177 = data_4$CD177, CSRNP1 = data_4$CSRNP1,
                 LUCAT1 = data_4$LUCAT1, OS = data_4$OS, OS.time = data_4$OS.time, row.names = rownames(data_4), id = rownames(data_4))
dat3 = rbind(dat1,dat2)
dat4 = data.frame(OS = dat3$OS, OS.time = dat3$OS.time,row.names = rownames(dat3))
dat4$OS.time=dat4$OS.time/365


# MS4A7  
# 最佳生存节点
library(survival)
library(survminer)
library(ggsci)
res.cut <- surv_cutpoint(dat3, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "MS4A7") #需要计算的数据列名
summary(res.cut) #查看数据最佳截断点及统计量
MS4A7 = ifelse(dat3$MS4A7<res.cut$cutpoint[,1],"Low","High")
MS4A7 = factor(MS4A7)
#KM
sfit <- survfit(Surv(OS.time, OS)~MS4A7, data=dat4)
p1 = ggsurvplot(sfit,
           palette = 'jco', 
           conf.int = T,conf.int.style='step', 
           pval = T,pval.method = T,
           # risk.table = T,risk.table.pos='in',
           legend=c(0.85,0.85),
           legend.title="MS4A7",
           legend.labs=c("High","Low"),
           title="Survival Curve for MS4A7", 
           xlab ="Time(Years)",
           surv.median.line = "hv",
           ggtheme = theme_bw(base_size = 12))
p1
ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')


# RETN  
# 最佳生存节点
res.cut <- surv_cutpoint(dat3, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "RETN") #需要计算的数据列名
summary(res.cut) #查看数据最佳截断点及统计量
RETN = ifelse(dat3$RETN<res.cut$cutpoint[,1],"Low","High")
RETN = factor(RETN)
#KM
sfit <- survfit(Surv(OS.time, OS)~RETN, data=dat4)
p2 = ggsurvplot(sfit,
                palette = 'jco', 
                conf.int = T,conf.int.style='step', 
                pval = T,pval.method = T,
                # risk.table = T,risk.table.pos='in',
                legend=c(0.85,0.85),
                legend.title="RETN",
                legend.labs=c("High","Low"),
                title="Survival Curve for RETN", 
                xlab ="Time(Years)",
                surv.median.line = "hv",
                ggtheme = theme_bw(base_size = 12))
p2
ggsave(filename = 'RETN_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')


# CXCR2  
# 最佳生存节点
res.cut <- surv_cutpoint(dat3, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "CXCR2") #需要计算的数据列名
summary(res.cut) #查看数据最佳截断点及统计量
CXCR2 = ifelse(dat3$CXCR2<res.cut$cutpoint[,1],"Low","High")
CXCR2 = factor(CXCR2)
#KM
sfit <- survfit(Surv(OS.time, OS)~CXCR2, data=dat4)
p3 = ggsurvplot(sfit,
                palette = 'jco', 
                conf.int = T,conf.int.style='step', 
                pval = T,pval.method = T,
                # risk.table = T,risk.table.pos='in',
                legend=c(0.85,0.85),
                legend.title="CXCR2",
                legend.labs=c("High","Low"),
                title="Survival Curve for CXCR2", 
                xlab ="Time(Years)",
                surv.median.line = "hv",
                ggtheme = theme_bw(base_size = 12))
p3
ggsave(filename = 'CXCR2_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')


# CD177  
# 最佳生存节点
res.cut <- surv_cutpoint(dat3, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "CD177") #需要计算的数据列名
summary(res.cut) #查看数据最佳截断点及统计量
CD177 = ifelse(dat3$CD177<res.cut$cutpoint[,1],"Low","High")
CD177 = factor(CD177)
#KM
sfit <- survfit(Surv(OS.time, OS)~CD177, data=dat4)
p4 = ggsurvplot(sfit,
                palette = 'jco', 
                conf.int = T,conf.int.style='step', 
                pval = T,pval.method = T,
                # risk.table = T,risk.table.pos='in',
                legend=c(0.85,0.85),
                legend.title="CD177",
                legend.labs=c("High","Low"),
                title="Survival Curve for CD177", 
                xlab ="Time(Years)",
                surv.median.line = "hv",
                ggtheme = theme_bw(base_size = 12))
p4
ggsave(filename = 'CD177_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')


# CSRNP1  
# 最佳生存节点
res.cut <- surv_cutpoint(dat3, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "CSRNP1") #需要计算的数据列名
summary(res.cut) #查看数据最佳截断点及统计量
CSRNP1 = ifelse(dat3$CSRNP1<res.cut$cutpoint[,1],"Low","High")
CSRNP1 = factor(CSRNP1)
#KM
sfit <- survfit(Surv(OS.time, OS)~CSRNP1, data=dat4)
p5 = ggsurvplot(sfit,
                palette = 'jco', 
                conf.int = T,conf.int.style='step', 
                pval = T,pval.method = T,
                # risk.table = T,risk.table.pos='in',
                legend=c(0.85,0.85),
                legend.title="CSRNP1",
                legend.labs=c("High","Low"),
                title="Survival Curve for CSRNP1", 
                xlab ="Time(Years)",
                surv.median.line = "hv",
                ggtheme = theme_bw(base_size = 12))
p5
ggsave(filename = 'CSRNP1_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')


# LUCAT1  
# 最佳生存节点
res.cut <- surv_cutpoint(dat3, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "LUCAT1") #需要计算的数据列名
summary(res.cut) #查看数据最佳截断点及统计量
LUCAT1 = ifelse(dat3$LUCAT1<res.cut$cutpoint[,1],"Low","High")
LUCAT1 = factor(LUCAT1)
#KM
sfit <- survfit(Surv(OS.time, OS)~LUCAT1, data=dat4)
p6 = ggsurvplot(sfit,
                palette = 'jco', 
                conf.int = T,conf.int.style='step', 
                pval = T,pval.method = T,
                # risk.table = T,risk.table.pos='in',
                legend=c(0.85,0.85),
                legend.title="LUCAT1",
                legend.labs=c("High","Low"),
                title="Survival Curve for LUCAT1", 
                xlab ="Time(Years)",
                surv.median.line = "hv",
                ggtheme = theme_bw(base_size = 12))
p6
ggsave(filename = 'LUCAT1_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')

#拼接图片
splots <- list()
splots[[1]] <- p1
splots[[2]] <- p5
splots[[3]] <- p6
splots[[4]] <- p2
splots[[5]] <- p3
splots[[6]] <- p4
# 将多个图合并一起
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 3, nrow = 2) #定义行数和列数


# scRNA 中看表达量
load(file = '../Step2-2/res0.4_NSCLC.Integrate_28cell.Rdata')
genes_to_check <- sort(c('CD177', 'CSRNP1', 'CXCR2', 'LUCAT1', 'MS4A7', 'RETN'))
library(ggplot2)
library(reshape2)
library(ggsci)
library(dplyr)
vln.df=as.data.frame(datSet[["RNA"]]@data[genes_to_check,])
vln.df[1:6,1:6]
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df ,id="gene") # 转化为包含gene cell exp 三列的数据框
colnames(vln.df)[c(2,3)]= c("CB", "exp")
head(vln.df)
vln.df[1:10,]
anno=data.frame(Cell = datSet@active.ident)
anno$CB <- rownames(anno)
vln.df=inner_join(vln.df ,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = genes_to_check)#为了控制画图的基因顺序

pdf(file = '../../Fig/Step4-4/6gene_scRNA.pdf',width = 12,height = 10)
vln.df %>% ggplot(aes(Cell,exp)) + geom_violin(aes(fill=gene),scale = "width") + 
  facet_grid(vln.df$gene~. ,scales = "free_y")+
  # scale_fill_brewer(palette = "set3",direction = 1) +
  scale_x_discrete("") + scale_y_continuous("")+
  theme_bw()+
  labs(title='scRNA Expression')+ # 添加标题，x轴，y轴内容
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1, face = "bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none", plot.title=element_text(hjust=0.5,face = "bold.italic",size = 14))+
  scale_fill_lancet()
dev.off()




#堆叠小提琴图
library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#配色方案
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
pdf(file = '../../Fig/Step4-4//6gene_scRNA.pdf',width = 12,height = 10)
StackedVlnPlot(datSet, genes_to_check, pt.size=0, cols=my36colors)
dev.off()

