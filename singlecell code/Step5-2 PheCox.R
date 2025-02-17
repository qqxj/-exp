

### Step5-2  探索预后模型的临床应用
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step5-2/')
rm(list = ls())
options(stringsAsFactors = F)


load(file = "../Step2-0/TCGA_Step1 output.Rdata")
load(file = "../Step2-0/TCGA_Step1 PheData.Rdata")
load(file = "../Step4-1/Train_Result.Rdata") 
df <- read.csv(file = '../Step5-1/TCGA_scale.csv')
rownames(df) <- df$X
df <- df[,-1]
df[1:6,1:6]
exprSet <- as.data.frame(t(df[rownames(gene_multi),]))
dim(exprSet)
exprSet$RiskScore <- exprSet[,rownames(gene_multi)[1]]*gene_multi[1,1]+exprSet[,rownames(gene_multi)[2]]*gene_multi[2,1]+exprSet[,rownames(gene_multi)[3]]*gene_multi[3,1]+exprSet[,rownames(gene_multi)[4]]*gene_multi[4,1]+exprSet[,rownames(gene_multi)[5]]*gene_multi[5,1]+exprSet[,rownames(gene_multi)[6]]*gene_multi[6,1]
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore) , "Low","High")
table(exprSet$RiskGroup)
exprSet$id <- rownames(exprSet)
table(substr(rownames(exprSet),14,15))
exprSet$Type <- ifelse(substr(rownames(exprSet),14,15) == 11, 'Normal','Tumor')
table(exprSet$Type)
# exprSet <- subset(exprSet,exprSet$Type == 'Tumor')
table(exprSet$RiskGroup)
rownames(dat) <- dat$Id
phe <- dat[rownames(exprSet),]
df2 <- data.frame(RiskScore = exprSet$RiskScore, RiskGroup = exprSet$RiskGroup, row.names = rownames(exprSet),Id = exprSet$id)
phe <- merge(phe,df2,by = 'Id')
rownames(phe) <- phe$Id

table(phe$OS)
phe$OS <- ifelse(phe$OS == "Alive",0 ,1)
table(phe$OS)
phe$OS.time <- ifelse(is.na(phe$`Days to Death`),phe$`Days to Last Follow`,phe$`Days to Death`)
table(is.na(phe$OS.time))
phe <- subset(phe,phe$OS.time != 'NA')
library(ggstatsplot)
library(tidyverse) 
library(ggpubr)
datSet <- phe
# Gender
table(datSet$Gender)
datSet$Gender <- ifelse(datSet$Gender == "female" , "Female", "Male")
p1 = ggbetweenstats(data = datSet, x = Gender,y = RiskScore,
               bf.message = FALSE,#去除贝叶斯相关的统计值
               type = "nonparamatric",#选择非参数检验
               package = "ggsci",#调用调色板
               palette = "lanonc_lancet", # nrc_npg
               title = 'Gender')+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
p1
ggsave(filename = 'Gender.pdf',plot = p1, width = 8,height = 6,path = '../../Fig/Step5-2/')

# Age
table(is.na(datSet$Age))
datSet$Age <- ifelse(datSet$Age <= 60 , '<=60','>60')
table(datSet$Age)
p2 = ggbetweenstats(data = datSet, x = Age,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric",#选择非参数检验
                    package = "ggsci",#调用调色板
                    palette = "lanonc_lancet",
                    title = 'Age')+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
p2
ggsave(filename = 'Age.pdf',plot = p2, width = 8,height = 6,path = '../../Fig/Step5-2/')

# M
colnames(datSet)[3:5] <- c("M Stage", "N Stage", "T Stage")
table(datSet$`M Stage` == '')
datSet2 <- subset(datSet,datSet$`M Stage` != '')
table(datSet2$`M Stage`)
datSet2$`M Stage` <- ifelse(datSet2$`M Stage` == 'MX','M0' ,
                            ifelse(datSet2$`M Stage` == 'M0', 'M0','M1'))
datSet2 <- subset(datSet2,datSet2$`M Stage` != 'MX')
table(datSet2$`M Stage`)
p3 = ggbetweenstats(data = datSet2, x = `M Stage`,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric", #选择非参数检验
                    grouping.var = `M Stage`,# 分组变量
                    package = "ggsci",#调用调色板
                    palette = "lanonc_lancet",
                    title = 'M State')+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
                    # outlier.tagging = TRUE,#标记异常值
                    # p.adjust.method = "fdr",
                    # pairwise.comparisons = TRUE)
    # ggpubr::stat_compare_means(comparisons=my_comparisons, label.y = c(0.08, 0.12, 0.18),
    #                          method = "wilcox.test",
    #                          label = 'p.signif')
p3
ggsave(filename = 'M Stage.pdf',plot = p3, width = 8,height = 6,path = '../../Fig/Step5-2/')

# N Stage
table(datSet$`N Stage` == '')
datSet3 <- subset(datSet,datSet$`N Stage` != '')
# 
# my_comparisons2 <-list(c("N0", "N1"),
#                        c("N1", "N2"),
#                        c("N0", "N2"))
p4 = ggbetweenstats(data = datSet3, x = `N Stage`,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric", #选择非参数检验
                    grouping.var = `N Stage`, # 分组变量
                    package = "ggsci", #调用调色板
                    palette = "lanonc_lancet",
                    title = 'N Stage',
                    pairwise.comparisons=TRUE,
                    pairwise.display ="all")+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
# outlier.tagging = TRUE,#标记异常值
# p.adjust.method = "fdr",
# pairwise.comparisons = TRUE)
  # ggpubr::stat_compare_means(comparisons=my_comparisons2, label.y = c(0.08, 0.12, 0.18),
  #                            method = "wilcox.test",
  #                            label = 'p.signif')
p4
ggsave(filename = 'N Stage.pdf',plot = p4, width = 8,height = 6,path = '../../Fig/Step5-2/')

# T Stage
table(datSet$`T Stage`)
datSet4 <- subset(datSet,datSet$`T Stage` != 'TX')
table(datSet4$`T Stage`)

# my_comparisons3 <- list(c("T1", "T2"),
#                         c("T2", "T3"),
#                         c("T3", "T4"),
#                         c("T1", "T3"),
#                         c("T1", "T4"),
#                         c("T2", "T4"))
p5 = ggbetweenstats(data = datSet4, x = `T Stage`,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric", #选择非参数检验
                    grouping.var = `T Stage`, # 分组变量
                    package = "ggsci", #调用调色板
                    palette = "lanonc_lancet",
                    title = 'T Stage',
                    pairwise.comparisons=TRUE,
                    pairwise.display ="all")+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
  # outlier.tagging = TRUE,#标记异常值
  # p.adjust.method = "fdr",
  # pairwise.comparisons = TRUE)
  # ggpubr::stat_compare_means(comparisons=my_comparisons3, label.y = c(0.08, 0.12, 0.18, 0.22, 0.26,
  #                                                                     0.30, 0.34, 0.38, 0.42, 0.52),
  #                            method = "wilcox.test",
  #                            label = 'p.signif') # p.signif p.format
p5
ggsave(filename = 'T Stage.pdf',plot = p5, width = 8,height = 6,path = '../../Fig/Step5-2/')

# Stage
table(datSet$Stage)
table(datSet$Stage != 'not reported')
datSet5 <- subset(datSet,datSet$Stage != 'not reported')
table(datSet5$Stage)
datSet5$Stage <- ifelse(datSet5$Stage == ' stage i','I' ,
                        ifelse(datSet5$Stage == ' stage ii' ,'II',
                               ifelse(datSet5$Stage == ' stage iii', 'III', 'IV')))
table(datSet5$Stage)

p6 = ggbetweenstats(data = datSet5, x = Stage,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric", #选择非参数检验
                    grouping.var = Stage, # 分组变量
                    package = "ggsci", #调用调色板
                    palette = "lanonc_lancet",
                    title = 'Stage',
                    pairwise.comparisons=TRUE,
                    pairwise.display ="all")+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
# outlier.tagging = TRUE,#标记异常值
# p.adjust.method = "fdr",
# pairwise.comparisons = TRUE)
# ggpubr::stat_compare_means(comparisons=my_comparisons3, label.y = c(0.08, 0.12, 0.18, 0.22, 0.26,
#                                                                     0.30, 0.34, 0.38, 0.42, 0.52),
#                            method = "wilcox.test",
#                            label = 'p.signif') # p.signif p.format
p6
ggsave(filename = 'Tumor Stage.pdf',plot = p6, width = 8,height = 6,path = '../../Fig/Step5-2/')

# p6 = ggplot(datSet4,aes(Stage,RiskScore))+
#   geom_violin(aes(fill=Stage),cex=0.8)+  #Age，其实用R默认颜色也不错，这里只是展示一下如何提取喜欢的图片颜色。
#   scale_fill_manual(values = c('#FB5554','#42F203','#579ABB','#B978AE'))+
#   geom_boxplot(width=0.2,cex=0.8)+
#   theme_classic(base_size = 20)+
#   theme(axis.text = element_text(color = 'black'),
#         legend.position = 'none')+
#   labs(title="Kruskal-Wallis test   p = 0.012")+ # 添加标题，x轴，y轴内容
#   theme(plot.title = element_text(size = 18,
#                                   colour = "black",
#                                   face = "bold",
#                                   vjust = 1.9,
#                                   hjust = 0.5))+
#   geom_signif(mapping=aes(x=Stage,y=RiskScore), # 不同组别的显著性
#               comparisons = list(c("I", "II"),
#                                  c("I", "III"),
#                                  c("I", "IV"),
#                                  c("II", "III"),
#                                  c("II", "IV"),
#                                  c("III", "IV")),# 哪些组进行比较
#               map_signif_level=T, # T显示显著性，F显示p value
#               tip_length=c(0.05,0.05, 0.05,0.05, 0.05,0.05, 0.05,0.05,0.05,0.05,
#                            0.05,0.05), # 修改显著性线两端的长短
#               y_position = c(4.2,4.6,5,
#                              5.4,5.8,
#                              6.4), # 设置显著性线的位置高度
#               size=1, # 修改线的粗细
#               textsize = 5, # 修改显著性标记的大小
#               test = "wilcox.test")# 检验的类型
# 
# p6

library(cowplot)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)




# 临床信息 单因素Cox
datSet0 <- datSet
table(datSet0$Gender)
datSet0$Gender <- ifelse(datSet0$Gender == 'Female',1,2) #女性为1，男性为2
table(datSet0$Gender)

table(datSet0$Age)
datSet0$Age <- ifelse(datSet0$Age == '<=60',1,2) # <60为1， >60为2
table(datSet0$Age)

table(datSet0$`M Stage`)
table(datSet0$`M Stage` != '')
datSet0 <- subset(datSet0,datSet0$`M Stage` != '')
datSet0 <- subset(datSet0,datSet0$`M Stage` != 'MX')
datSet0$`M Stage` <- ifelse(datSet0$`M Stage` == "M0",1,
                            ifelse(datSet0$`M Stage` == "M1",2,3)) # M0-0  M1-1  MX-2
table(datSet0$`M Stage`)

table(datSet0$`N Stage`)
datSet0$`N Stage` <- ifelse(datSet0$`N Stage` == "N0",1,
                            ifelse(datSet0$`N Stage` == "N1",2,3))
table(datSet0$`N Stage`)

table(datSet0$`T Stage`)
datSet0 <- subset(datSet0,datSet0$`T Stage` != 'TX')
datSet0$`T Stage` <- ifelse(datSet0$`T Stage` == "T1",1,
                            ifelse(datSet0$`T Stage`== "T2",2,
                                   ifelse(datSet0$`T Stage` == "T3",3,
                                          ifelse(datSet0$`T Stage` == "T4",4,5))))
table(datSet0$`T Stage`)

table(datSet0$Stage)
datSet0 <- subset(datSet0,datSet0$Stage != 'not reported')
datSet0$Stage <- ifelse(datSet0$Stage == " stage i",1,
                        ifelse(datSet0$Stage == " stage ii",2,
                               ifelse(datSet0$Stage == " stage iii",3,4)))
table(datSet0$Stage)

table(datSet0$RiskGroup)
datSet0$RiskGroup <- ifelse(datSet0$RiskGroup == 'Low',1,2)
table(datSet0$RiskGroup)

colnames(datSet0)
colnames(datSet0)[3:5] <- c("MStage","NStage","TStage")
covariates <- c("Gender","Age","MStage","NStage","TStage","Stage","RiskGroup")

library(survival)
library(survminer)
#分别对每一个变量，构建生存分析的公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time,OS)~', x)))

univ_formulas
#对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = datSet0)})
univ_models
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)",           "wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
str(univ_results)
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)


rownames(res)
res$Names <- c("Gender","Age","M Stage","N Stage","T Stage","Stage","Risk Score")
colnames(res)
res<-res[, c("Names","HR (95% CI for HR)","wald.test","p.value","beta")]

#############################################################
#对HR (95% CI for HR)做处理，得到HR和low .95和high .95
#当然也可以改计算univ_results这一步的代码，不要将HR和CI贴起来
############################################################
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)

HR1 <- HR
HR1 <- as.data.frame(t(HR1))
res$HR <- as.numeric(HR1$V1)
res$upper <- as.numeric(HR1$V3)
res$lower <- as.numeric(HR1$V2)
res$p.value <- as.numeric(res$p.value)
res$beta <- as.numeric(res$beta)
rownames(res) <- res$Names
str(res)

res$'P value' = ifelse(
  res$p.value < 0.001,
  "<0.001 ***",
  ifelse(
    res$p.value < 0.01,
    "<0.01  **",
    ifelse(
      res$p.value < 0.05,
      paste(round(res$p.value, 2), " *"),
      round(res$p.value, 2)
    )
  )
)


# # 作图1
# #左边和右边边距稍微留多一点来写变量名称，pvalue和HR
# par(mar=c(5,6,4,13))
# #先用小方块画出HR
# plot(as.numeric(HR[1,]),1:dim(HR)[2],
#      pch=15,cex=2,col="blue",bty='n',yaxt='n',ylab=NA,xlab="Hazard Ratio",
#      xlim=range(as.numeric(unlist(HR)))
# )
# #添加中线
# abline(v=1,col="blue",lwd=2,lty=2)
# 
# for(i in 1:ncol(HR)){
#   x=as.numeric(HR[2:3,i])
#   #循环画出CI
#   lines(x,c(i,i),col="blue")
#   #添加变量名
#   text(0.2,i,rownames(res)[i],xpd=T,adj = c(0,0))
#   #添加HR和CI
#   text(1.9,i,as.character(res[i,2]),xpd=T,adj = c(0,0))
#   #添加p值
#   text(2.8,i,as.numeric(res[i,9]),xpd=T,adj = c(0,0))
# }
# #添加标题
# text(1.9,ncol(HR)+1,"HR(95% CI)",xpd=T,adj = c(0,0),font=2)
# text(2.8,ncol(HR)+1,"P value",xpd=T,adj = c(0,0),font=2)
# text(0.2,ncol(HR)+1,"Names",xpd=T,adj = c(0,0),font=2)
# text(0.8,ncol(HR)+2,"Univariate Cox Regression Analysis",xpd=T,adj = c(0,0),font = 4)
# 
# lines(x = c(0.2,2.8),y=c(7,7),col="blue")


# 作图2
tabletext <- cbind(c("Variable",res$Names),
                   c("HR(95% CI)",res$`HR (95% CI for HR)`),
                   c("P Value",res$`P value`))
str(tabletext)
class(tabletext)
tabletext0 <- data.frame(mean = res$HR,lower =res$lower, upper = res$upper,row.names = rownames(res))


# jpeg(file = "results_Value_2.jpg",width =2000,height = 1800,units = "px",res =300) #结果保存
library(forestplot)
forestplot(tabletext,  #显示的文本
           mean=c(NA,tabletext0$mean),
           lower=c(NA,tabletext0$lower), upper=c(NA,tabletext0$upper),
           # c(NA,tabletext0$mean), #误差条的均值(此处为差值的中值)
           # c(NA,tabletext0$lower), #误差条的下界(此处为差值的25%分位数)
           # c(NA,tabletext0$upper), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           graph.pos=2,
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.3, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 3,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.8), cex = 1), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Hazard Ratio(HR)",#x轴的标题
           #xticks = F,
           is.summary = c(T, rep(F, 7)),
           align = "l",
           hrzl_lines = list(
             "1" = gpar(lty=1),
             "2" = gpar(lty=1),
             "9"= gpar(lty=1)),
           colgap = unit(6, 'mm'),
           title="Univariate Cox Regression Analysis",
)
# dev.off()




# 多因素Cox
  library(glmnet)
  library(survival)
  library(survminer)

  formula_for_multivarirate <- as.formula(paste0('Surv(OS.time, OS)~',paste(covariates,sep = '',collapse = "+")))
  multi_varirate_cox <- coxph(formula_for_multivarirate, data = datSet0)
  ph_hypo_multi <- cox.zph(multi_varirate_cox)
  ph_hypo_table <-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),] ######
  #formula_for_multivarirate <- as.formula(paste0('Surv(os.time, os)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = " ",collapse = "+")))
  #multi_varirate_cox <- coxph(formula_for_multivarirate,data = sv1)
  multiCoxSum <- summary(multi_varirate_cox)
  multi_cox<- as.data.frame(multi_varirate_cox$coefficients)
  multi_cox
  
  #correlation <- cor(sv1[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],method = 'pearson')
  # correlation <- cor(datSet0[,rownames(ph_hypo_table)],method = 'pearson')
  # correlation <- cor(sv1[,rownames(ph_hypo_table)],method = 'pearson')
  # install.packages('GGally')
  # library('GGally')
  # ggpairs(datSet0[,rownames(ph_hypo_table)],
  #         axisLabels = "show")+
  #   theme_bw()+
  #   theme(panel.background = element_rect(colour = 'red',size = 1,fill = "white"),
  #         panel.grid = element_blank())
  
  # # vifϵ??
  # library("rms")
  # vif <- rms::vif(multi_varirate_cox)
  # sqrt(vif) < 2
  # library(survival)
  # library(survminer)
  # ggforest(model = multi_varirate_cox,data = datSet0, main =  "Hazard",fontsize = 1)
  # C_index <- multi_varirate_cox$concordance["concordance"]
  # if(C_index >= 0.9){print("high accuracy")
  # }else{
  #   if(C_index <0.9 & C_index >= 0.7){
  #     print("Medium accuracy")
  #   }else{print('low accuracy')
  #   }
  # }
  
  out_multi <- cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  
  out_multi
  class(out_multi)
  out_multi <- as.data.frame(out_multi)
  
  #HR和它的置信区间
  dat2 = as.data.frame(round(multiCoxSum$conf.int[, c(1, 3, 4)], 2))
  dat2 = tibble::rownames_to_column(dat2, var = "Variable")
  colnames(dat2)[2:4] = c("HR", "lower", "upper")
  #需要在图上显示的HR文字和p值
  dat2$'HR(95% CI)' = paste0(dat2[, 2], "(", dat2[, 3], "-", dat2[, 4], ")")
  dat2$'P Value' = out_multi$pvalue
  
  dat2$`P Value` = ifelse(
    dat2$`P Value` < 0.001,
    "<0.001 ***",
    ifelse(
      dat2$`P Value` < 0.01,
      "<0.01  **",
      ifelse(
        dat2$`P Value` < 0.05,
        paste(round(dat2$`P Value`, 2), " *"),
        round(dat2$`P Value`, 2)
      )
    )
  )
  dat2$Variable <- c("Gender","Age", "M Stage","N Stage","T Stage","Stage","Risk Score")
  str(dat2)
  
  dat3 <- as.data.frame(t(dat2))
  dat3$V8 <- rownames(dat3)
  dat3 <- dat3[,c(8,1:7)]
  dat3 <- as.data.frame(t(dat3))
  rownames(dat3) <- NULL
  colnames(dat3) <- NULL
  dat4 <- dat3[,c(1,5,6)]
  
  
  library(forestplot)
  forestplot(dat4,  #显示的文本
             mean=c(NA,dat2$HR),
             lower=c(NA,dat2$lower), upper=c(NA,dat2$upper),
             # c(NA,tabletext0$mean), #误差条的均值(此处为差值的中值)
             # c(NA,tabletext0$lower), #误差条的下界(此处为差值的25%分位数)
             # c(NA,tabletext0$upper), #误差条的上界(此处为差值的75%分位数)
             zero = 1, #显示y=0的垂直线
             xlog=FALSE, #x轴的坐标不取对数
             graph.pos=2,
             fn.ci_norm = fpDrawCircleCI, #误差条显示方式
             boxsize = 0.3, ##误差条中的圆心点大小
             col=fpColors(line = "#CC79A7", #误差条的线的颜色
                          box="#D55E00"), #误差条的圆心点的颜色
             lty.ci = 7,   # 误差条的线的线型
             lwd.ci = 3,   # 误差条的线的宽度
             ci.vertices.height = 0.15, # # 误差条末端的长度
             txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.8), cex = 1), #文本大小设置
             lineheight = "auto", #线的高度
             xlab="Hazard Ratio(HR)",#x轴的标题
             #xticks = F,
             is.summary = c(T, rep(F, 7)),
             align = "l",
             hrzl_lines = list(
               "1" = gpar(lty=1),
               "2" = gpar(lty=1),
               "9"= gpar(lty=1)),
             colgap = unit(6, 'mm'),
             title="Multivariate Cox Regression Analysis"
  )

