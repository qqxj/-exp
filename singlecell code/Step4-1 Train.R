

### Step4-1  构建预后模型
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step4-1/')

rm(list = ls())
options(stringsAsFactors = F)



#Step1  差异分析
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
library(dplyr)
library(DESeq2)

dat <- NSCLCcount
# dat <- na.omit(dat)
dat <- 2^dat -1 # FPKM转位Count
dat <- round(dat,digits = 0) # 取整

dat[dat<0] <- 0
# 排查异常值：37819行 22378253076
dat <- dat[-37819,]
All_dat <- dat

Group <- as.data.frame(group_list[,1])
Group$`group_list[, 1]` <- ifelse(Group[,1] == 1,"Tumor","Normol")
names(Group) <- "Sample"



#构建dds矩阵
countData <- All_dat
condition <- factor(Group$Sample)
head(condition)
countData[1:4,1:4]


# 排查异常值：37819行 22378253076
#  countData2 <- countData[37819:37819,]
#   countData2 <- as.data.frame(t(countData2))
# # countData2[countData2 > 1249869] <- 0

# dds <- DESeqDataSetFromMatrix(countData, colData = Group, design= ~ condition )
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
head(dds)
dim(dds)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)

DEG_INFO <- as.data.frame(res)

# 提取差异分析结果
# 获取padj（p值经过多重校验校正后的值）小于0.01，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
table(res$padj<0.05) #取P值小于0.01的结果
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)# 所有差异基因

Up_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange > 2)) 
Up_gene_deseq2 <- row.names(Up_gene_deseq2) # 上调差异基因

Down_gene_deseq2 <-  subset(res,padj < 0.05 & (log2FoldChange <  -2))
Down_gene_deseq2 <- row.names(Down_gene_deseq2) # 下调差异基因


resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)


save(resdata,diff_gene_deseq2,Up_gene_deseq2,Down_gene_deseq2,file = "DEseq2.Rdata")

# # 得到csv格式的差异表达分析结果
# write.csv(resdata,file= "output/control_vs_akap95.cvs",row.names = F)


# 火山图
# library(ggplot2)
# for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
#                           'padj' = res$padj,
#                           'trend' = rep('No', length(res$log2FoldChange)))
# up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 2), which(for_volcano$padj < 0.05))
# down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -2), which(for_volcano$padj < 0.05))
# for_volcano[up_sig_indices,'trend'] <- 'Up'
# for_volcano[down_sig_indices,'trend'] <- 'Down'
# for_volcano$trend <- as.factor(for_volcano$trend)
# for_volcano$padj <- -log10(for_volcano$padj)
# 
# p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = trend))+
#   geom_point(size = I(0.7))+
#   scale_color_manual(values = c('No'='black', 'Up' = 'red', 'Down' = 'blue'))+
#   geom_vline(xintercept = c(2, -2), lty=2, size=I(0.4), colour = 'grey11')+
#   geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2, size=I(0.1),colour = 'grey11')+
#   theme_bw()+
#   theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'),
#         panel.grid = element_blank())+
#   labs(x='log2FoldChange', y = '-log10Pvalue')
# 
# p
# ggsave(filename = "../../Fig/Step4-1/DEgene.png",plot = p,width = 8,height = 6)
# 
# 
# 
# # 2022.6.9
# load(file = "DEseq2.Rdata")
res = resdata
rownames(res) = res$Row.names
res = res[,-1]
library(ggplot2)
for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                          'padj' = res$padj,
                          'State' = rep('No', length(res$log2FoldChange)),
                          row.names = rownames(res))
up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 2), which(for_volcano$padj < 0.05))
down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -2), which(for_volcano$padj < 0.05))
for_volcano[up_sig_indices,'State'] <- 'Up'
for_volcano[down_sig_indices,'State'] <- 'Down'
for_volcano$State <- as.factor(for_volcano$State)
for_volcano$padj <- -log10(for_volcano$padj)

this_tile <- paste0('Cutoff for logFC is 2',
                    '\nThe number of Up gene is ',nrow(for_volcano[for_volcano$State =='Up',]) ,
                    '\nThe number of Down gene is ',nrow(for_volcano[for_volcano$State =='Down',]))


p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = State))+
  geom_point(size = I(1))+
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

ggsave(filename = 'DEG_vol.png',plot = p,width = 8,height = 6,path = "../../Fig/Step4-1/")





#Step2 交集基因
load(file = "../Step3-1/Neutrophil_Statemarker.Rdata")
load(file = "DEseq2.Rdata")

Statemarker0 <- unique(Statemarker$gene)
Intergene <- intersect(Statemarker0,diff_gene_deseq2)

library(VennDiagram)
venn_list <- list(NDRGs = Statemarker0,
                  Up_DEGs = Up_gene_deseq2,
                  Down_DEGs = Down_gene_deseq2)
venn.diagram(venn_list, filename = '../../Fig/Step4-1/IntersectGene.png', imagetype = 'png', 
             fill = c('red', 'green','purple'), alpha = 0.50, 
             cat.col = c('red', 'green','purple'), cat.cex = 1.2, cat.fontfamily = 'serif',cat.default.pos = "outer",cat.pos = c(-5, 0, 0),  # 位置，用圆的度数
             cat.dist = c(0.055, 0.055, 0.035),  # 位置，离圆的距离
             col = c('red', 'green','purple'), cex = 1.5, fontfamily = 'serif')
inter <- get.venn.partitions(venn_list) # 提取交集基因



# Step3 构建模型

# 测试集 
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
load(file = "../Step2-0/TCGA_Step1 PheData.Rdata")

##  7:3 拆分数据
colnames(dat)[7] <- "OS.time"
data = dat

# 设置OS.time为datatime+last fllowe time
data2 = as.data.frame(data$OS.time)
data2$Follwe <- data$`Days to Last Follow`
rownames(data2) <- data$Id
data2$time <- ifelse(is.na(data2$`data$OS.time`),
                     data2$Follwe,
                     data2$`data$OS.time`)
table(! is.na(data2$time))
data3 <- as.data.frame(row.names = rownames(data2),data2$time)
data3 <- na.omit(data3)
names(data3) <- "OS.time"
rownames(data) <- data$Id
data <- data[rownames(data3),]
data[,7] <- data3$OS.time


# 2/3 为测试集
train_sub_d = sample(nrow(data),2/3*nrow(data))
data9 <- data[train_sub_d,]
save(data9,train_sub_d,file = "Ver_2.3_d.Rdata")



# 开始训练
load(file = "Ver_2.3_d.Rdata") # 目前最好

#  临床数据整理
phe <- data9
# rownames(phe) <- phe$Id
# phe$OS.time <- ifelse(is.na(phe$`Days to Death`),phe$`Days to Last Follow`,phe$`Days to Death`)
# table(is.na(phe$OS.time))
# phe <- subset(phe, phe$OS.time != "NA")

table(phe$OS)
phe$OS <- ifelse(phe$OS == "Alive",0,1)

table(phe$OS.time == 0)
phe <- subset(phe, phe$OS.time != 0)

table(phe$OS.time <= 4000)
phe <- phe[phe$OS.time <= 4000,]


# 表达数据整理
expr <- as.data.frame(t(NSCLCcount[Intergene,rownames(phe)]))
colnames(expr)
refGene <- c("ACTB","GAPDH","TFRC","TUBB")
refExpr <- as.data.frame(t(NSCLCcount[refGene,rownames(phe)]))
Basal_CT <- cbind(refExpr,expr)
for (i in Intergene) {
  Basal_CT[,i] = ((Basal_CT[,i] - Basal_CT$ACTB)+(Basal_CT[,i]-Basal_CT$GAPDH)+(Basal_CT[,i]-Basal_CT$TFRC)+(Basal_CT[,i]-Basal_CT$TUBB))/4
  
}
expr0 <- Basal_CT[,Intergene]
expr0 <- scale(expr0)
expr0 <- as.data.frame(expr0)
expr0$Id = rownames(expr0)
# sv1 <- cbind(phe, expr0)
sv1 = inner_join(phe, expr0)
rownames(sv1) = sv1$Id
colnames(sv1)[1:11]
covariates <- colnames(sv1)[12:ncol(sv1)]
# TrainSet <- sv1
# save(TrainSet,file = "TrainSet_5gene.Rdata")

df <- sv1

# 弹性网络
{
library(glmnet)
library(survival)
library(survminer)

covariates2 <- covariates
# 弹性网络 测试alphe值

sv2 <- df[,covariates2]
sv2$OS.time <- df$OS.time
sv2$OS <- df$OS

x <- as.matrix(sv2[,covariates2])
table(is.na(x))
rownames(x) <- NULL
# x[x<0]<- 0
y <- as.matrix(survival::Surv(sv2$OS.time, sv2$OS))

}

# 使用caret包自动选择最佳的调整参数alpha和lambda
# 最佳的alpha和lambda值是那些使交叉验证误差最小的值
library(caret)
set.seed(68) 
{
  model <- train(
    OS.time ~.+ OS , data = sv2, method = "glmnet",
    trControl = trainControl("cv", number = 10),
    tuneLength = 10 # 测试10个不同的alpha和lambda值的组合
  )
  ap = model$bestTune
  ap[1,1]
  
  
  x <- as.matrix(sv2[ ,gsub(covariates,pattern = '-', replacement = '_')])
  fit <- glmnet(x,Surv(sv2$OS.time, sv2$OS),family = 'cox',alpha =ap[1,1])
  plot(fit,label = T,lwd=2,xlab="Log lambda") + abline(v=ap[1,1])
  box(lwd=4)
  
  y <- as.matrix(survival::Surv(sv2$OS.time, sv2$OS))
  lasso_fit <- cv.glmnet(x, y, family='cox', type.measure = 'deviance')
  plot(lasso_fit)
  box(lwd=4)
  
  coefficient <- coef(lasso_fit, s=lasso_fit$lambda.min)
  Active.Index <- which(as.numeric(coefficient) != 0)
  active.coefficients <- as.numeric(coefficient)[Active.Index]
  sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
  
  
  sig_gene_multi_cox

active.coefficients
length(sig_gene_multi_cox)
gene <- data.frame(sig_gene_multi_cox,active.coefficients)
rownames(gene) <- gene$sig_gene_multi_cox
gene
rownames(gene)
}



# 多因素Cox
{
  library(glmnet)
  library(survival)
  library(survminer)
  covariates3 <- rownames(gene)
  formula_for_multivarirate <- as.formula(paste0('Surv(OS.time, OS)~',paste(covariates3,sep = '',collapse = "+")))
  multi_varirate_cox <- coxph(formula_for_multivarirate, data = sv2)
  ph_hypo_multi <- cox.zph(multi_varirate_cox)
  ph_hypo_table <-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),] ######

  multiCoxSum <- summary(multi_varirate_cox)
  multi_cox<- as.data.frame(multi_varirate_cox$coefficients)
  multi_cox
  
  correlation <- cor(sv2[,rownames(ph_hypo_table)],method = 'pearson')

  library('GGally')
  ggpairs(sv2[,rownames(ph_hypo_table)],
          axisLabels = "show")+
    theme_bw()+
    theme(panel.background = element_rect(colour = 'red',size = 1,fill = "white"),
          panel.grid = element_blank())
  

  library("rms")
  vif <- rms::vif(multi_varirate_cox)
  sqrt(vif) < 2
  library(survival)
  library(survminer)
  ggforest(model = multi_varirate_cox,data = sv2, main =  "Hazard",fontsize = 1)
  
  
  C_index <- multi_varirate_cox$concordance["concordance"]
  if(C_index >= 0.9){print("high accuracy")
  }else{
    if(C_index <0.9 & C_index >= 0.7){
      print("Medium accuracy")
    }else{print('low accuracy')
    }
  }
  
  out_multi <- cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  
  out_multi
  class(out_multi)
  out_multi <- as.data.frame(out_multi)
  gene_multi <- subset(out_multi,out_multi$pvalue < 0.05)
  
  gene_multi
  rownames(gene_multi)
  
}

sv2$exp <- sv2[,rownames(gene_multi)[1]]*gene_multi[1,1]+sv2[,rownames(gene_multi)[2]]*gene_multi[2,1]+sv2[,rownames(gene_multi)[3]]*gene_multi[3,1]+sv2[,rownames(gene_multi)[4]]*gene_multi[4,1]+sv2[,rownames(gene_multi)[5]]*gene_multi[5,1]+sv2[,rownames(gene_multi)[6]]*gene_multi[6,1]
save(sv1,sv2,gene_multi,file = "Train_Result.Rdata")


load(file = "Train_Result.Rdata")
dat1<- sv2
# ROC 1,3,5
{library(survivalROC)
  cutoff <- 1*365
  ROC1= survivalROC(Stime=dat1$OS.time,##生存时间
                    status=dat1$OS,## 终止事件    
                    marker =dat1$exp, ## marker value    
                    predict.time = cutoff## 预测时间截点
                    ,method="KM")##span,NNE法的namda
  str(ROC1)## list结构
  
  cutoff <- 3*365
  ROC3= survivalROC(Stime=dat1$OS.time,##生存时间
                    status=dat1$OS,## 终止事件    
                    marker =dat1$exp, ## marker value    
                    predict.time = cutoff## 预测时间截点
                    ,method="KM")##span,NNE法的namda
  str(ROC3)## list结构
  
  cutoff <- 5*365
  ROC5= survivalROC(Stime=dat1$OS.time,##生存时间
                    status=dat1$OS,## 终止事件    
                    marker = dat1$exp, ## marker value    
                    predict.time = cutoff## 预测时间截点
                    ,method="KM")##span,NNE法的namda
  str(ROC5)## list结构
  
  
  plot(ROC5$FP, ROC5$TP, ## x=FP,y=TP
       type="l",col="blue",lwd=2,##线条设置
       xlim=c(0,1), ylim=c(0,1),   
       xlab=paste( "False positive rate"), ##连接
       ylab="True positive rate",
       main="ROC of Training Set"
       )## \n换行符
  # lines(ROC12$FP,ROC12$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC10$FP,ROC10$TP, type="l",col="yellowgreen",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC9$FP,ROC9$TP, type="l",col="pink",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC7$FP,ROC7$TP, type="l",col="red",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC5$FP,ROC5$TP, type="l",col="salmon",xlim=c(0,1), ylim=c(0,1),lwd=2)
  lines(ROC3$FP,ROC3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1),lwd=2)
  lines(ROC1$FP,ROC1$TP, type="l",col="red",xlim=c(0,1), ylim=c(0,1),lwd=2)
  
  legend(0.5,0.4,c(paste("AUC of 1 years =",round(ROC1$AUC,3)),
                   paste("AUC of 3 years =",round(ROC3$AUC,3)),
                   paste("AUC of 5 years =",round(ROC5$AUC,3))),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,col=c("red","green","blue"),
         bty = "n",# bty框的类型
         seg.len=1,cex=0.8)# 
  abline(0,1,col="black",lty=1,lwd=2)##线条颜色
  box(lwd=2)
}


##### 生存分析km
{
  # factGene <- c("CYGB","CLVS1","CBY3","CPNE6")
  
  fac_mix <- dat1
  library(survival)
  library(survminer)
  library(ggplotify)
  library(cowplot)
  library(Hmisc)
  library(pheatmap)
  library(gridExtra)
  #s = as.formula(paste('Surv(OS.time, OS)~', noquote(paste(factGene,collapse = ' + '))))
  #model <- coxph(s, data = fac_mix )
  #summary(model,data=fac_mix)
  # RiskScore <- predict(model,type = "risk")
  RiskScore <- fac_mix$exp
  #names(RiskScore) = rownames(fac_mix)
  fp <- RiskScore
  phe <- fac_mix
  phe$id <- rownames(phe)
  
  # 生存图 
  {
    #最佳节点
    # res.cut <- surv_cutpoint(dat1, #数据集
    #                          time = "OS.time", #生存时间
    #                          event = "OS", #生存状态
    #                          variables = "exp") #需要计算的数据列名
    # 
    # summary(res.cut) #查看数据最佳截断点及统计量
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    dat1$OS.time=dat1$OS.time/365
    sfit <- survfit(Surv(OS.time, OS)~RiskGroup, data=dat1)
    # ggsurvplot(sfit, pval=TRUE,xlab ="Time(Years)",surv.median.line = "hv",pval.method = T)
    # ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
    #            risk.table =TRUE,pval =TRUE,
    #            conf.int =TRUE,xlab ="Time in years",
    #            ggtheme =theme_light(),
    #            ncensor.plot = TRUE)
    ggsurvplot(sfit,
               palette = 'jco',
               conf.int = T,conf.int.style='step', 
               pval = T,pval.method = T,
               # risk.table = T,risk.table.pos='in',
               legend=c(0.85,0.85),
               legend.title="Risk Group",
               legend.labs=c("High","Low"),
               title="Survival Curve for Training Set", 
               xlab ="Time(Years)",
               surv.median.line = "hv",
               ggtheme = theme_bw(base_size = 12))
  }
  
  # 散点+热图 中位数
  {
    
    fp_dat=data.frame(patientid=phe$id,fp=phe$exp)
    # fp_dat$RiskGroup= ifelse(fp_dat$fp>= res.cut$cutpoint[,1],'High','Low')
    fp_dat$RiskGroup= ifelse(fp_dat$fp>= median(fp_dat$fp),'High','Low')
    library(dplyr)
    fp_dat=arrange(fp_dat,fp)
    
    sur_dat=data.frame(patientid=phe$id,time=phe$OS.time/365,Status=phe$OS)
    sur_dat$Status=ifelse(sur_dat$Status==0,'Alive','Dead')
    sur_dat$Status=factor(sur_dat$Status,levels = c("Dead","Alive"))
    sur_dat$time <- sur_dat$time
    rownames(sur_dat)=sur_dat$patientid
    sur_dat=sur_dat[fp_dat$patientid,]
    
    #exp_dat <- fac_mix[names(sort(fp)),c(factGene)]
    # exp_dat=sv1[,c("CYGB","CLVS1","CBY3","CPNE6")]
    exp_dat=dat1[,rownames(gene_multi)]
    rownames(exp_dat)=phe$id
    exp_dat=exp_dat[fp_dat$patientid,]
    fp_dat$patientid=1:length(fp)
    sur_dat$patientid=1:length(fp)
    rownames(exp_dat)=1:length(fp)
    rownames(sur_dat)=1:length(fp)
    
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp_dat$fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    
    
    # p1 <- ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=RiskGroup))+
    #   scale_colour_manual(values = c("#FF6666","#00CCCC"))+
    #   theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
    #   geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
    #   # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
    #   geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
    #              colour="black", linetype="dotted",size=0.8) +
    #   theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
    #         axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
    #         legend.text=element_text(size=14),legend.title =element_text(size=14))
    # p1
    
    library(ggsci)
    library(scales)
    
    fp_dat$catpo <- rep(median(fp_dat$fp),length(fp_dat$patientid))
    p1 <- ggplot(fp_dat,aes(patientid))+
      geom_ribbon(aes(ymin = catpo, ymax = fp, fill = RiskGroup), alpha = 1)+
      scale_fill_jco()+
      # scale_colour_manual(values = c("#FF6666","#00CCCC"))+
      theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
      geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
      # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
                 colour="black", linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p1
    
    
    pal <- pal_jco('default')(10)
    p2 <- ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=Status),size = 3)+theme_bw()+
      scale_colour_manual(values = c(pal[1],pal[2]))+
      labs(x="Patient ID(Increasing Risk Score)",y="Survival Time(Years)")+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),colour="black", 
                 linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p2
      
    # mycolors <- colorRampPalette(c("#64b5f6", "#fffde7", "#ff5252"), bias = 1.2)(100)
    # mycolors <- colorRampPalette(c("#FF6666", "White", "#00CCCC"), bias = 1.2)(100)
    mycolors <- colorRampPalette(c(pal[1], "White", pal[2]), bias = 1.2)(100)
    tmp=t(scale(exp_dat))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
                show_rownames = T,cluster_rows = F,
                fontsize_row = 14,fontsize = 14)
    p3
    plots = list(p1,p2,as.ggplot(as.grob(p3)))
    lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7)))
    grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
    # plots = list(p1,as.ggplot(as.grob(p3)))
    # lay1 = rbind(c(rep(1,7)),c(rep(2,7)))
    # grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10))
  }
}


# nomaogram
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2) 
library(survival)
library(rms)

train_phe = sv1[,c(2,6,7,9)]
train_phe$id = rownames(train_phe)
risk = data.frame(row.names = rownames(sv2), RiskScore = sv2$exp, id = rownames(sv1))
train_phe = merge(train_phe, risk , by = 'id')
rownames(train_phe) = train_phe$id
train_phe = train_phe[,-1]

{

  Nmdat <- train_phe
  # Nmdat <- Nmdat[,c('RiskScore','OS','OS.time')]
  # Nmdat$OS.time <- as.numeric(Nmdat$OS.time)*30
  dd <- datadist(Nmdat)
  options(datadist="dd")
  multivarl <- as.formula(paste0('Surv(OS.time,OS)~', 
                                 paste('RiskScore', sep = '', collapse = '+')))
  
  coxm_1 <- cph(formula = multivarl,data=Nmdat,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  nomo <- nomogram(coxm_1,fun = list(surv1,surv3,surv5),lp = T,
                   funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
                   maxscale = 100,fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  pdf(file = '../../Fig/Step4-1/nomaogram.pdf',width = 6,height = 5)
  plot(nomo,lplabel = 'Linear Preadictor',
       xfrac = .35,varname.label = T,varname.label.sep = '=',ia.space = .2,
       tck = NA,tcl = 0.2,lmgp = 0.3,
       points.label = 'Points',total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
  dev.off()
  
  ## 校准曲线
  ## 参数说明：
  ## 1、绘制校正曲线前需要在模型函数中添加参数x=T, y=T，详细内容参考帮助
  ## 2、u需要与之前模型中定义好的time.inc一致，即365或730；
  ## 3、m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）
  ## 而m代表每组的样本量数，因此m*3应该等于或近似等于样本量；
  ## 4、b代表最大再抽样的样本量
  f1 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=243, B=1000)
  f3 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=243, B=1000)
  f5 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=243, B=1000)
  
  pdf(file = '../Figer/nomaogram_calibrate.pdf',width = 6,height = 5)
  plot(cal1,lwd=1,lty=1, cex.axis = 1,cex.lab=1,
       errbar.col = '#666699',
       xlab='Nomogram-Predicted Probability',
       ylab='Actual',
       col = '#666699',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal3,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#339933',
       col = '#339933',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal5,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#FF0033',
       col = '#FF0033',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  abline(0,1,lty=1,lwd=1)
  legend("bottomright",legend=c("1 - year","3 - year","5 - year"), 
         col=c("#666699","#339933","#FF0033"),
         lty= 1 ,lwd= 4,
         bty = "n",
         seg.len=1,cex=1)
  dev.off()
}

