

### Step4-3  预后模型的验证
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step4-3/')

rm(list = ls())
options(stringsAsFactors = F)


library(GEOquery)
eSet <- getGEO('GSE81089', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b = eSet[[1]]
phe=pData(b)
GSE81089_phe <- phe[,c(1:2,48:50,54:57)]
colnames(GSE81089_phe)
names(GSE81089_phe) <- c("title","geo_accession","Age","OS","Gender","Stage","Surgery data","TorN","Vital data")
write.csv(GSE81089_phe,file = "GSE81089_phe.csv")# 导出数据集，计算OS.time
a <- read.csv(file = "GSE81089_phe.csv")
colnames(a)[11] <- "OS.time"
GSE81089_phe$OS.time <- a$OS.time
dat2 <- read.table(gzfile("GSE81089_readcounts_featurecounts.tsv.gz"),header = T)
rownames(dat2) <- dat2[,1]

# ID转换
library(biomaRt)
library(stringr)
listMarts()                                                     
plant <- useMart("ensembl")                                     
listDatasets(plant)                                             
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) 
listFilters(mart)                                               
inputid <-  rownames(dat2)
inputid[1:12]
hg_symbols <- getBM(attributes=c('hgnc_symbol',
                                 'transcript_biotype','ensembl_gene_id'),
                    filters= 'ensembl_gene_id', values = inputid, mart = mart) 
table(hg_symbols$transcript_biotype)
library(dplyr)
colnames(hg_symbols)
colnames(dat2)[1]
colnames(hg_symbols)[3] <- "Ensembl_gene_id"

ids <- data.frame(hg_symbols$Ensembl_gene_id,hg_symbols$hgnc_symbol)
names(ids) <- c("probe_id","symbol")
dat <- dat2
ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息

dat <- dat[,-1]
range(dat)
dat <- log2(dat+1)
range(dat)
GSE81089_expr <- dat

a <- GSE81089_phe[,1:2]
b <- as.data.frame(t(GSE81089_expr))
b$title <- rownames(b)
d <- intersect(rownames(b),a$title)
library(dplyr)
c <- inner_join(a,b,by = "title")
c <-as.data.frame(t(c))
colnames(c) <- c[2,]
c <- c[-c(1:2),]

GSE81089_expr <- c
GSE81089_phe <- GSE81089_phe[colnames(c), ]
table(GSE81089_phe$OS)
GSE81089_phe <- GSE81089_phe[GSE81089_phe$OS != "n/a",]
GSE81089_expr = GSE81089_expr[ ,rownames(GSE81089_phe)]
save(GSE81089_expr,GSE81089_phe,file = "GSE81089.Rdata")



# 验证集
load(file = "GSE81089.Rdata")
load(file = "../Step4-1/Train_Result.Rdata")
gene_multi


Gene=rownames(gene_multi)
expr1 <- as.data.frame(t(GSE81089_expr[Gene,]))
phe1 <- GSE81089_phe
phe1$OS.time <- as.numeric(phe1$OS.time)
phe1$OS <- as.numeric(phe1$OS)
phe1 <- subset(phe1,phe1$OS.time != 0)
expr1 <- expr1[rownames(phe1),]


refGene <- c("ACTB","GAPDH","TFRC","TUBB")
refExpr <- as.data.frame(t(GSE81089_expr[refGene,rownames(phe1)])) 
Basal_CT <- cbind(refExpr,expr1)
for (i in 1:length(colnames(Basal_CT))) {
  Basal_CT[,i] <- as.numeric(Basal_CT[,i])
}
str(Basal_CT)
for (i in Gene) {
  Basal_CT[,i] = ((Basal_CT[,i] - Basal_CT$ACTB)+(Basal_CT[,i]-Basal_CT$GAPDH)+(Basal_CT[,i]-Basal_CT$TFRC)+(Basal_CT[,i]-Basal_CT$TUBB))/4
  
}
expr0 <- Basal_CT[,Gene]
expr0 <- scale(expr0)
expr0 <- as.data.frame(expr0)
sv3 <- cbind(expr0,phe1)
sv3$exp <- sv3[,rownames(gene_multi)[1]]*gene_multi[1,1]+sv3[,rownames(gene_multi)[2]]*gene_multi[2,1]+sv3[,rownames(gene_multi)[3]]*gene_multi[3,1]+sv3[,rownames(gene_multi)[4]]*gene_multi[4,1]+sv3[,rownames(gene_multi)[5]]*gene_multi[5,1]+sv3[,rownames(gene_multi)[6]]*gene_multi[6,1]
# train_sub_6 = sample(nrow(sv3),80)
data_6 <- sv3[train_sub_6,]
# train_sub_7 = sample(nrow(sv3),90)
data_7 <- sv3[train_sub_7,]
save(train_sub_6,data_6,train_sub_7,data_7,file = 'Verify_GEO80_90.Rdata')


load("Train_Result.Rdata")
load(file = 'Verify_GEO80_90.Rdata')
dat1<- data_6

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
       main="ROC of External Validation Set")
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


##### 生存分析
{
  
  
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
    dat1$OS.time=dat1$OS.time/365
    #最佳节点
    res.cut <- surv_cutpoint(dat1, #数据集
                             time = "OS.time", #生存时间
                             event = "OS", #生存状态
                             variables = "exp") #需要计算的数据列名
    
    summary(res.cut) #查看数据最佳截断点及统计量
    RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    # RiskGroup = ifelse(fp<median(fp),"Low","High")
    RiskGroup = factor(RiskGroup)
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
               title="Survival Curve for External Validation Set", 
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
    
    #mycolors <- colorRampPalette(c("#64b5f6", "#fffde7", "#ff5252"), bias = 1.2)(100)
    # mycolors <- colorRampPalette(c("#FF6666", "White", "#00CCCC"), bias = 1.2)(100)
    mycolors <- colorRampPalette(c(pal[1], "White", pal[2]), bias = 1.2)(100)
    tmp=t(scale(exp_dat))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,border_color=NA,
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

