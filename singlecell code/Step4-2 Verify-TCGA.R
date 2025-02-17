

### Step4-2  预后模型的验证
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step4-2/')


rm(list = ls())
options(stringsAsFactors = F)

load("../Step4-1/Train_Result.Rdata")


# 测试集 全集
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
load(file = "../Step2-0/TCGA_Step1 PheData.Rdata")
load(file = "../Step4-1/Ver_2.3_d.Rdata")

# 1/3 验证集
rownames(dat) <- dat$Id
Veriftphe <- dat[-train_sub_d,]
Veriftphe$OS.time <- ifelse(is.na(Veriftphe$`Days to Death`),Veriftphe$`Days to Last Follow`,Veriftphe$`Days to Death`)
table(is.na(Veriftphe$OS.time))
Veriftphe <- subset(Veriftphe, Veriftphe$OS.time != "NA")
table(Veriftphe$OS.time == 0)
Veriftphe <- subset(Veriftphe, Veriftphe$OS.time != 0)
table(Veriftphe$OS.time < 4000)
Veriftphe <- subset(Veriftphe, Veriftphe$OS.time <= 4000)
table(Veriftphe$OS)
Veriftphe$OS <- ifelse(Veriftphe$OS == "Alive",0,1)
table(Veriftphe$OS)
Veriftexpr <- as.data.frame(t(NSCLCcount[rownames(gene_multi),rownames(Veriftphe)]))
refGene <- c("ACTB","GAPDH","TFRC","TUBB")
refExpr <- as.data.frame(t(NSCLCcount[refGene,rownames(Veriftphe)]))
Basal_CT <- cbind(refExpr,Veriftexpr)
inter <- rownames(gene_multi)
for (i in inter) {
  Basal_CT[,i] = ((Basal_CT[,i] - Basal_CT$ACTB)+(Basal_CT[,i]-Basal_CT$GAPDH)+(Basal_CT[,i]-Basal_CT$TFRC)+(Basal_CT[,i]-Basal_CT$TUBB))/4
  
}
expr0 <- Basal_CT[,inter]
expr0 <- scale(expr0)
expr0 <- as.data.frame(expr0)
sv3 <- cbind(Veriftphe,expr0)
str(sv3)
sv3$exp <- sv3[,rownames(gene_multi)[1]]*gene_multi[1,1]+sv3[,rownames(gene_multi)[2]]*gene_multi[2,1]+sv3[,rownames(gene_multi)[3]]*gene_multi[3,1]+sv3[,rownames(gene_multi)[4]]*gene_multi[4,1]+sv3[,rownames(gene_multi)[5]]*gene_multi[5,1]+sv3[,rownames(gene_multi)[6]]*gene_multi[6,1]
train_sub_4 = sample(nrow(sv3),3/4*nrow(sv3))
data_4 <- sv3[train_sub_4,]
save(train_sub_4,data_4,file = 'Verify_TCGA283.Rdata')



# 验证
load("../Step4-1/Train_Result.Rdata")
load(file = "Verify_TCGA283.Rdata")

dat1<- data_4
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
       main="ROC of Internal Validation Set")
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
               title="Survival Curve for Internal Validation Set", 
               xlab ="Time(Years)",
               surv.median.line = "hv",
               ggtheme = theme_bw(base_size = 12))
  }
    
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
    table(RiskGroup)
    
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


