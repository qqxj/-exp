


### Step4-0  合并TCGA数据
### 整理时间： 2022/7/22
### 作者： 庞建宇
setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step2-0/')


rm(list=ls())
options(stringsAsFactors = F)

# 加载TCGA计数矩阵跟临床信息
LUADcount <- read.table(file = "TCGA-LUAD.htseq_counts.tsv",
                        header = T,row.names = 1)
LUADphe <- read.delim(file = "TCGA-LUAD.GDC_phenotype.tsv")

LUSCcount <- read.table(file = "TCGA-LUSC.htseq_counts.tsv",
                        header = T,row.names = 1)
LUSCphe <- read.delim(file = "TCGA-LUSC.GDC_phenotype.tsv")

dim(LUADcount)
dim(LUSCcount)

#合并LUAD跟LUSC
NSCLCcount <- cbind(LUADcount,LUSCcount)
NSCLCphe <- rbind(LUADphe,LUSCphe)
NSCLCcount$id <- rownames(NSCLCcount)


# ID转换
ids <- read.delim(file = "gencode.v22.annotation.gene.probeMap")
ids1 <- data.frame(ids$id,ids$gene) #取出ID跟Gene对应关系
names(ids1) <- c("id","gene")
NSCLCc = NSCLCcount
# NSCLCc$id <- NSCLCcount[,1] #取出Ensenmble Id这列
# NSCLCc <- NSCLCc[,-1]
head(ids1)
table(ids1$id %in% NSCLCc$Ensembl_ID)
#TRUE 
#60483 
library(dplyr)
NSCLCc <- inner_join(ids1,NSCLCc, by="id")
length(NSCLCc$gene)
length(unique(NSCLCc$gene)) # 有重复Gene
NSCLCc <- NSCLCc[!duplicated(NSCLCc$gene),] # 去除重复gene
rownames(NSCLCc) <- NSCLCc$gene
NSCLCc <- NSCLCc[,-c(1,2)]
NSCLCcount = NSCLCc

#除去表达矩阵中没有的临床样本
NSCLCphe$submitter_id.samples <- gsub("-",".",NSCLCphe$submitter_id.samples)#将-修改为.，与表达矩阵中的点与临床数据相对应，好进行取交集
NSCLCcount <- NSCLCcount[,colnames(NSCLCcount) %in% NSCLCphe$submitter_id.samples]#取交集,改变表达矩阵的列于临床数据对应
NSCLCphe <- NSCLCphe[NSCLCphe$submitter_id.samples %in% colnames(NSCLCcount),]#取交集改变临床的列与表达矩阵对应

# write.table(NSCLCcount,file = "NSCLCcount.txt")# 保存为txt 去跑CIBERSORT
# save(NSCLCcount,file = "NSCLC_FPKM.Rdata")

# 多数据整合，需批次校正
library(sva)
library(bladderbatch)

dat = NSCLCcount
class(dat)
# dat = as.matrix(dat)
# na.omit(dat)
batch = paste0("batch",rep(c(1,2),c(585,550)))
combat_edata <- ComBat(dat = dat, batch = batch)

NSCLCcount <- as.data.frame(combat_edata)



# 添加分组信息
rownames(NSCLCphe)<-NSCLCphe[,1]
group_text <- NSCLCphe
group_text<-group_text[c(-1)]   # #改变行名
pd <- data.frame(group_text$submitter_id,group_text$sample_type_id.samples)
table(pd$group_text.sample_type_id.samples)
pd[pd==2] <- 1
table(pd$group_text.sample_type_id.samples)

Group<-factor(pd$group_text.sample_type_id.samples,levels=c('1','11'))  ## 11 -normal,  1-toumal
design<-model.matrix(~0+Group)
colnames(design)<-c('Tumor','Normal')
rownames(design)<-rownames(group_text)

a = as.data.frame(colnames(NSCLCcount))
group_list = design[a$`colnames(NSCLCcount)`,]

library(stringr)
group = str_split(as.character(phe$sample_type.samples),' ',simplify = T)[,1]   #有两个是不一样的。

phe <- group_text
dat <- data.frame(rownames(phe),phe$age_at_diagnosis.diagnoses,phe$pathologic_M,phe$pathologic_N,phe$pathologic_T,phe$gender.demographic,phe$days_to_death.demographic,phe$days_to_last_follow_up.diagnoses,phe$vital_status.demographic,phe$tumor_stage.diagnoses)
names(dat) <- c("Id","Age","M","N","T","Gender","Days to Death","Days to Last Follow","OS","Stage")
table(dat$M)
dat[dat=="M1a"] <- "M1"
dat[dat=="M1b"] <- "M1"
table(dat$M)
table(dat$N)
dat[dat=="N3"] <- "N2"
dat[dat=="NX"] <- "N2"
table(dat$N)
table(dat$T)
dat[dat=="T1a"] <- "T1"
dat[dat=="T1b"] <- "T1"
dat[dat=="T2a"] <- "T2"
dat[dat=="T2b"] <- "T2"
table(dat$T)
table(dat$OS)
table(dat$Stage)
dat[dat=="stage i"] <- " stage i"
dat[dat=="stage ia"] <- " stage i"
dat[dat=="stage ib"] <- " stage i"
dat[dat=="stage ii"] <- " stage ii"
dat[dat=="stage iia"] <- " stage ii"
dat[dat=="stage iib"] <- " stage ii"
dat[dat=="stage iii"] <- " stage iii"
dat[dat=="stage iiia"] <- " stage iii"
dat[dat=="stage iiib"] <- " stage iii"
table(dat$Stage)
dat$Age <- phe$age_at_initial_pathologic_diagnosis
table(Group)
table(group)
group[group=="Recurrent"] <- "Primary"
dat$type <- ifelse(group=="Primary","Tumor","Normal")
table(dat$type)

save(dat,file = "TCGA_Step1 PheData.Rdata")
save(group_list,NSCLCcount,NSCLCphe,file = "TCGA_Step1 output.Rdata")

