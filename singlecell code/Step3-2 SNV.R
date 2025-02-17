

### Step3-2  NDRGs体细胞突变
### 整理时间： 2022/7/22
### 作者： 庞建宇

rm(list =ls())
options(stringsAsFactors = F)

setwd('/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step3-2/')
# BiocManager::install('maftools')
library(maftools)

# 加载突变数据
# luad.snv = read.delim(gzfile('/opt/disease/TCGA/LUAD/TCGA-LUAD.muse_snv.tsv.gz'))
# lusc.snv = read.delim(gzfile('/opt/disease/TCGA/LUSC/TCGA-LUSC.muse_snv.tsv.gz'))
lusc.snv = read.delim(gzfile('TCGA-LUSC.mutect2_snv.tsv.gz'))
luad.snv = read.delim(gzfile('TCGA-LUAD.mutect2_snv.tsv.gz'))

# 加载中性粒细胞分化相关基因
load(file = '../Step3-1/Neutrophil_Statemarker.Rdata')
gene = unique(Statemarker$gene) # 提取分化基因
load(file = '../Step2-0/TCGA_Step1 PheData.Rdata')

dat$Id2 = gsub('[.]','-',dat$Id)
rownames(dat) = dat$Id2
tmp = rbind(lusc.snv,luad.snv)
length(unique(tmp$gene))
length(unique(tmp$Sample_ID))
phe = subset(dat,dat$Id2 %in% unique(tmp$Sample_ID))
head(tmp)   
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(tmp$Variant_Type )
tmp2 = subset(tmp, tmp$Tumor_Sample_Barcode  %in% dat$Id2 ) # 添加临床信息
phe$age = ifelse(phe$Age <= 60,'<=60','>60')
tmp2$age = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$age, 'Unknow')
table(tmp2$age)
tmp2$Gender = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$Gender, 'Unknow')
table(tmp2$Gender)
tmp2$`Survive State` = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$OS, 'Unknow')
table(tmp2$`Survive State`)


neu.state = subset(tmp2,tmp2$Hugo_Symbol %in% gene)
# tcga.nsclc = read.maf(maf = tmp,
#                      vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))

tcga.neu.state = read.maf(maf = neu.state,
                      vc_nonSyn=names(tail(sort(table(neu.state$Variant_Classification )))))

# 添加临床信息
a = tcga.neu.state@clinical.data
a$age = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$age,'Unkown')
a$age[is.na(a$age)] = 'Unkown'
table(a$age)
a$Gender = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$Gender,'Unkown')
table(a$Gender)
a$`Survive State` = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$`Survive State`,'Unkown')
table(a$`Survive State`)
tcga.neu.state@clinical.data = a
# oncoplot(maf = tcga.nsclc) # 高频突变的前10个基因

oncoplot(maf = tcga.neu.state,top = 10) # 高频突变的前10个基因


# tcga.nsclc2 = read.maf('./SNP/fd9cb958-7f0f-4d33-b419-f161ffd292b9.wxs.aliquot_ensemble_masked.maf.gz')

getFields(tcga.neu.state)
getClinicalData(tcga.neu.state) #查看临床信息
getSampleSummary(tcga.neu.state)#查看每个样品发生突变的情况，此处就可以计算tumor mutation load,TML=Missense_Mutation/外显子数。


# plotmafSummary(maf = tcga.nsclc, rmOutlier = TRUE, 
#                addStat = 'median', dashboard = TRUE,
#                titvRaw=FALSE)#绘制整体的突变情况

plotmafSummary(maf = tcga.neu.state, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)#绘制整体的突变情况

# plotmafSummary(maf = tcga.nsclc2, rmOutlier = TRUE, 
#                addStat = 'median', dashboard = TRUE,
#                titvRaw=FALSE)#绘制整体的突变情况


#waterfall plot
#We will draw oncoplots for top ten mutated genes.
# oncoplot(maf = tcga.nsclc)
#绘制前20个突变基因的瀑布图。oncoplot()中参数gene=c()可以指定基因名，绘制感兴趣的基因的瀑布图

# oncoplot(maf = tcga.nsclc, top = 30)

# oncoplot(maf = tcga.neu.state, top = 24) # 展示突变率在5%以上的基因

library(ggsci)
col = pal_lancet()(6)
# col = RColorBrewer::brewer.pal(n = 6, name = 'Paired')
table(tcga.neu.state@data[["effect"]])
names(col) = c('frameshift_variant','synonymous_variant','missense_variant', 'intron_variant', '3_prime_UTR_variant' , 'stop_gained')

tcga.neu.state
getClinicalData(x = tcga.neu.state)

# anocol = pal_aaas()(3)
# names(anocol) = c("<=60", ">60", "Unkown")
# # anocol = list(age=anocol)
# 
# anocol2 = pal_nejm()(2)
# names(anocol2) = c("female", "male")
# # anocol = list(Gender=anocol2)
# 
# anocol3 = pal_d3()(2)
# names(anocol3) = c("Alive", "Dead")
# 
# fabcolors = list(age = anocol,Gender = anocol2, `Survive State` = anocol3)

pdf(file  = '../../Fig/Step3-2/neu_oncoplot.pdf',width = 10,height = 6)
oncoplot(maf = tcga.neu.state, top = 30, colors = col,
         draw_titv = T,
         clinicalFeatures = c('age','Gender','Survive State'),
         sortByAnnotation = TRUE,
         # annotationColor = c(anocol,anocol2,anocol3)
         )
dev.off()

# 转换 颠倒 
pdf(file  = '../../Fig/Step3-2/neu_titv.pdf',width = 10,height = 6)
titv(tcga.neu.state, useSyn = FALSE, plot = TRUE, file = NULL)
dev.off()

# 不同State突变情况
# State1 
table(Statemarker$cluster)
State1 = subset(Statemarker,Statemarker$cluster == 'State1')
table(unique(tmp2$Hugo_Symbol) %in% unique(State1$gene))
length(unique(State1$gene))
State1_pre = (276/303)*100
# neu.state1 = subset(tmp,tmp$Hugo_Symbol %in% State1$gene)
# table(neu.state1$effect)


# State2
table(Statemarker$cluster)
State2 = subset(Statemarker,Statemarker$cluster == 'State2')
table(unique(tmp2$Hugo_Symbol) %in% unique(State2$gene))
length(unique(State2$gene))
State2_pre = (141/152)*100
# neu.state2 = subset(tmp,tmp$Hugo_Symbol %in% State2$gene)


# State3
table(Statemarker$cluster)
State3 = subset(Statemarker,Statemarker$cluster == 'State3')
table(unique(tmp2$Hugo_Symbol) %in% unique(State3$gene))
length(unique(State3$gene))
State3_pre = (342/369)*100
# neu.state3 = subset(tmp,tmp$Hugo_Symbol %in% State3$gene)


# State4
table(Statemarker$cluster)
State4 = subset(Statemarker,Statemarker$cluster == 'State4')
table(unique(tmp2$Hugo_Symbol) %in% unique(State4$gene))
length(unique(State4$gene))
State4_pre = (307/334)*100
# neu.state4 = subset(tmp,tmp$Hugo_Symbol %in% State4$gene)


df = data.frame(Type = c('State1','State1','State2','State2','State3','State3','State4','State4'), 
                Mut = c('Mutation','Unmutated','Mutation','Unmutated','Mutation','Unmutated','Mutation','Unmutated'),
                Freq = c(276,303-276, 141,152-141, 342,369-342, 307,334-309),
                percent = c(State1_pre, 100-State1_pre, State2_pre, 100-State2_pre, State3_pre, 100-State3_pre, State4_pre, 100-State4_pre),
                lable = c(paste0(round(State1_pre, digits = 2),'%'),paste0(round(100-State1_pre, digits = 2),'%'),
                          paste0(round(State2_pre, digits = 2),'%'),paste0(round(100-State2_pre, digits = 2),'%'),
                          paste0(round(State3_pre, digits = 2),'%'),paste0(round(100-State3_pre, digits = 2),'%'),
                          paste0(round(State4_pre, digits = 2),'%'),paste0(round(100-State4_pre, digits = 2),'%')))

pvalue <- chisq.test(c(df$Freq,ncol=4))$p.value #卡方检验
library(plyr)
library(ggplot2)
ggplot(df,aes(Type,percent,fill=Mut))+
  geom_bar(stat="identity",position = position_stack())+
  # scale_color_igv()+
  scale_fill_manual(values = c("#DB423E","#008ECA"),label=c("Mutation","Unmutated"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="State",y="Percent Weidght",
       fill="")+
  geom_text(aes(label=lable),vjust=1.5,size=6,color="black")+
  annotate(geom = "text",
           cex=6,
           x=2.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("Chi-Squared Test, P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = 'neu_SNP.pdf',width = 8,height = 6,path = '../../Fig/Step3-2/')




# #计算每个基因出现的个数
# library(dplyr)
# neu.state2 = merge(lusc.snv,luad.snv) 
# neu.state2 = subset(neu.state2, gene %in% gene)
# 
# mut2 <- neu.state2 %>% filter(effect %in% c('frameshift_variant','missense_variant', 'intron_variant', '3_prime_UTR_variant' , 'stop_gained','synonymous_variant')) %>%
#   select(Sample_ID,gene) %>% 
#   group_by(gene) %>% 
#   summarise(Freq = n()) %>% 
#   arrange(desc(Freq))
# 
# head(mut2)
# 
# ####绘制基因词云#####
# library(wordcloud2)
# #绘制频次大于等于5的
# da <- subset(mut2,Freq >= 5) #、
# wordcloud2(da)
