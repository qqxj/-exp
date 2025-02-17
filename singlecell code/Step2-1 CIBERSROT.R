

### Step2-1  Bulk RNA-seq Analyze：CIBERSORT
### 整理时间： 2022/7/22
### 作者： 庞建宇


rm(list = ls())
options(stringsAsFactors = F)
gc()


setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step2-1/")

# NSCLC
# C0,C1,C10,C15.C19 : T Cell                                 Marker:"TRBC1","CD3D","CD3E"
# C2 : B cell                                                Marker:"CD79A","MS4A1","CD79B"
# C3 : Monocyte                                              Marker:"LST1","FCN1"
# C4,C5,C13,C20 : Dendritic Cell                             Marker:C1QB,C1QC,C1QA   C13:"MMP12","PLIN2","IGSF6"
# C6,C7,C8,C9,C12,C17,C23,C24,C26 : Epithelial/Cancer cell   Marker:KRT19,KRT18,KRT8
# C11,C21  Fibroblasts                                       Marker:"DCN","FBLN1","LUM","COL1A1" ,"COL3A1"
# C14 : Mast Cell                                            Marker:"TPSB2","TPSAB1","CPA3","CTSG"
# C16,C22  Plasma_cell                                       Marker:"IGLC2","IGHG1","IGHG3"
# C18,C27: Endothelial cell                                  Marker:"GNG11", "RAMP2" ,"VWF"
# C25  Oligodendrocytes                                      Marker:"PLP1", "MBP", "GPM6B", "APLP1"


load(file = "../Step1-1/res0.4_NSCLC.Integrate.Rdata")
library(dplyr)
library(Seurat)


##Step1   提取所有细胞Marker
table(Idents(NSCLC.Integrate))

FibroblastMarker <- FindMarkers(NSCLC.Integrate, ident.1 = "Fibroblast", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
FibroblastMarker_top10 <-FibroblastMarker  %>% top_n(n = 10, wt = avg_log2FC)

EndothelialMarker <- FindMarkers(NSCLC.Integrate, ident.1 = "Endothelial Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
EndothelialMarker_top10 <-EndothelialMarker  %>% top_n(n = 10, wt = avg_log2FC)

MastCellMarker <- FindMarkers(NSCLC.Integrate, ident.1 = "Mast Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
MastCellMarker_top10 <-MastCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

EpiCancerMarker <- FindMarkers(NSCLC.Integrate, ident.1 = "Epithelial/Cancer Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
EpiCancerMarker_top10 <-EpiCancerMarker  %>% top_n(n = 10, wt = avg_log2FC)

OligodendrocytesMarker <- FindMarkers(NSCLC.Integrate, ident.1 = "Oligodendrocytes", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
OligodendrocytesMarker_top10 <- OligodendrocytesMarker  %>% top_n(n = 10, wt = avg_log2FC)

rm(NSCLC.Integrate)
gc()


# TBcell-Sub
# C0,C9 : Natural killer (NK) Cell            Marker:"FGFBP2","GNLY","FCGR3A","XCL1","XCL2","KLRC1"
# C1,C5 : CD4+ T Memory Cell                  Marker:"NOSIP","CCR7"
# C2,C8 : CD8+ T Cell                         Marker:"GZMK","CD8A","CD8B"
# C3,C13 : B Cell                             Marker:"MS4A1","CD79A","VPREB3"
# C4,C11 : Natural killer T (NKT) Cell        Marker:"LMNA","RGCC","TUBA1B","TUBB","HIST1H4C","TYMS"
# C7 : Regulatory T (Treg) cell               Marker:"TIGIT","BATF","CORO1B"
# C10,C12 : Plasma Cell                       Marker:"IGLC2","IGKC","JCHAIN","IGHG4","IGHG1"
load(file = "../Step1-2/res0.4_TBcell.Rdata")
table(Idents(TBcell))
NKCellMarker <- FindMarkers(TBcell, ident.1 = "NK Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
NKCellMarker_top10 <-NKCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

CD4TMemoryCellMarker <- FindMarkers(TBcell, ident.1 = "CD4+ T Memory Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
CD4TMemoryCellMarker_top10 <-CD4TMemoryCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

CD8TCellMarker <- FindMarkers(TBcell, ident.1 = "CD8+ T Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
CD8TCellMarker_top10 <-CD8TCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

BCellMarker <- FindMarkers(TBcell, ident.1 = "B Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
BCellMarker_top10 <- BCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

NKTCellMarker <- FindMarkers(TBcell, ident.1 = "NKT Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
NKTCellMarker_top10 <- NKTCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

TregCellMarker <- FindMarkers(TBcell, ident.1 = "Regulatory T (Treg) Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
TregCellMarker_top10 <- TregCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

PlasmaCellMarker <- FindMarkers(TBcell, ident.1 = "Plasma Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
PlasmaCellMarker_top10 <- PlasmaCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

rm(TBcell)
gc()


# MonDC-Sub
#C0,C2,C5,C6,C7,C9,C11：Monocyte       Marker:"FBP1" , "MARCO" ,"RETN"
#C1,C12：Macrophage                    Marker:"FOLR2"   ,"PLTP" , "SEPP1" ,"CCL18"
#C3,C4：Dendritic Cell                 Marker:"CD48" , "FCN1"  , "IL1B" ,"LST1"
#C8：Granulocyte-Monocyte Progenitor   Marker:"KIAA0101" , "UBE2C"  , "PTTG1"
#10：Plasmacytoid Dendritic Cell       Marker:"GZMB" , "JCHAIN"  , "ITM2C","C12orf75"
#C13：Granulosa cell                   Marker:"BIRC3" , "CCL17"  , "CCL19"
#C14：Neutrophil                       Marker:"PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R"
load(file = "../Step1-3/res0.4_MonDC.Rdata")
table(Idents(MonDC))
MonocyteMarker <- FindMarkers(MonDC, ident.1 = "Monocyte", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
MonocyteMarker_top10 <- MonocyteMarker  %>% top_n(n = 10, wt = avg_log2FC)

MacrophageMarker <- FindMarkers(MonDC, ident.1 = "Macrophage", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
MacrophageMarker_top10 <- MacrophageMarker  %>% top_n(n = 10, wt = avg_log2FC)

DendriticCellMarker <- FindMarkers(MonDC, ident.1 = "Dendritic Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
DendriticCellMarker_top10 <- DendriticCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

GranulocyteMonocyteProgenitorMarker <- FindMarkers(MonDC, ident.1 = "Granulocyte-Monocyte Progenitor", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
GranulocyteMonocyteProgenitorMarker_top10 <- GranulocyteMonocyteProgenitorMarker  %>% top_n(n = 10, wt = avg_log2FC)

PlasmacytoidDendriticCellMarker <- FindMarkers(MonDC, ident.1 = "Plasmacytoid Dendritic Cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
PlasmacytoidDendriticCellMarker_top10 <- PlasmacytoidDendriticCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

GranulosaCellMarker <- FindMarkers(MonDC, ident.1 = "Granulosa cell", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
GranulosaCellMarker_top10 <- GranulosaCellMarker  %>% top_n(n = 10, wt = avg_log2FC)

NeutrophilMarker <- FindMarkers(MonDC, ident.1 = "Neutrophil", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
NeutrophilMarker_top10 <- NeutrophilMarker  %>% top_n(n = 10, wt = avg_log2FC)
# NeutrophilMarker <- FindMarkers(MonDC, ident.1 = "Neutrophil", only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)

# NeutrophilMarker <- FindMarkers(MonDC, ident.1 = "Neutrophil", only.pos = TRUE, logfc.threshold = 0.25)

rm(MonDC)
gc()


# Cancer Cell
load(file = "../Step1-4/res0.4_Epi.Rdata")
Idents(object = Epi) <- "copykat.pred"
CancerCellmarkers <- FindMarkers(Epi, ident.1 = "aneuploid", min.pct = 0.5,logfc.threshold = 0.5)
CancerCellmarkers_top10 <-CancerCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

rm(Epi)
gc()


# Normpl_Epi
#C0,C3,C4,C5,C7,C10,C11,C13：Basal Cell    Marker:"KRT17","GPX2",    "NTS","NIPSNAP2",    "TFPI2",   "NMU",   "ASS1", "BIRC3",     "LCN2","SLK","SELENOW",  "LTF"  , "DEFB1","FGB"
#C1：Pulmonary Alveolar Type II Cell       Marker:"PEBP4 ,"PLA2G1B" 
#C2：FOXN4+ Cell                           Marker:"CDKN2A" 
#C6：Luminal Epithelial Cell               Marker:"TCN1" 
#C8,C9：SLC16A7+ Cell                      Marker:"XAF1" ,"AGER"
#C12:Ionocyte Cell                         Marker:"PLCG2"
#C14:Langerhans Cell                       Marker:"RGS1","SRGN","CORO1A"
#C15:Unkonw
#C16:Ciliated Cell                         Marker:"C9orf24","C20orf85","C1orf194","FAM183A"
#C17,C18:Secretory Cell                    Marker:"ZG16B","SLC26A2"
load(file = "../Step1-4/res0.4_Normal_Epi.Rdata")

table(Idents(object = Normal_Epi))
BasalCellmarkers <- FindMarkers(Normal_Epi, ident.1 = "Basal Cell", min.pct = 0.5,logfc.threshold = 0.5)
BasalCellmarkers_top10 <-BasalCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

PulmonaryAlveolarTypeIICellmarkers <- FindMarkers(Normal_Epi, ident.1 = "Pulmonary Alveolar Type II Cell", min.pct = 0.5,logfc.threshold = 0.5)
PulmonaryAlveolarTypeIICellmarkers_top10 <-PulmonaryAlveolarTypeIICellmarkers %>% top_n(n = 10, wt = avg_log2FC)

FOXN4_Cellmarkers <- FindMarkers(Normal_Epi, ident.1 = "FOXN4+ Cell", min.pct = 0.5,logfc.threshold = 0.5)
FOXN4_Cellmarkers_top10 <-FOXN4_Cellmarkers %>% top_n(n = 10, wt = avg_log2FC)

LuminalEpithelialCellmarkers <- FindMarkers(Normal_Epi, ident.1 = "Luminal Epithelial Cell", min.pct = 0.5,logfc.threshold = 0.5)
LuminalEpithelialCellmarkers_top10 <-LuminalEpithelialCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

SLC16A7_ellmarkers <- FindMarkers(Normal_Epi, ident.1 = "SLC16A7+ Cell", min.pct = 0.5,logfc.threshold = 0.5)
SLC16A7_Cellmarkers_top10 <-SLC16A7_ellmarkers %>% top_n(n = 10, wt = avg_log2FC)

IonocyteCellmarkers <- FindMarkers(Normal_Epi, ident.1 = "SLC16A7+ Cell", min.pct = 0.5,logfc.threshold = 0.5)
IonocyteCellmarkers_top10 <-IonocyteCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

LangerhansCellmarkers <- FindMarkers(Normal_Epi, ident.1 = "Langerhans Cell", min.pct = 0.5,logfc.threshold = 0.5)
LangerhansCellmarkers_top10 <-LangerhansCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

CiliatedCellmarkers <- FindMarkers(Normal_Epi, ident.1 = "Ciliated Cell", min.pct = 0.5,logfc.threshold = 0.5)
CiliatedCellmarkers_top10 <-CiliatedCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

SecretoryCellmarkers <- FindMarkers(Normal_Epi, ident.1 = "Secretory Cell", min.pct = 0.5,logfc.threshold = 0.5)
SecretoryCellmarkers_top10 <-SecretoryCellmarkers %>% top_n(n = 10, wt = avg_log2FC)

rm(Normal_Epi)
gc()


# load(file = "res0.4_AllMarkergene.Rdata")
# load(file = "res0.4_All-Top10-Markergene.Rdata")
save(CancerCellmarkers,FibroblastMarker,EndothelialMarker,MastCellMarker,NKCellMarker,CD4TMemoryCellMarker,
     CD8TCellMarker,BCellMarker,NKTCellMarker,TregCellMarker,PlasmaCellMarker,MonocyteMarker,
     MacrophageMarker,DendriticCellMarker,GranulocyteMonocyteProgenitorMarker,PlasmacytoidDendriticCellMarker,
     GranulosaCellMarker,NeutrophilMarker,BasalCellmarkers,PulmonaryAlveolarTypeIICellmarkers,FOXN4_Cellmarkers,
     LuminalEpithelialCellmarkers,CiliatedCellmarkers,SecretoryCellmarkers,SLC16A7_ellmarkers,IonocyteCellmarkers,LangerhansCellmarkers,OligodendrocytesMarker,
     file = "res0.4_AllMarkergene.Rdata")

save(CancerCellmarkers_top10,FibroblastMarker_top10,EndothelialMarker_top10,MastCellMarker_top10,NKCellMarker_top10,CD4TMemoryCellMarker_top10,
     CD8TCellMarker_top10,BCellMarker_top10,NKTCellMarker_top10,TregCellMarker_top10,PlasmaCellMarker_top10,MonocyteMarker_top10,
     MacrophageMarker_top10,DendriticCellMarker_top10,GranulocyteMonocyteProgenitorMarker_top10,PlasmacytoidDendriticCellMarker_top10,
     GranulosaCellMarker_top10,NeutrophilMarker_top10,BasalCellmarkers_top10,PulmonaryAlveolarTypeIICellmarkers_top10,FOXN4_Cellmarkers_top10,
     LuminalEpithelialCellmarkers_top10,CiliatedCellmarkers_top10,SecretoryCellmarkers_top10,SLC16A7_Cellmarkers_top10,IonocyteCellmarkers_top10,LangerhansCellmarkers_top10,OligodendrocytesMarker_top10,
     file = "res0.4_All-Top10-Markergene.Rdata")




##Step2  提取细胞总表达量
rm(list = ls())
options(stringsAsFactors = F)
gc()

# 提取单独的细胞簇基因的表达量之和
library(Seurat)

# TCcell
load(file = "../Step1-2/res0.4_TBcell.Rdata")

table(Idents(TBcell))
PlasmaCell  = TBcell[,TBcell@active.ident %in% "Plasma Cell"] #提取Plasma Cell簇
x = as.matrix(PlasmaCell@assays[["RNA"]]@data) # 提取基因表达数据
y = rowSums(x) # 计算该簇中所有基因的总和
head(y)
class(y)
y = as.data.frame(y)
colnames(y) = "Plasma Cell"


TregCell  = TBcell[,TBcell@active.ident %in% "Regulatory T (Treg) Cell"] #提取TregCell簇
x1 = as.matrix(TregCell@assays[["RNA"]]@data) # 提取基因表达数据
y1 = rowSums(x1) # 计算该簇中所有基因的总和
head(y1)
class(y1)
y$TregCell = y1


NKTCell  = TBcell[,TBcell@active.ident %in% "NKT Cell"] #提取NKTCell簇
x2 = as.matrix(TregCell@assays[["RNA"]]@data) # 提取基因表达数据
y2 = rowSums(x2) # 计算该簇中所有基因的总和
head(y2)
y$NKTCell = y2


BCell = TBcell[,TBcell@active.ident %in% "B Cell"] #提取B Cell簇
x3 = as.matrix(BCell@assays[["RNA"]]@data) # 提取基因表达数据
y3 = rowSums(x3) # 计算该簇中所有基因的总和
head(y3)
y$BCell = y3


CD8TCell  = TBcell[,TBcell@active.ident %in% "CD8+ T Cell"] #提取CD8+ T Cell簇
x4 = as.matrix(CD8TCell@assays[["RNA"]]@data) # 提取基因表达数据
y4 = rowSums(x4) # 计算该簇中所有基因的总和
head(y4)
y$CD8TCell = y4


CD4TMemoryCell  = TBcell[,TBcell@active.ident %in% "CD4+ T Memory Cell"] #提取CD4+ T Memory Cell簇
x5 = as.matrix(CD4TMemoryCell@assays[["RNA"]]@data) # 提取基因表达数据
y5 = rowSums(x5) # 计算该簇中所有基因的总和
head(y5)
y$CD4TMemoryCell = y5


NKCell  = TBcell[,TBcell@active.ident %in% "NK Cell"] #提取NK Cell簇
x6 = as.matrix(NKCell@assays[["RNA"]]@data) # 提取基因表达数据
y6 = rowSums(x6) # 计算该簇中所有基因的总和
head(y6)
y$NKCell = y6

Genedata = y
colnames(Genedata)
names(Genedata) = c("Plasma Cell","Treg Cell","NKT Cell","B Cell","CD8+ T Cell","CD4+ T Memory Cell", "NK Cell")
save(Genedata, file = "Genedata.Rdata")




# 总细胞类群
rm(list = ls())
options(stringsAsFactors = F)
gc()
load(file = "../Step1-1/res0.4_NSCLC.Integrate.Rdata")
load(file = "Genedata.Rdata")

# NSCLC
# C0,C1,C10,C15.C19 : T Cell                                 Marker:"TRBC1","CD3D","CD3E"
# C2 : B cell                                                Marker:"CD79A","MS4A1","CD79B"
# C3 : Monocyte                                              Marker:"LST1","FCN1"
# C4,C5,C13,C20 : Dendritic Cell                             Marker:C1QB,C1QC,C1QA   C13:"MMP12","PLIN2","IGSF6"
# C6,C7,C8,C9,C12,C17,C23,C24,C26 : Epithelial/Cancer cell   Marker:KRT19,KRT18,KRT8
# C11,C21  Fibroblasts                                       Marker:"DCN","FBLN1","LUM","COL1A1" ,"COL3A1"
# C14 : Mast Cell                                            Marker:"TPSB2","TPSAB1","CPA3","CTSG"
# C16,C22  Plasma_cell                                       Marker:"IGLC2","IGHG1","IGHG3"
# C18,C27: Endothelial cell                                  Marker:"GNG11", "RAMP2" ,"VWF"
# C25  Oligodendrocytes                                      Marker:"PLP1", "MBP", "GPM6B", "APLP1"

table(Idents(NSCLC.Integrate))
Fibroblast  = NSCLC.Integrate[,NSCLC.Integrate@active.ident %in% "Fibroblast"] #提取Fibroblast簇
x = as.matrix(Fibroblast@assays[["RNA"]]@data) # 提取基因表达数据
y = rowSums(x) # 计算该簇中所有基因的总和
head(y)
Genedata$Fibroblast = y


MastCell  = NSCLC.Integrate[,NSCLC.Integrate@active.ident %in% "Mast Cell"] #提取Mast Cell簇
x1 = as.matrix(MastCell@assays[["RNA"]]@data) # 提取基因表达数据
y1 = rowSums(x1) # 计算该簇中所有基因的总和
head(y1)
Genedata$MastCell = y1


EndothelialCell  = NSCLC.Integrate[,NSCLC.Integrate@active.ident %in% "Endothelial Cell"] #提取Endothelial Cell簇
x2 = as.matrix(EndothelialCell@assays[["RNA"]]@data) # 提取基因表达数据
y2 = rowSums(x2) # 计算该簇中所有基因的总和
head(y2)
Genedata$EndothelialCell = y2


Oligodendrocytes  = NSCLC.Integrate[,NSCLC.Integrate@active.ident %in% "Oligodendrocytes"] #提取Oligodendrocytes簇
x3 = as.matrix(Oligodendrocytes@assays[["RNA"]]@data) # 提取基因表达数据
y3 = rowSums(x3) # 计算该簇中所有基因的总和
head(y3)
Genedata$Oligodendrocytes = y3

names(Genedata) = c("Plasma Cell","Treg Cell","NKT Cell","B Cell","CD8+ T Cell","CD4+ T Memory Cell", "NK Cell",
                    "Fibroblast","Mast Cell","Endothelial Cell")
save(Genedata, file = "Genedata.Rdata")




# MonDC
rm(list = ls())
options(stringsAsFactors = F)
gc()
#C0,C2,C5,C6,C7,C9,C11：Monocyte       Marker:"FBP1" , "MARCO" ,"RETN"
#C1,C12：Macrophage                    Marker:"FOLR2"   ,"PLTP" , "SEPP1" ,"CCL18"
#C3,C4：Dendritic Cell                 Marker:"CD48" , "FCN1"  , "IL1B" ,"LST1"
#C8：Granulocyte-Monocyte Progenitor   Marker:"KIAA0101" , "UBE2C"  , "PTTG1"
#10：Plasmacytoid Dendritic Cell       Marker:"GZMB" , "JCHAIN"  , "ITM2C","C12orf75"
#C13：Granulosa cell                   Marker:"BIRC3" , "CCL17"  , "CCL19"
#C14：Neutrophil                       Marker:"PI3" , "IL1R2"  , "LRG1","SLC25A37","CSF3R"
load(file = "Genedata.Rdata")
load(file = "../Step1-3/res0.4_MonDC.Rdata")

table(Idents(MonDC))
Neutrophil  = MonDC[,MonDC@active.ident %in% "Neutrophil"] #提取Neutrophil簇
x1 = as.matrix(Neutrophil@assays[["RNA"]]@data) # 提取基因表达数据
y1 = rowSums(x1) # 计算该簇中所有基因的总和
head(y1)
Genedata$Neutrophil = y1

GranulosaCell  = MonDC[,MonDC@active.ident %in% "Granulosa cell"] #提取Granulosa cell簇
x2 = as.matrix(GranulosaCell@assays[["RNA"]]@data) # 提取基因表达数据
y2 = rowSums(x2) # 计算该簇中所有基因的总和
head(y2)
Genedata$GranulosaCell = y2

PlasmacytoidDendriticCell  = MonDC[,MonDC@active.ident %in% "Plasmacytoid Dendritic Cell"] #提取Plasmacytoid Dendritic Cell簇
x3 = as.matrix(PlasmacytoidDendriticCell@assays[["RNA"]]@data) # 提取基因表达数据
y3 = rowSums(x3) # 计算该簇中所有基因的总和
head(y3)
Genedata$PlasmacytoidDendriticCell = y3

GranulocyteMonocyteProgenitor  = MonDC[,MonDC@active.ident %in% "Granulocyte-Monocyte Progenitor"] #提取Granulocyte-Monocyte Progenitor簇
x4 = as.matrix(GranulocyteMonocyteProgenitor@assays[["RNA"]]@data) # 提取基因表达数据
y4 = rowSums(x4) # 计算该簇中所有基因的总和
head(y4)
Genedata$GranulocyteMonocyteProgenitor = y4

DendriticCell  = MonDC[,MonDC@active.ident %in% "Dendritic Cell"] #提取Dendritic Cell簇
x5 = as.matrix(DendriticCell@assays[["RNA"]]@data) # 提取基因表达数据
y5 = rowSums(x5) # 计算该簇中所有基因的总和
head(y5)
Genedata$DendriticCell = y5

Macrophage  = MonDC[,MonDC@active.ident %in% "Macrophage"] #提取Macrophage簇
x6 = as.matrix(Macrophage@assays[["RNA"]]@data) # 提取基因表达数据
y6 = rowSums(x6) # 计算该簇中所有基因的总和
head(y6)
Genedata$Macrophage = y6

Monocyte  = MonDC[,MonDC@active.ident %in% "Monocyte"] #提取Monocyte簇
x7 = as.matrix(Monocyte@assays[["RNA"]]@data) # 提取基因表达数据
y7 = rowSums(x7) # 计算该簇中所有基因的总和
head(y7)
Genedata$Monocyte = y7

names(Genedata) = c("Plasma Cell","Treg Cell","NKT Cell","B Cell","CD8+ T Cell","CD4+ T Memory Cell", "NK Cell",
                    "Fibroblast","Mast Cell","Endothelial Cell","Neutrophil","Granulosa Cell","Plasmacytoid Dendritic Cell",
                    "Granulocyte-Monocyte Progenitor","Dendritic Cell","Macrophage","Monocyte")
save(Genedata, file = "Genedata.Rdata")




# 上皮细胞类群
rm(list = ls())
options(stringsAsFactors = F)
gc()
load(file = "Genedata.Rdata")
load(file = "../Step1-4/res0.4_Epi.Rdata.Rdata")
CancerCell  = Epi[,Epi@meta.data$copykat.pred %in% "aneuploid"] #提取CancerCell簇
x = as.matrix(CancerCell@assays[["RNA"]]@data) # 提取基因表达数据
y = rowSums(x) # 计算该簇中所有基因的总和
head(y)
class(y)
Genedata$CancerCell = y
save(Genedata, file = "Genedata.Rdata")



#Noemal_Epi
load(file = "../Step1-4/res0.4_Normal_Epi.Rdata")

table(Idents(Normal_Epi))
BasalCell  = Normal_Epi[,Normal_Epi@active.ident %in% "Basal Cell"] #提取Basal Cell簇
x1 = as.matrix(BasalCell@assays[["RNA"]]@data) # 提取基因表达数据
y1 = rowSums(x1) # 计算该簇中所有基因的总和
head(y1)
Genedata$BasalCell = y1

PulmonaryAlveolarTypeIICell  = Normal_Epi[,Normal_Epi@active.ident %in% "Pulmonary Alveolar Type II Cell"] #提取Pulmonary Alveolar Type II Cell簇
x2 = as.matrix(PulmonaryAlveolarTypeIICell@assays[["RNA"]]@data) # 提取基因表达数据
y2 = rowSums(x2) # 计算该簇中所有基因的总和
head(y2)
Genedata$PulmonaryAlveolarTypeIICell = y2

FOXN4_Cell  = Normal_Epi[,Normal_Epi@active.ident %in% "FOXN4+ Cell"] #提取FOXN4+ Cell簇
x3 = as.matrix(FOXN4_Cell@assays[["RNA"]]@data) # 提取基因表达数据
y3 = rowSums(x3) # 计算该簇中所有基因的总和
head(y3)
Genedata$FOXN4_Cell = y3

LuminalEpithelialCell  = Normal_Epi[,Normal_Epi@active.ident %in% "Luminal Epithelial Cell"] #提取Luminal Epithelial Cell簇
x4 = as.matrix(LuminalEpithelialCell@assays[["RNA"]]@data) # 提取基因表达数据
y4 = rowSums(x4) # 计算该簇中所有基因的总和
head(y4)
Genedata$LuminalEpithelialCell = y4

SLC16A7_Cell  = Normal_Epi[,Normal_Epi@active.ident %in% "SLC16A7+ Cell"] #提取SLC16A7+ Cell簇
x5 = as.matrix(SLC16A7_Cell@assays[["RNA"]]@data) # 提取基因表达数据
y5 = rowSums(x5) # 计算该簇中所有基因的总和
head(y5)
Genedata$SLC16A7_Cell = y5

IonocyteCell  = Normal_Epi[,Normal_Epi@active.ident %in% "Ionocyte Cell"] #提取Ionocyte Cell簇
x6 = as.matrix(IonocyteCell@assays[["RNA"]]@data) # 提取基因表达数据
y6 = rowSums(x6) # 计算该簇中所有基因的总和
head(y6)
Genedata$IonocyteCell = y6

LangerhansCell  = Normal_Epi[,Normal_Epi@active.ident %in% "Langerhans Cell"] #提Langerhans Cell簇
x7 = as.matrix(LangerhansCell@assays[["RNA"]]@data) # 提取基因表达数据
y7 = rowSums(x7) # 计算该簇中所有基因的总和
head(y7)
Genedata$LangerhansCell = y7

CiliatedCell  = Normal_Epi[,Normal_Epi@active.ident %in% "Ciliated Cell"] #提Ciliated Cell簇
x8 = as.matrix(CiliatedCell@assays[["RNA"]]@data) # 提取基因表达数据
y8 = rowSums(x8) # 计算该簇中所有基因的总和
head(y8)
Genedata$CiliatedCell = y8

SecretoryCell  = Normal_Epi[,Normal_Epi@active.ident %in% "Secretory Cell"] #提Secretory Cell簇
x9 = as.matrix(SecretoryCell@assays[["RNA"]]@data) # 提取基因表达数据
y9 = rowSums(x9) # 计算该簇中所有基因的总和
head(y9)
Genedata$SecretoryCell = y9

colnames(Genedata)
names(Genedata) = c("Plasma Cell","Treg Cell","NKT Cell","B Cell","CD8+ T Cell","CD4+ T Memory Cell", "NK Cell",
                    "Fibroblast","Mast Cell","Endothelial Cell","Neutrophil","Granulosa Cell","Plasmacytoid Dendritic Cell",
                    "Granulocyte-Monocyte Progenitor","Dendritic Cell","Macrophage","Monocyte","Cancer Cell","Basal Cell",
                    "Pulmonary Alveolar TypeII Cell","FOXN4+ Cell", "Luminal Epithelial Cell", "SLC16A7+ Cell", "Ionocyte Cell",
                    "Langerhans Cell", "Ciliated Cell", "Secretory Cell")
save(Genedata, file = "Genedata.Rdata")




##Step3  CIBERSORT
rm(list = ls())
options(stringsAsFactors = F)
setwd("/home/pjy/NSCLC/NSCLC/Rproject/Rdata/Step2-1/")

# 表达数据处理
library(limma)
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
exp = NSCLCcount
# exp1 = NSCLCcount
# 整理一下行名，列名，删除表达特别低的基因
exp=as.matrix(exp)
# rownames(exp)=exp[,1]
# exp=exp[,2:ncol(exp)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

# 去除负值
data[data<0] = NA
data = na.omit(data)
# exp1[exp1<0] = NA
# exp1 = na.omit(exp1)
is.recursive(data)
is.atomic(data)

#把准备输入CIBERSORT的数据保存一下
# v <-voom(data, plot = F, save.plot = F) # log2转换
# out=v$E
# out=rbind(ID=colnames(out),out)
write.table(data,file="TCGA_1135_ready.txt",sep = "\t") 

# write.table(out,file="TCGA_100_ready.txt",sep="\t",quote=F,col.names=F)



# 细胞注释信息处理
load(file = "res0.4_AllMarkergene.Rdata")
load(file = "Genedata.Rdata")


y1 = rownames(BCellMarker)
y2 = rownames(CD4TMemoryCellMarker)
y3 = rownames(CD8TCellMarker)
y4 = rownames(NKCellMarker)
y5 = rownames(NKTCellMarker)
y6 = rownames(PlasmaCellMarker)
y7 = rownames(TregCellMarker)

y8 = rownames(EndothelialMarker)
y9 = rownames(FibroblastMarker)
y10 = rownames(MastCellMarker)

y11 = rownames(DendriticCellMarker)
y12 = rownames(GranulocyteMonocyteProgenitorMarker)
y13 = rownames(GranulosaCellMarker)
y14 = rownames(MacrophageMarker)
y15 = rownames(MonocyteMarker)
y16 = rownames(NeutrophilMarker)
y17 = rownames(PlasmacytoidDendriticCellMarker)

y18 = rownames(FOXN4_Cellmarkers)
y19 = rownames(BasalCellmarkers)
y20 = rownames(CancerCellmarkers)
y21 = rownames(CiliatedCellmarkers)
y22 = rownames(LuminalEpithelialCellmarkers)
y23 = rownames(SecretoryCellmarkers)
y24 = rownames(PulmonaryAlveolarTypeIICellmarkers)
y25 = rownames(SLC16A7_ellmarkers)
y26 = rownames(IonocyteCellmarkers)
y27 = rownames(LangerhansCellmarkers)

y28 = rownames(OligodendrocytesMarker)

x = list(c(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,
           y20,y21,y22,y23,y24,y25,y26,y27,y28))

MarkerGene = Genedata[x[[1]],]
colnames(MarkerGene)
MarkerGene$`Gene Symble` <- rownames(MarkerGene)
MarkerGene <- MarkerGene[,c(29,1:28)]

# names(MarkerGene) = c("Gene Symble","Treg Cell","NKT Cell","B Cell","CD8+ T Cell","CD4+ T Memory Cell", "NK Cell",
#                       "Fibroblast","Mast Cell","Endothelial Cell","Neutrophil","Granulosa Cell","Plasmacytoid Dendritic Cell",
#                       "Granulocyte-Monocyte Progenitor","Dendritic Cell","Macrophage","Monocyte","Cancer Cell","Basal Cell",
#                       "Pulmonary Alveolar TypeII Cell","FOXN4+ Cell", "Luminal Epithelial Cell", "SLC16A7+ Cell", "Ionocyte Cell",
#                       "Langerhans Cell", "Ciliated Cell", "Secretory Cell","Plasma Cell")
# 


length(MarkerGene$`Gene Symble`)
length(unique(MarkerGene$`Gene Symble`)) # 有重复Gene
MarkerGene <- MarkerGene[!duplicated(MarkerGene$`Gene Symble`),] # 去除重复gene
rownames(MarkerGene) = NULL

save(MarkerGene,file = "MarGene.Rdata")
write.table(MarkerGene,file="MarGene.txt",sep = "\t")# EXCEL打开，删除第一列


#top10
# y1 = rownames(BCellMarker_top10)
# y2 = rownames(CD4TMemoryCellMarker_top10)
# y3 = rownames(CD8TCellMarker_top10)
# y4 = rownames(NKCellMarker_top10)
# y5 = rownames(NKTCellMarker_top10)
# y6 = rownames(PlasmaCellMarker_top10)
# y7 = rownames(TregCellMarker_top10)
# 
# y8 = rownames(EndothelialMarker_top10)
# y9 = rownames(FibroblastMarker_top10)
# y10 = rownames(MastCellMarker_top10)
# 
# y11 = rownames(DendriticCellMarker_top10)
# y12 = rownames(GranulocyteMonocyteProgenitorMarker_top10)
# y13 = rownames(GranulosaCellMarker_top10)
# y14 = rownames(MacrophageMarker_top10)
# y15 = rownames(MonocyteMarker_top10)
# y16 = rownames(NeutrophilMarker_top10)
# y17 = rownames(PlasmacytoidDendriticCellMarker_top10)
# 
# y18 = rownames(FOXN4_Cellmarkers_top10)
# y19 = rownames(BasalCellmarkers_top10)
# y20 = rownames(CancerCellmarkers_top10)
# y21 = rownames(CiliatedCellmarkers_top10)
# y22 = rownames(LuminalEpithelialCellmarkers_top10)
# y23 = rownames(SecretoryCellmarkers_top10)
# y24 = rownames(PulmonaryAlveolarTypeIICellmarkers_top10)
# y25 = rownames(SLC16A7_Cellmarkers_top10)
# y26 = rownames(IonocyteCellmarkers_top10)
# y27 = rownames(LangerhansCellmarkers_top10)
# 
# y28 = rownames(OligodendrocytesMarker_top10)
# 
# x = list(c(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,
#            y20,y21,y22,y23,y24,y25,y26,y27,y28))
# 
# MarkerGene_top10 = Genedata[x[[1]],]
# colnames(MarkerGene_top10)
# MarkerGene_top10$`Gene Symble` <- rownames(MarkerGene_top10)
# MarkerGene_top10 <- MarkerGene_top10[,c(29,1:28)]
# 
# length(MarkerGene_top10$`Gene Symble`)
# length(unique(MarkerGene_top10$`Gene Symble`)) # 有重复Gene
# MarkerGene_top10 <- MarkerGene_top10[!duplicated(MarkerGene_top10$`Gene Symble`),] # 去除重复gene
# rownames(MarkerGene_top10) = NULL
# 
# save(MarkerGene_top10,file = "MarGene_top10.Rdata")
# write.table(MarkerGene_top10,file="MarGene_top10.txt",sep = "\t")# EXCEL打开，删除第一列
# 
# 
# load(file = "Genedata_Before.Rdata")
# a = list(c(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y28))
# 
# MarkerGene2 = Genedata_Before[a[[1]],]
# colnames(MarkerGene2)
# MarkerGene2$PlasmaCell = MarkerGene2[,1]
# MarkerGene2$`Plasma Cell` = rownames(MarkerGene2)
# names(MarkerGene2) = c("Gene Symble","Treg Cell","NKT Cell","B Cell","CD8+ T Cell","CD4+ T Memory Cell", "NK Cell",
#                       "Fibroblast","Mast Cell","Endothelial Cell","Neutrophil","Granulosa Cell","Plasmacytoid Dendritic Cell",
#                       "Granulocyte-Monocyte Progenitor","Dendritic Cell","Macrophage","Monocyte","Epithelium/Cancer Cell","Plasma Cell")
# rownames(MarkerGene2) = NULL
# 
# 
# 
# save(MarkerGene,file = "MarGene.Rdata")
# save(MarkerGene2,file = "MarGene_Before.Rdata")
# write.table(MarkerGene,file="MarGene.txt",sep = "\t")# EXCEL打开，删除第一列
# write.table(MarkerGene2,file="MarGene_Before.txt",sep = "\t")# EXCEL打开，删除第一列


# top 10
# load(file = 'res0.4_NSCLC.Integrate_27cell.Rdata')
# Marker <- FindAllMarkers(datSet, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5) # 运行时间较长
# Marker_top10 <-Marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# save(Marker,Marker_top10,file = 'res0.4_NSCLC.Integrate_27cell_Marker.Rdata')
# 
# 
# load(file = 'res0.4_NSCLC.Integrate_27cell_Marker.Rdata')
# load(file = "Genedata.Rdata")
# # load(file = "MarGene.Rdata")
# 
# # Marker_data <-  Genedata[Marker_top10$gene,]
# Marker_data <-  Genedata[Marker$gene,]
# Marker_data$symble <- rownames(Marker_data)
# Marker_data <- Marker_data[,c(28,1:27)]
# rownames(Marker_data) <- NULL
# # write.table(Marker_data,file="top10MarGenedata.txt",sep = "\t")# EXCEL打开，删除第一列
# write.table(Marker_data,file="MarGenedata.txt",sep = "\t")# EXCEL打开，删除第一列

# 运行Cibersort
source("Cibersort.R")
library(preprocessCore)
library(e1071)
library(parallel)


# result1 <- CIBERSORT("LM22.txt", "exp.txt", perm = 1000, QN = T) # 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result1,file = "Cibersort.Rdata")
# 
# result2 <- CIBERSORT("LM22.txt", "TCGA_1135_ready.txt", perm = 1000, QN = T) # 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result2,file = "Cibersort2.Rdata")

# result3 <- CIBERSORT("MarGene.txt", "exp.txt", perm = 1000, QN = T)# 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result3,file = "Cibersort3.Rdata")
# 
# result4 <- CIBERSORT("MarGene.txt", "TCGA_1135_ready.txt", perm = 1000, QN = T)# 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result4,file = "Cibersort_Final.Rdata")
# write.table(result4,file = "Cibersort_Final.txt",sep = "\t")


# result5 <- CIBERSORT("MarGene_Before.txt", "TCGA_1135_ready.txt", perm = 1000, QN = T)# 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result5,file = "Cibersort_Before.Rdata")
# write.table(result5,file = "Cibersort_Before.txt",sep = "\t")


# result_top10 <- CIBERSORT("28cellMarGene_top10.txt", "TCGA_1135_ready.txt", perm = 1000, QN = T)# 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result_top10,file = "Cibersort_top10.Rdata")
# write.table(result_top10,file = "Cibersort_top10.txt",sep = "\t")


# result_Allmarker <- CIBERSORT("28cellMarGene.txt", "TCGA_1135_ready.txt", perm = 1000, QN = T)# 正常运行
# # perm置换次数=1000，QN分位数归一化=TRUE
# # 文件名可以自定义
# save(result_Allmarker,file = "Cibersort_Allmarker.Rdata")
# write.table(result_Allmarker,file = "Cibersort_Allmarker.txt",sep = "\t")


result_Allmarker <- CIBERSORT("MarGene.txt", "TCGA_1135_ready.txt", perm = 1000, QN = T)# 正常运行
# perm置换次数=1000，QN分位数归一化=TRUE
# 文件名可以自定义
save(result_Allmarker,file = "Cibersort_Allmarker.Rdata")
write.table(result_Allmarker,file = "Cibersort_Allmarker.txt",sep = "\t")


# # 结果可视化
# # 热图
# load(file = "Cibersort_Final.Rdata")
# # load(file = "Cibersort_top10.Rdata")
# load(file = "../TCGA-NSCLC-Cibersort/Step1 output.Rdata")
# re <- result_top10[,-(29:31)]
# library(pheatmap)
# k <- apply(re,2,function(x) {sum(x == 0) < nrow(result_top10)/2})
# table(k)
# re2 <- as.data.frame(t(re[,k]))
# an = data.frame(group = group_list,
#                 row.names = colnames(NSCLCcount))
# png(filename = "top10_pheatmap.png",width = 800,height = 1200)
# pheatmap(re2,scale = "row",
#          show_colnames = F,
#          # annotation_col = group_list,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# dev.off()
# 
# # 直方图
# library(RColorBrewer)
# library(dplyr)
# library(tidyr)
# library(tibble)
# library(ggplot2)
# mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
# 
# dat <- re %>% as.data.frame() %>%
#   rownames_to_column("Sample") %>%
#   gather(key = Cell_type,value = Proportion,-Sample)
# 
# 
# ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) +
#   geom_bar(stat = "identity") +
#   labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom") +
#   scale_y_continuous(expand = c(0.01,0)) +
#   scale_fill_manual(values = mypalette(28))
# ggsave(filename = "top10_Est_Celltype.png",width = 16,height = 10)
# 
# 
# 
# 
# # 箱线图
# 
# ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
#   geom_boxplot(outlier.shape = 21,color = "black") +
#   theme_bw() +
#   labs(x = "Cell Type", y = "Estimated Proportion") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom") +
#   scale_fill_manual(values = mypalette(28))
# ggsave(filename = "top10_Est_Celltype(2).png",width = 16,height = 10)
# 
# 
# 
# # 让箱线图有顺序
# a = dat %>%
#   group_by(Cell_type) %>%
#   summarise(m = median(Proportion)) %>%
#   arrange(desc(m)) %>%
#   pull(Cell_type)
# 
# dat$Cell_type = factor(dat$Cell_type,levels = a)
# 
# 
# ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
#   geom_boxplot(outlier.shape = 21,color = "black") +
#   theme_bw() +
#   labs(x = "Cell Type", y = "Estimated Proportion") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "bottom") +
#   scale_fill_manual(values = mypalette(28))
# 
# ggsave(filename = "top10_Est_Celltype(3).png",width = 16,height = 10)
# 
# 
# 
# # 根据分组画箱图
# library(stringr)
# dat$Group = ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"tumor","normal")
# library(ggpubr)
# 
# ggplot(dat,aes(Cell_type,Proportion,fill = Group)) +
#   geom_boxplot(outlier.shape = 21,color = "black") +
#   theme_bw() +
#   labs(x = "Cell Type", y = "Estimated Proportion") +
#   theme(legend.position = "top") +
#   theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
#   scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
# 
# 
# ggsave(filename = "top10_Est_Celltype(4).png",units = "cm",width = 20,height = 14)




# All
# 结果可视化 
# 热图
load(file = "Cibersort_Allmarker.Rdata")
load(file = "../Step2-0/TCGA_Step1 output.Rdata")
re <- result_Allmarker[,-(29:31)]
library(pheatmap)
k <- apply(re,2,function(x) {sum(x == 0) < nrow(result_Allmarker)/2})
table(k)
re2 <- as.data.frame(t(re[,k]))
an = data.frame(group = group_list,
                row.names = colnames(NSCLCcount))
pdf(file = "../../Fig/Step2-1/Cibersort_pheatmap.pdf",width = 12,height = 8)
pheatmap(re2,scale = "row",
         show_colnames = F,
         # annotation_col = group_list,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

# 直方图
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

library(ggsci)
mypalette = pal_jco()(10)
mypalette1 = pal_ucscgb()(10)
mypalette2 = pal_lancet()(8)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Cibersort Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = c(mypalette1,mypalette,mypalette2))
ggsave(filename = "../../Fig/Step2-1/Cibersort_bar.pdf",width = 16,height = 10)




# 箱线图

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(28))
ggsave(filename = "../../Fig/Step2-1/Cibersort_box.pdf",width = 16,height = 10)



# 让箱线图有顺序
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)


ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(28))

ggsave(filename = "../../Fig/Step2-1/Cibersort_box2.pdf",width = 16,height = 10)



# 根据分组画箱图
library(stringr)
dat$Group = ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"tumor","normal")
library(ggpubr)
library(ggsci)
mypalette = pal_ucscgb()(10)
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Cibersort  Proportion") +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "t.test")


ggsave(filename = "../../Fig/Step2-1/Cibersort_box2.pdf",units = "cm",width = 25,height = 14)


# ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
#   geom_boxplot(outlier.shape = 21,color = "black") + 
#   theme_bw() + 
#   labs(x = "Cell Type", y = "Cibersort  Proportion") +
#   theme(legend.position = "top") + 
#   theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
#   scale_fill_manual(values = mypalette[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
# 
# ggsave(filename = "../../Fig/Step2-1/Cibersort_box3.pdf",units = "cm",width = 20,height = 14)




