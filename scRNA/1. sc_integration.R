## scRNA成纤维细胞整合
## jingsiyu-20220325

library(dplyr)
library(Seurat)
library(harmony)


# all cell type-------------------------------------------------------------------------
HCC <- readRDS("GSE156625_HCCF_integrated.RData")
# 我需要的数据
barcode <- HCC@meta.data %>% 
  filter(#condition %in%c("Adj Normal","Tumor"),
    #PatientID != "HN",
    nFeature_RNA>200) %>%
  row.names(.)

HCC <- subset(HCC,cells=barcode)
summary(as.factor(HCC@meta.data$CellType))
DimPlot(HCC,group.by="CellType")


# fibroblast整合：harmony ----------------------------------------------------


HCC<-readRDS("GSE156625_HCCF_integrated.RData") #seuratobj from Rdata
HCC.list <- SplitObject(HCC, split.by = "PatientID")
HCC.list$Fw16=NULL
HCC.list$Fw21=NULL
HCC.list$Fw18=NULL
HCC.list$Fw14=NULL
HCC.list$HN=NULL

#fibroblast
fibroblast=subset(x = HCC, subset=CellType == "Fibroblast")
pateint.fb.list <- SplitObject(fibroblast, split.by = "PatientID")
pateint.fb.list$Fw14=NULL
pateint.fb.list$Fw16=NULL
pateint.fb.list$Fw18=NULL
pateint.fb.list$Fw21=NULL
pateint.fb.list$HN=NULL
#merge
human.list<-pateint.fb.list[c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13")]
#HCC.list<-lapply(X=HCC.list,FUN=FindVariableFeatures)
mymerg<-merge(pateint.fb.list$P14,human.list)
#ln(counts)转counts
DefaultAssay(mymerg) <- "RNA"
#raw_integrated
mymerg <- NormalizeData(mymerg) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
mymerg <- RunHarmony(mymerg, group.by.vars = "PatientID")
mymerg <- RunUMAP(mymerg, reduction = "harmony", dims = 1:30)
mymerg <- FindNeighbors(mymerg, reduction = "harmony", dims = 1:30) 
mymerg <- FindClusters(mymerg,resolution =c(0.2,0.4,0.5,0.6,0.8,1))

saveRDS(mymerg,"fibro_harmony.rds")

# 最后分亚型了之后有一些是内皮
Idents(fibroblast) <- 'new_group'
subset(fibroblast,idents=c('N_3','T_2'),invert=T)



# 画图 -------------------------------------------------------------------------
# scRNA t-sne
HCC <- readRDS("GSE156625_HCCF_integrated.RData") 
Idents(HCC) <- "CellType"
DimPlot(HCC)

# sc 亚型差异基因 tsne-------------------------------------------------------------------------
sc <- readRDS("fibroblast.rds")

p<- FeaturePlot(sc,features="COL1A2",pt.size=1,order=T,cols=c("lightgray","red"),
            min.cutoff = "q60",max.cutoff="q98")