## 将ST-fibroblast ~ sc-fibroblast
## 2022.05.06

library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)

# 1-- sc:得到 marker genes -------------------------------------------------------------------------
sc <- readRDS("fibroblast.rds") 

Idents(sc) <- "RNA_snn_res.0.2"
sc <- subset(sc,idents = c(0,1,2,3,4,5))) # 去除6，7
sc <- SCTransform(sc,assay = "RNA",verbose=FALSE,
                  vars.to.regress=c("percent.mito","nFeatures_RNA","PatientID"),
                  variable.features.n = 2000) 

DefaultAssay(sc) <- "SCT"

saveRDS(sc,file="fibroblast.rds")

all.markers <- FindAllMarkers(sc,only.pos = T)
save(all.markers,file="fibroblast_marker.RData")

markers <- all.markers %>% 
  filter(p_val_adj < 0.05)
summary(as.factor(markers$cluster))


# pathway -------------------------------------------------------------------------

sc <- readRDS("spacial_transcriptome/RESULT_1026/2.sc-Fs/0510/fibroblast.rds")
Idents(sc) <- "RNA_snn_res.0.2"
# 平均表达
expr <- AverageExpression(sc,assays="SCT",slot="data")[[1]]
expr <- as.matrix( expr[rowSums(expr)>0,] )
dim(expr) #13852,6

# gsva
gsva_GO <- gsva(expr, gsGO, min.sz = 10)
gsva_KEGG <- gsva(expr, gsKEGG, min.sz = 10)
save(gsva_GO,gsva_KEGG)

# 可视化
plot_df <- 
  
peatmap(plot_df)

# 2-- ST: 计算score -------------------------------------------------------------------------
setwd("/csb2/project/HeYuFei_TME_2021/spacial_transcriptome/RESULT_1026/RESULT/9.marker_gene/1202/")
load("seurat_marker.RData") # ST CAF
load("sc_fibroblast_marker.RData")  # sc all.markers

markers <- all.markers %>% 
  filter(p_val_adj < 0.05)

# zscore
seurat_marker <- subset(seurat_marker,idents=c(0,1,2,3,4,5))
my_table <- seurat_marker@meta.data
expr <- seurat_marker$SCT@data

for (i in 1:6){
  genes <- intersect(rownames(expr), markers %>% filter(cluster == i) %>% pull(gene))
  counts <- length(genes)
  #my_table[,paste0("sc_",i)] <- colSums(expr[genes,])/counts
  my_table[,paste0("sc_zscore",i)] <- scale(colSums(expr[genes,]))
  
}

# Tsne 可视化
seurat_marker@meta.data <- my_table
FeaturePlot(seurat_marker,features="sc_zscore3",pt.size=1,order=T,cols=c("lightgray","red"),
            min.cutoff = "q70",max.cutoff="q99")
FeaturePlot(seurat_marker,features="sc_zscore5",pt.size=1,order=T,cols=c("lightgray","red"),
            min.cutoff = "q60",max.cutoff="q99")


# 3-- 画图：气泡图
plot_df <- my_table[,c(13:19)]

plot_df2  <- plot_df %>% 
  group_by(seurat_clusters) %>%
  summarise(sc1 = median(sc_zscore1),
            sc2 = median(sc_zscore2),
            sc3 = median(sc_zscore3),
            sc4 = median(sc_zscore4),
            sc5 = median(sc_zscore5),
            sc6 = median(sc_zscore6)) 

plot_df3  <- plot_df %>% 
  group_by(seurat_clusters) %>%
  summarise(sc1 = mean(sc_zscore1),
            sc2 = mean(sc_zscore2),
            sc3 = mean(sc_zscore3),
            sc4 = mean(sc_zscore4),
            sc5 = mean(sc_zscore5),
            sc6 = mean(sc_zscore6)) 

# 4-- 查看 gene overlap -------------------------------------------------------------------------

# F5 marker gene
ST_all.markers$cluster <- as.factor(as.numeric(ST_all.markers$cluster)) #1-7
gene_ST <- ST_all.markers %>% filter(cluster == 5,p_val_adj < 0.05) %>% pull(gene)

all.markers$cluster <- as.factor(as.numeric(all.markers$cluster)) #1-6
# sc4 marker gene
gene_sc4 <- all.markers %>% filter(cluster == 4,p_val_adj < 0.05) %>% pull(gene)

# sc5 marker gene
gene_sc5 <- all.markers %>% filter(cluster == 5,p_val_adj < 0.05) %>% pull(gene)

intersect(gene_ST,gene_sc4) #2 
intersect(gene_ST,gene_sc5) #49
intersect(gene_sc4,gene_sc5) #0

# AUCell-------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(stringr)

setwd("/csb2/project/HeYuFei_TME_2021/spacial_transcriptome/RESULT_1026/RESULT/9.marker_gene/1202/")
load("seurat_marker.RData")

# ST的pancerCAF 打分
AUCell <- read.table("AUCell/1.AUCell_Pancancer_marker.txt") %>% t() %>% as.data.frame() %>% scale()
# ST的sc1-6打分
AUCell <- read.table("AUCell/2.AUCell_sc_fibType.txt") %>% t() %>% as.data.frame() %>% scale()

rownames(AUCell) <- str_replace(rownames(AUCell),"\\.","-")
AUCell_table <- merge(seurat_marker@meta.data[,c("position","seurat_clusters")],
                      AUCell, by=0) %>% as.data.frame()