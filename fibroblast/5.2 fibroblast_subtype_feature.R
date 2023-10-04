## 成纤维细胞群体后续分析
## 12.4

library(Seurat)
library(dplyr)

#1. 富集在T中的cluster: fisher.test-------------------------------------------------------------------------
load("seurat_combined0.5")
position_df <- metadata %>% 
  group_by(integrated_snn_res.0.5,position) %>%
  summarise(n=n())

colnames(position_df) <- c("res0.5","position","n")

# T&!T  (同理 T&N,TI&N)
for (i in 1:7){
  print(i)
  position_df %>% filter(res0.5 == i,position != "N") %>% pull(n) %>% sum() -> a
  position_df %>% filter(res0.5 == i,position == "N") %>% pull(n) %>% sum() -> c
  position_df %>% filter(res0.5 != i,position != "N") %>% pull(n) %>% sum() -> b
  position_df %>% filter(res0.5 != i,position == "N") %>% pull(n) %>% sum() -> d
  
  x <- matrix(c(a,c,b,d),ncol=2,nrow=2)
  p <- fisher.test(x,alternative = "greater")$p.value
  print(p)
}

#2. cluster marker gene -------------------------------------------------------------------------

load(seurat_marker)
all.markers <- FindAllMarkers(seurat_marker, only.pos=T, assay="SCT",slot = "data")
save(all.markers)

all.markers <- all.markers %>% filter(p_val_adj < 0.05)

# 筛marker gene
# i <- 0
genes <- c()
for (i in c(0,1,3:6)){
  tmp1 <- all.markers %>% filter(cluster==i)
  if (nrow(tmp1) >10){
    tmp <- tmp1 %>% 
      filter(avg_log2FC > 0.7) %>%
      pull(gene) 
    print(length(tmp))
    genes <- unique(c(genes,tmp))
  }else{
    tmp <- tmp1 %>%
      pull(gene) 
    print(length(tmp))
    genes <- unique(c(genes,tmp))
  }
}


# cluster 5 通路富集-------------------------------------------------------------------------

library(GSVA)
library(Seurat)
library(dplyr)

load("spacial_transcriptome/RESULT_1026/RESULT/gene_set/geneset_GOkegg.RData")
load("spacial_transcriptome/RESULT_1026/RESULT/gene_set/gene_base_hallmarkers.RData")   
load(all.markers)
expr <- readRDS("C:/Users/DELL/Desktop/expr_SCT4.rds")

genes <- all.markers %>% 
  filter(p_val_adj < 0.05, cluster==4) %>%
  pull(gene)

# 所有亚群功能富集-------------------------------------------------------------------------
library(GSVA)
library(Seurat)
library(dplyr)

load("geneset_GOkegg.RData")

# F1-6
load(seurat_marker)
Idents(seurat_marker) <- "seurat_clusters"
# 平均表达
expr <- AverageExpression(seurat_marker,assays="SCT",slot="data")[[1]]
expr <- as.matrix( expr[rowSums(expr)>0,] )
dim(expr) #19892,7

# gsva
gsva_GO <- gsva(expr, gsGO, min.sz = 10)
gsva_KEGG <- gsva(expr, gsKEGG, min.sz = 10)
save(gsva_GO,gsva_KEGG,file="subset_pathway.RData")
