## 10.26 nichenet 运行
## jingsiyu
## 2022.20.26


library(Seurat)
library(nichenetr)

setwd("/csb2/project/HeYuFei_TME_2021/spacial_transcriptome/RESULT_1026/RESULT/14.nichenet/")

# 1--run: 加载先验数据库 -------------------------------------------------------------------------

ligand_target_matrix <- readRDS("ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]

lr_network <- readRDS("lr_network.rds")
head(lr_network)

weighted_networks <- readRDS("weighted_networks.rds") # list
weighted_networks_lr <- weighted_networks$lr_sig %>% 
  inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
head(weighted_networks$lr_sig)

# 加载 seurat -------------------------------------------------------------------------
my_table <- readRDS("my_table.rds")
load( paste0(DATA_DIR,"RESULT_09/1_4sample/1_cell_type/seurat_combined0.5.RData") )

seurat_combined@meta.data <- my_table

Idents(seurat_combined) <- "nichenet_2type"
seurat_obj <- subset(seurat_combined, idents=c( "cancer","Fibroblast" ))
summary(as.factor(seurat_obj@meta.data$condition))
summary(as.factor(seurat_obj@meta.data$nichenet_2type))

# run: nichenet-------------------------------------------------------------------------
# fibroblast -> cancer
nichenet_output <- nichenet_seuratobj_aggregate(seurat_obj = seurat_obj, 
                                                top_n_ligands = 40,
                                                receiver = "cancer",
                                                sender = "Fibroblast",
                                                condition_colname = "condition", 
                                                condition_oi = "stemness", 
                                                condition_reference = "general", 
                                                assay_oi= "SCT", #不然会用Integrated
                                                ligand_target_matrix = ligand_target_matrix, 
                                                lr_network = lr_network, 
                                                weighted_networks = weighted_networks)

saveRDS(nichenet_output,"nichenet_output_ST_f2c.rds")
