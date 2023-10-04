## irGSEA计算 module score
# env：

library(irGSEA)
library(Seurat)
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)

setwd("13.module_score_irGSEA/")

# 准备数据 -------------------------------------------------------------------------
# 加载基因集
geneset_cancer <- read_xlsx("module_gene.xlsx",sheet = 1)
geneset_immune <- read_xlsx("module_gene.xlsx",sheet = 2)
geneset_stem <- read_xlsx("module_gene.xlsx",sheet = 3)

# 加载Seurat 数据
my_table <- readRDS("1.seurat_combined_metadata2.rds" )
load( "seurat_combined0.5.RData" )

# 将分组信息加入 metadata
load("immuneBarcode_3type.RData")
load("cancerBarcode_3type.RData")

my_table %<>% mutate(cancer_3type = ifelse(spot %in% fib_neHepa,"fib_neHepa", # bg 和 fib有重合，还没有改
                                    ifelse(spot %in% CAF_neHepa,"CAF_neHepa",
                                    ifelse(spot %in% bg_Hepa,"bg_Hepa","non_cancer"))),
                     immune_3type = ifelse(spot %in% fib_neStro,"fib_neStro",
                                    ifelse(spot %in% CAF_neStro,"CAF_neStro",
                                    ifelse(spot %in% bg_Stro,"bg_Stro","non_stroma")))
                       )

seurat_combined@meta.data <- my_table

# irGSEA: immune -------------------------------------------------------------------------
Idents(seurat_combined) <- "immune_3type"
gs <- lapply(geneset_immune, na.omit) # 转为List 并删NA
seurat_obj <- seurat_combined

seurat_obj <- irGSEA.score(object = seurat_obj, 
                           assay = "SCT", slot = "data", 
                           new.assay='ssgsea',
                           geneset = gs, 
                           seeds = 123, ncores = 10,
                           method = c("ssgsea"),  #"AUCell", "UCell", "ssgsea","singscore"
                           # aucell.MaxRank = 3457,
                           aucell.MaxRank = ceiling(0.2*nrow(seurat_obj)), 
                           ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian', 
                           custom = T)
Seurat::Assays(seurat_obj)

AUCell <- seurat_obj@assays$AUCell@data
UCell <- seurat_obj@assays$UCell@data
singscore <- seurat_obj@assays$singscore@data
ssgsea <- seurat_obj@assays$ssgsea@data

# 计算差异gene set
result.dge <- irGSEA.integrate(object = seurat_obj,  # 不知道为啥取subset后会报错
                               group.by = "immune_3type",
                               metadata = NULL, col.name = NULL,
                               method = c("UCell","singscore",
                                          "ssgsea"))

save(UCell,ssgsea,singscore,result.dge, file="1.ModuleScore_immune.RData")


# irGSEA: cancer -------------------------------------------------------------------------
Idents(seurat_combined) <- "cancer_3type"
gs <- lapply(geneset_cancer, na.omit) # 转为List 并删NA
seurat_obj <- seurat_combined

seurat_obj <- irGSEA.score(object = seurat_obj, 
                           assay = "SCT", slot = "data", 
                           geneset = gs, 
                           seeds = 123, ncores = 10,
                           method = c("UCell", "singscore", "ssgsea"), 
                           ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian', 
                           custom = T)
Seurat::Assays(seurat_obj)

UCell <- seurat_obj@assays$UCell@data
singscore <- seurat_obj@assays$singscore@data
ssgsea <- seurat_obj@assays$ssgsea@data

# 计算差异gene set
result.dge <- irGSEA.integrate(object = seurat_obj,  # 不知道为啥取subset后会报错
                               group.by = "cancer_3type",
                               metadata = NULL, col.name = NULL,
                               method = c("UCell","singscore",
                                          "ssgsea"))

save(UCell,ssgsea,singscore,result.dge, file="1.ModuleScore_cancer.RData")

# 画图-------------------------------------------------------------------------

irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)

ggsave(irGSEA.heatmap.plot, file="irGSEA.heatmap.plot.pdf")