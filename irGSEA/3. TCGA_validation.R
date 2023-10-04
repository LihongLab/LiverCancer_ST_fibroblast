##TCGA 验证
## 20221024
## jingsiyu

library(irGSEA)
library(Seurat)
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)

setwd("/csb2/project/HeYuFei_TME_2021/spacial_transcriptome/RESULT_1026/RESULT/13.module_score_irGSEA/")

# 准备数据 -------------------------------------------------------------------------
# 加载基因集
geneset_cancer <- read_xlsx("module_gene.xlsx",sheet = 1)[["stem"]] %>% na.omit()
geneset_immune <- read_xlsx("module_gene.xlsx",sheet = 2)[["M2"]] %>% na.omit()
gs <- list(stem = geneset_cancer,
           M2 = geneset_immune,
           TSF = c("COL4A1","COL1A2","FSTL1","CTGF","COL4A2"))

# 加载 TCGA 数据
TCGA_rna <- read.table("mixture_file_TCGA_LIHC_Primary_FPKM.txt")

# TCGA 病人对 TSF 和 stemness, M2 打分-------------------------------------------------------------------------
TCGA_score <- irGSEA.score(object = TCGA_rna, 
                           geneset = gs, 
                           custom = T,
                           seeds = 123, ncores = 10,
                           method = c("UCell", "ssgsea","singscore"),  #"AUCell", "UCell", "ssgsea","singscore"
                           ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian', 
                           min.cells = 3, 
                           min.feature = 0)

Seurat::Assays(TCGA_score)
UCell <- TCGA_score@assays$UCell@data
singscore <- TCGA_score@assays$singscore@data
ssgsea <- TCGA_score@assays$ssgsea@data

save(UCell,ssgsea,singscore, file="3.ModuleScore_TCGA.RData")

# TCGA: 相关性 -------------------------------------------------------------------------
# env sc10
library(dplyr)
library(ggplot2)
library(Seurat)
library(ggpubr)
setwd("/csb2/project/HeYuFei_TME_2021/spacial_transcriptome/RESULT_1026/RESULT/13.module_score_irGSEA/")

load("3.ModuleScore_TCGA.RData")
irgsea_methods <- c("singscore","ssgsea","UCell")

for (irmethod in irgsea_methods){
  assay <- get(irmethod) %>% as.matrix() %>% t() %>% as.data.frame() #循环三种方法的结果
  module_names <- colnames(assay)
  
  for (module in c("stem","M2")){
    p <- ggplot(assay, aes(x=TSF, y=assay[[module]])) + 
      geom_point()+ geom_smooth(method = 'lm', se = T, color = 'red') +
      labs(y=module, title=irmethod ) + 
      stat_cor(method = "pearson") + 
      theme_classic()
    ggsave(p, file=paste0("TCGA_plot/",irmethod,"_",module,".pdf"))
  }
}