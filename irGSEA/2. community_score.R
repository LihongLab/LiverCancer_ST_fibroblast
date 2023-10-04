## irGSEA 输出的score matrix ,计算 community score
## 2022.10.24
## jingsiyu

library(Seurat)
library(dplyr)
library(magrittr)

setwd("/csb2/project/HeYuFei_TME_2021/spacial_transcriptome/RESULT_1026/RESULT/13.module_score_irGSEA/")

irgsea_methods <- c("singscore","ssgsea","UCell")
my_table <- readRDS( "RESULT_1026/1.cell_type/1.seurat_combined_metadata2.rds" )
nei_info <- readRDS("nei_info_4sample.rds")


# 计算 community score: cancer -------------------------------------------------------------------------

load("1.ModuleScore_cancer.RData") # irGSEA的结果环境, pathway_num * 17381

for (irmethod in irgsea_methods){
  assay <- get(irmethod) %>% as.matrix() #循环三种方法的结果
  
  # 以cancer为例,stem 同
  plot_df <- my_table %>%  # all fib&TSF 的 TSF score df, 722
    filter(TSF_type %in% c("Fibroblast","TSF_high")) %>%  # cell_type=="Fibroblast"
    select(spot,TSF_score)
  community_df <- nei_info %>%  # neighbor is cancer cell. 959 pairs
    filter(First_spot %in% (plot_df %>% pull(spot)),
           Second_spot %in% (my_table %>% filter(TSF_type=="cancer Hepatocyte") %>% pull(spot)))
  
  Second_moduleScore <- assay[, community_df[["Second_spot"]]] %>% t()
  tmp <- rownames(Second_moduleScore)
  Second_moduleScore %<>% as.data.frame() %>% mutate(Second_spot = tmp)
  
  community_df <- left_join(community_df,Second_moduleScore,by=c("Second_spot"="Second_spot"))
  
  niche_score_sum <- community_df %>% # calculate community score: sum or max
    group_by(First_spot) %>% 
    summarise_at(vars(contains(rownames(assay))), sum, na.rm=TRUE)  # 356*20
  
  plot_df <- left_join(plot_df,niche_score_sum,by=c("spot"="First_spot")) # 722*20
  
  assign( paste0(irmethod,"_cancer_sum"), plot_df) # 保存为变量
}

save(list=ls()[str_detect(ls(),"cancer_")], # 变量
     file="2.communityScore_cancer.RData" )