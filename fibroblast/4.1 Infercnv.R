## ABCD,HCC1-4,CHC 的成纤维细胞cnv事件推测，以及
## ABCD肝细胞cnv事件推测
## 2021.11.20

library(infercnv)
library(ggplot2)

# 做标签：ABCD用免疫细胞 cluster8 -------------------------------------------------------------------------

load(paste0(HOME,"2.T_N_CAF/ord/sample_anno/A.RData"))

load("/spacial_transcriptome/RESULT_09/1_4sample/1_cell_type/seurat_combined0.5.RData")
lym_id <- rownames(seurat_combined@meta.data)[which(seurat_combined@meta.data$integrated_snn_res.0.5 == 8)] # 979
tmp <- data.frame(lym = lym_id)
tmp$patient <- apply(tmp,1,function(x){strsplit(x[[1]],"_")[[1]][2]})
tmp$spot <- apply(tmp,1,function(x){strsplit(x[[1]],"_")[[1]][1]})

info <- read.csv(file="/T_N_H_S/D_4position.csv")
info <- info[-which(is.na(info$X4_position)),]
info$position <- apply(info,1,function(x){strsplit(x[[2]],"_")[[1]][1]}) 
info$spot <- paste0(info$Barcode,"_4")
info <- info[which(info$position=="N"),3:4]

info_final <- rbind(info_final,info)
save(info_final,file="info_final.RData")

# 做标签：ABCD 肝细胞spot -------------------------------------------------------------------------
load(paste0(HOME,"2.T_N_CAF/sample_anno/A.RData"))
load("/spacial_transcriptome/RESULT_09/1_4sample/1_cell_type/seurat_combined0.5.RData")
hepe_id <- rownames(seurat_combined@meta.data)[which(seurat_combined@meta.data$integrated_snn_res.0.5 %in% c(0,1,2,5,6,9,11))] 

# 本地得到T,N位置
info <- read.csv(file="D_4position.csv")
info <- info[-which(is.na(info$X4_position)),]
info$position <- apply(info,1,function(x){strsplit(x[[2]],"_")[[1]][1]}) 
info$spot <- paste0(info$Barcode,"_4") #改
info <- info[,3:4]

info_final <- rbind(info_final,info)
info_final$patient <- apply(info_final,1,function(x){strsplit(x[[2]],"_")[[1]][2]})
info_final$barcode <- apply(info_final,1,function(x){strsplit(x[[2]],"_")[[1]][1]})
hepe_info <- hepe_info[which(hepe_info$spot %in% hepe_id),]
save(hepe_info,file="hepe_info.RData")

# CAF的cnv-------------------------------------------------------------------------
load("/picb/bigdata/project/CancerSysBio/HeYuFei_Project2/spacial_transcriptome/RESULT_09/1_4sample/1_cell_type/seurat_combined0.5.RData")
load("/picb/bigdata/project/CancerSysBio/HeYuFei_Project2/spacial_transcriptome/RESULT_1026/RESULT/non_CAF_info.RData")

# HCC1-4：用N中非CAF细胞作为reference
for(i in c("HCC1","HCC2","HCC3","HCC4")){
  HOME="/picb/bigdata/project/CancerSysBio/HeYuFei_Project2/spacial_transcriptome/RESULT_1026/RESULT/"
  load(paste0(HOME,"1.seurat_patient/seurat_",i,".RData")) # expr
  
  load(paste0(HOME,"2.T_N_CAF/patient_anno/",i,"_CAF.RData")) #注释
  test <- test[-which(rownames(test)%in% non_CAF_info$spot),]
  
  out_dir=paste0(HOME,'8.infercnv/',i,'/')
  
  #制作标签
  cell_type_annotation <- data.frame(cells = rownames(seurat_combined@meta.data))
  cell_type_annotation$sample <- apply(cell_type_annotation,1,function(x){
    
    strsplit(x[[1]],"_")[[1]][2]})  #列 patient
  cell_type_annotation[which(cell_type_annotation$sample == 1),"type"] = 0  #远端N
  rownames(cell_type_annotation) <- cell_type_annotation$cells
  
  cell_type_annotation[test[which(test$position == "N"),"spot"],"type"] = 1 # N和I_N中CAF
  cell_type_annotation[test[which(test$position == "I"),"spot"],"type"] = 2
  cell_type_annotation[test[which(test$position == "T"),"spot"],"type"] = 3
  cell_type_annotation <- cell_type_annotation[-which(is.na(cell_type_annotation$type)),]
  cell_type_annotation <- cell_type_annotation[3]
  cell_type_annotation$type <- as.factor(cell_type_annotation$type)
  
  # raw_mat
  raw_mat <- seurat_combined$Spatial@counts[,rownames(cell_type_annotation)]
  dim(raw_mat)
  # gene标注
  gene_ordering_table <- read.table(paste0(HOME,'4.infercnv0/gene_ordering_file.txt'),row.names=1)
  
  # 创建对象
  infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=raw_mat, 
                                                 annotations_file=cell_type_annotation, 
                                                 gene_order_file=gene_ordering_table,  
                                                 ref_group_names= '0', 
                                                 chr_exclude=c("chrX","chrY","chrM","KI270734.1")) 
  
  #带过滤的运行
  infercnv_obj2 <- infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir=out_dir,
                                 cluster_by_groups=T,
                                 denoise=T,
                                 scale_data= T, # 参考和推测的细胞类型不是同一种
                                 sd_amplifier=1.5,
                                 HMM=T)
  infercnv_obj2 <- infercnv::apply_median_filtering(infercnv_obj2)
}
