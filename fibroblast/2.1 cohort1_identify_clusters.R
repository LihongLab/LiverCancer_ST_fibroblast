## 4样本原始数据整合，聚类,确定marker gene与细胞类型
## jingsiyu-2021.10


library(Seurat)
library(magrittr)
library(dplyr)
library(patchwork)
library(ggplot2)
# 重新合并文件
setwd("/spacial_transcriptome/RAW_DATA/")
folder <- c('1.Count/A/',
            '1.Count/B/',
            '1.Count/C/',
            '1.Count/D/')



# 1 --得到seurat.list --------------------------------------------------------------------
seurat.list <- sapply(folder,function(folder_tmp){
        img <- Read10X_Image(paste0(folder_tmp,'spatial'))
        Seurat::DefaultAssay(object = img) <- 'Spatial'
        
        data=Read10X_h5(paste0(folder_tmp,'filtered_feature_bc_matrix.h5')) 
        seurat_single=Seurat::CreateSeuratObject(data,assay='Spatial') #没有QC!!
        seurat_single[['image']] <- img
        img <- img[colnames(x=seurat_single)]
        seurat_single[["percent.mt"]] <- PercentageFeatureSet(seurat_single, pattern = "^MT-")
        seurat_single=SCTransform(seurat_single,assay = "Spatial",verbose=FALSE,
                                  vars.to.regress=c("percent.mt","nFeature_Spatial"),
                                  variable.features.n = 3500)
        return(seurat_single)
})

save(seurat.list,file = "/spacial_transcriptome/JSY/4_21_result/seurat.list.RData")

# 2 --合并way1:pca ----------------------------------------------------------------

seurat.list <- lapply(seurat.list, FUN = function(x) {
        x <- RunPCA(x, verbose = FALSE) })

seurat_anchors <- FindIntegrationAnchors(object.list = seurat.list, reduction = "rpca")

# select all genes for integration
feature_integrate <- unique(rapply(lapply(seurat_anchors@object.list,function(x){
        x@assays[["SCT"]]@data@Dimnames[[1]]
}),function(y){return(y)},how='unlist'))

seurat_combined <- IntegrateData(anchorset = seurat_anchors, dims = 1:20, 
                                 normalization.method = "SCT",
                                 features.to.integrate=feature_integrate)

# 2 --合并way2 少数样本不pca,即CCA ----------------------------------------------------------------

save(seurat.list,file = "./seurat.list.RData")

features <- SelectIntegrationFeatures(seurat.list,nfeatures=3000)  #对feature排序,得到3000个基因
seurat.list <- PrepSCTIntegration(object.list=seurat.list, anchor.features=features) # normalization等准备工作
seurat_anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",
                                         anchor.features = features,
                                         k.filter = 100, k.score = 30) # find anchors

# select all genes for integration
feature_integrate <- unique(rapply(
        lapply(seurat_anchors@object.list,function(x){x@assays[["SCT"]]@data@Dimnames[[1]]} ),
        function(y){return(y)},how='unlist'))

seurat_combined <- IntegrateData(anchorset = seurat_anchors,dims = 1:20, 
                                 normalization.method = "SCT",
                                 features.to.integrate=feature_integrate,
                                 k.weight=50)

# 3 --clustering --------------------------------------------------------------

DefaultAssay(seurat_combined) <- "integrated"

seurat_combined <- ScaleData(seurat_combined, features=feature_integrate,verbose = FALSE)
features <- SelectIntegrationFeatures(seurat.list,nfeatures=3000)  #对feature排序,得到3000个基因
seurat_combined <- RunPCA(seurat_combined, features=features, npcs = 30, verbose = FALSE)
seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:30)
seurat_combined <- FindNeighbors(seurat_combined, reduction = "pca", dims = 1:30)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)

save(seurat_combined, file = "seurat_combined0.5.RData")

# 4 -- find markers ----------------------------------------------------------------------

all.markers <- FindAllMarkers(object = seurat_combined,only.pos=T)
markers <- FindMarkers(seurat_combined, subset.ident = "2") #如果只查看部分cluster

my_table <- data.frame()
id <- levels(all.markers$cluster)
for (i in id){
        idex <- which(all.markers$cluster == i)
        table <- all.markers[idex,]
        table <- table[order(table[,"p_val_adj"])[1:100],]
        
        my_table <- rbind(my_table,table)
}
write.table(my_table, file="./table.txt")