## 所有成纤维细胞群体的分析
## 2021.12.4

library(Seurat)
#1. 整合并提取所有的成纤维细胞 -------------------------------------------------------------------------

seurat.list <- list()
sample_id <- 1
for (folder_tmp in folder){  # 不需要图像了，按照单细胞分析
  img <- Read10X_Image(paste0(folder_tmp,'spatial'))
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  data=Read10X_h5(paste0(folder_tmp,'filtered_feature_bc_matrix.h5')) 
  seurat_single=Seurat::CreateSeuratObject(data,assay='Spatial')
  barcode = info[which(info$sample == sample_id),"Barcode"]
  seurat_single=subset(seurat_single,cells=barcode) # 取成纤维细胞
  seurat_single[["percent.mt"]] <- PercentageFeatureSet(seurat_single, pattern = "^MT-")
  seurat_single@meta.data$sample_id <- sample_id
  seurat_single=SCTransform(seurat_single,assay = "Spatial",verbose=FALSE,
                            vars.to.regress=c("percent.mt","nFeature_Spatial"),
                            variable.features.n = 3500)
  
  seurat.list[[sample_id]] <- seurat_single
  sample_id <- sample_id+1
}

# merge 10号
merge1 <- merge(seurat.list[[10]], y=seurat.list[[13]])
merge1 =SCTransform(merge1,assay = "Spatial",verbose=FALSE,
                    vars.to.regress=c("percent.mt","nFeature_Spatial"),
                    variable.features.n = 3500)            

seurat.list <- seurat.list[c(-10,-13)]
seurat.list[24] <- merge1 # 将merge 对象放到后面

# 整合开始
features <- SelectIntegrationFeatures(seurat.list)
seurat.list <- PrepSCTIntegration(object.list=seurat.list, anchor.features=features)
seurat_anchors <- FindIntegrationAnchors(object.list = seurat.list,anchor.features = features,
                                         dims = 1:20, reduction="cca",
                                         k.anchor=5,   k.filter=20,  k.score=20)  #改
feature_integrate <- unique(rapply(
  lapply(seurat_anchors@object.list,function(x){x@assays[["SCT"]]@data@Dimnames[[1]]} ),
  function(y){return(y)},how='unlist'))

seurat_combined <- IntegrateData(anchorset = seurat_anchors, features.to.integrate=feature_integrate,
                                 dims = 1:20,
                                 k.weight=10)  # 容易错，改
save(seurat_combined_ord.RData)

#2. 聚类 -------------------------------------------------------------------------

DefaultAssay(seurat_combined) <- "integrated"
seurat_combined <- ScaleData(seurat_combined, features=feature_integrate,verbose = FALSE)
features <- SelectIntegrationFeatures(seurat.list, nfeatures=3000)  #对feature排序,得到3000个基因
seurat_combined <- RunPCA(seurat_combined, features=features, npcs = 30, verbose = FALSE)
seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:30)
seurat_combined <- FindNeighbors(seurat_combined, reduction = "pca", dims = 1:30)

# 更改metadata,匹配position信息
load(info)
seurat_combined@meta.data$spot <- rownames(seurat_combined@meta.data)
seurat_combined@meta.data$X <- apply(seurat_combined@meta.data,1,function(x){
  strsplit(x[["spot"]],"_")[[1]][1]})
seurat_combined@meta.data$info_spot <- paste(seurat_combined@meta.data$X,
                                             seurat_combined@meta.data$sample_id, sep="_")
tmp <- merge(seurat_combined@meta.data,info[,c("spot","position")],by.x="info_spot",by.y="spot")
rownames(tmp) <- tmp$spot
seurat_combined@meta.data <- tmp[rownames(seurat_combined@meta.data),] #得保证顺序一样

# 亚型
seurat_combined0.5 <- FindClusters(seurat_combined, resolution = 0.5)
save(seurat_combined0.5.RData)

#3. 得到表达数据 -------------------------------------------------------------------------

seurat_marker <- SCTransform(seurat_combined0.5,assay = "Spatial",verbose=FALSE,
                             vars.to.regress=c("percent.mt","nFeature_Spatial"),
                             variable.features.n = 3500) 
save(seurat_marker.R)