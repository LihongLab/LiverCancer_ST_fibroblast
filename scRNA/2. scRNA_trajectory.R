## scRNA 伪时序分析
## jingsiyu-2022-0506

library(dplyr)
library(Seurat)
library(monocle)

fibroblast <- readRDS("fibroblast.rds")
# 不考虑T5
fibroblast <- subset(fibroblast,idents=c(6,7),invert=T) #剩下1836个

## 构建cds
expr_matrix <- as(as.matrix(fibroblast@assays$RNA@counts), 'sparseMatrix')
p_data <- fibroblast@meta.data 
f_data <- data.frame(gene_short_name = row.names(fibroblast),
                     row.names = row.names(fibroblast))

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

## 评估
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## 过滤gene
cds <- detectGenes(cds, min_expr = 0.1) #得到num_cells_expressed
print(head(fData(cds))) #此时有13714个基因
expressed_genes <- subset(fData(cds),num_cells_expressed >= 10) %>% row.names(.)


# seurat 差异gene-------------------------------------------------------------------------
Idents(fibroblast) <- "RNA_snn_res.0.2"
deg.cluster <- FindAllMarkers(fibroblast)
deg.cluster <- subset(deg.cluster,p_val_adj<0.05)
write.table(deg.cluster,file="monocle_trainDEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

express_genes <- deg.cluster %>% pull(gene)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)


# 降维和排序-------------------------------------------------------------------------

cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)

## 可视化
plot_cell_trajectory(cds,color_by="RNA_snn_res.0.2", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)


library(ggsci)  #自定义颜色
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
plot_cell_trajectory(cds, color_by = "new_group2")  + scale_color_manual(values = colour)

saveRDS(cds,"monocle_cds.rds")



# 空转伪时序分析 -------------------------------------------------------------------------

library(dplyr)
library(Seurat)
library(monocle)

load('/picb/bigdata/project/CancerSysBio/Project/HeYuFei_Project2/spacial_transcriptome/RESULT_1026/RESULT/9.marker_gene/1202/seurat_combined0.5.RData')
fibroblast <- seurat_combined0.5
fibroblast <- subset(fibroblast,idents=6,invert=T) #不考虑F6

## 构建cds

expr_matrix <- as.matrix(fibroblast@assays$Spatial@counts)
p_data <- fibroblast@meta.data 
f_data <- data.frame(gene_short_name = row.names(expr_matrix),
                     row.names = row.names(expr_matrix))

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

## 评估
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## 过滤gene
cds <- detectGenes(cds, min_expr = 0.1) #得到num_cells_expressed
print(head(fData(cds))) #此时有13714个基因
expressed_genes <- subset(fData(cds),num_cells_expressed >= 10) %>% row.names(.)


# seurat 差异gene
Idents(fibroblast) <- "seurat_clusters"
deg.cluster <- FindAllMarkers(fibroblast)
deg.cluster <- subset(deg.cluster,p_val_adj<0.05)
write.table(deg.cluster,file="/picb/bigdata/project/CancerSysBio/Project/HeYuFei_Project2/spacial_transcriptome/RESULT_1026/RESULT/9.marker_gene/1202/monocle_trainDEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

express_genes <- deg.cluster %>% pull(gene)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)


# 降维和排序
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)

## 可视化
plot_cell_trajectory(cds, size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)


library(ggsci)  #自定义颜色
cluster_col <- c("#E9C46A","#F4A261","#a8dadc","#e76f51",
                 "#457b9d","#264653","#2A9D8F")

plot_cell_trajectory(cds)  + scale_color_manual(values = cluster_col)

saveRDS(cds,"monocle_cds.rds")

