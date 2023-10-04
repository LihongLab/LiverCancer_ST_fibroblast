library(SPOTlight)
library(Seurat)
library(dplyr)

# 单细胞数据处理 -------------------------------------------------------------------------
sc <- readRDS("/singlecell/sc_1500_sample100.rds")
summary(as.factor(sc@meta.data$CellType))
colnames(sc@meta.data)
DefaultAssay(sc) <- "integrated"

# 计算marker gene
cluster_markers_all <- FindAllMarkers(object = sc,assay = "integrated",slot = "data",
                                      verbose = TRUE, only.pos = TRUE, logfc.threshold = 1,
                                      min.pct = 0.7)
cluster_markers_all$gene <- rownames(cluster_markers_all)

# spotlight -------------------------------------------------------------------------

# 加载空转数据
load("F:/Ajingsiyu/project/HeYuFei-TME/DATA/数据-空转/singlecell/seurat_combined.rds") #all spots

seurat_CAF <- SCTransform(seurat_CAF, assay = "Spatial", verbose = FALSE)

# 反卷积
spotlight_ls <- spotlight_deconvolution(se_sc = sc,
                                        counts_spatial = seurat_CAF@assays$Spatial@counts,
                                        clust_vr = "CellType",
                                        cluster_markers = cluster_markers_all,
                                        cl_n = 50,
                                        hvg = 3000,
                                        ntop = NULL,
                                        transf = "uv",
                                        method = "nsNMF",
                                        min_cont = 0.09)
# 结果
decon_mtrx <- spotlight_ls[[2]]  # 细胞占比矩阵
# 行名和空转的 metadata 一致
rownames(decon_mtrx) <- rownames(seurat_CAF@meta.data)
saveRDS(decon_mtrx,file="F:/Ajingsiyu/PROJECT/HeYuFei-TME/RESULT/figure/fig 2-2.cell_type/decon_mtrx.rds")

# 后续分析数据生成 -------------------------------------------------------------------------
library(dplyr)
library(ggplot2)

load("info.RData")
deconv <- readRDS("/decon_mtrx_1500_sample100.rds")
deconv <- as.data.frame(deconv)
summary(deconv$Fibroblast)
colnames(deconv)[7] <- "others"  # 残差就是其他未知的

# rank
df_rank <- apply(deconv[,1:7],1, rank )
df_rank <- as.data.frame(t(df_rank))
df_rank <- floor(8-df_rank)
summary(as.factor(df_rank$Fibroblast)) # fib在6种细胞中的排名

deconv[rownames(df_rank),"fi_rank"] <- df_rank$Fibroblast
deconv[rownames(df_rank),"bin"] <- df_rank$bin

# 位置注释
info <- rbind(CAF21[,c("spot","position","patient")],CAF4[,c("spot","position","patient")])
df_rank$spot <- rownames(df_rank)
df_rank <- merge(df_rank,info,by="spot")
rownames(df_rank) <- df_rank$spot

save(df_rank,file="C:/Users/DELL/Desktop/df_rank.RData")
save(deconv,file="C:/Users/DELL/Desktop/deconv.RData")

# 统计 -------------------------------------------------------------------------
# 统计绝对百分比
fi_percent <- round(x=deconv$Fibroblast,digits=1) #离散化
df_rank[rownames(deconv),"bin"] <- fi_percent  # fib 占比的Bin

df_rank %>%
  group_by(Fibroblast) %>%
  count(bin) -> test

ggplot() + geom_bar(data =test, aes(x=Fibroblast,y=n,fill=as.factor(bin)),stat = "identity",position = "fill")+
  labs(title="percentage of infer fibroblast with \n different ranks in fibroblast spot")+
  theme(plot.title = element_text(size=15),axis.text = element_text(size=14))

# 大体成纤维细胞的占比情况：分为1-6
library(scatterpie)
deconv %>%
  group_by(fi_rank) %>%
  summarise(Bcells = mean(Bcells), Endothelial=mean(Endothelial),Fibroblast=mean(Fibroblast),
            Hepatocytes=mean(Hepatocytes), Immune=mean(Immune),Myeloid=mean(Myeloid),
            others=mean(others)) -> df

df$x <- rep(c(2,4,6),times=2)  # 排图的坐标
df$y <- rep(c(2,4),each=3)
ggplot() + geom_scatterpie(data=df,
                aes(x,y,group=fi_rank,r=0.9),
                cols = c("Bcells", "Endothelial", "Fibroblast", "Hepatocytes", "Immune", "Myeloid", "others"))+
  labs(title="percentage of infer fibroblast with different ranks in fibroblast spot")

# 分为bin
deconv %>%
  group_by(bin) %>%
  summarise(Bcells = mean(Bcells), Endothelial=mean(Endothelial),Fibroblast=mean(Fibroblast),
            Hepatocytes=mean(Hepatocytes), Immune=mean(Immune),Myeloid=mean(Myeloid),
            others=mean(others)) -> df

df$x <- rep(c(2,4,6,8,10),times=2)  # 排图的坐标
df$y <- rep(c(2,4),each=5)
p <- ggplot() + 
  geom_scatterpie(data=df,
    aes(x,y,group=bin,r=0.9),
    cols = c("Bcells", "Endothelial", "Fibroblast", "Hepatocytes", "Immune", "Myeloid", "others"))+
  labs(title="percentage of infer fibroblast with different levels in fibroblast spot")
ggsave(p, file="C:/Users/DELL/Desktop/fib_bin.pdf",
       height = 3,width = 9)

# 各bin的数目
deconv %>% 
  group_by(bin) %>%
  summarise(n=n()) -> plot_df
plot_df[10,2] <- 7 # 原为1，画图看不出来
ggplot() + geom_bar(data =plot_df, aes(x=bin,y=n),stat = "identity")
