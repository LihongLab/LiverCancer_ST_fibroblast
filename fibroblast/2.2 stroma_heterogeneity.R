## 间质异质性分析
## 2022.05.05

readRDS("merge.rds") # 所有间质spot

MVGs <- VariableFeatures(merge)
expr <- as.matrix( merge$SCT@scale.data[MVGs,] )
cor <- cor(expr ,method = "pearson")

saveRDS(cor,file="cor.rds")

# 所有S 聚类 -------------------------------------------------------------------------

#聚类
cor <- readRDS("cor.rds")
my_table <- readRDS("margeS_matadata.rds")

my_table <- my_table[,c(8,10)] #？
my_table$spot <- rownames(my_table) #？

d2 <- as.dist(1-cor) #距离矩阵
tree <- hclust(d2, method = "complete", members = NULL) # 聚类结果

# 画图linux
order <- tree[["order"]]  # 聚类顺序
order_barcode <- rownames(cor)[order]

rownames(df) <- df$spot
df <- df[order_barcode,]
df_tmp1 <- df[1:3226,] %>% arrange(`4_position`)
df2 <- rbind(df_tmp1,df[3227:7362,])
order_barcode2 <- rownames(df2)

plot_cor <- cor[order_barcode2,order_barcode2]

annotation_row <- df %>% select(`4_position`)  # 注释位置
#annotation_row <- df %>% select(sample_id)  # 注释病人

color_high <- colorRampPalette(c("#faf9f9","#AB0002"))(50)
color_low <- colorRampPalette(c("#000093","#faf9f9"))(50)
my_color <- c(color_low,color_high)
ann_colors = list(`4_position`= c(T_S="#00BFC4",N_S="#F8766D"))

# gap
df <- df[rownames(plot_cor),]
summary(as.factor(df$cuttree_3))  # c(3226,2044,2092)
summary(as.factor(df$cuttree_4))  # c(3226,2044,818,1274)

pheatmap( pmax(pmin(plot_cor,0.5),-0.5),
          color = my_color,
          cluster_rows = FALSE, cluster_cols = FALSE,
          show_colnames = FALSE, show_rownames = FALSE,
          cellwidth =0.05, 
          cellheight = 0.05,
          gaps_col = c(3226,5270),
          gaps_row = c(3226,5270),
          annotation_row = annotation_row,
          annotation_col = annotation_row,
          annotation_colors = ann_colors,
          #annotation_colors = ann_colors,
          filename = "./test.png")


# mean pearson distance --------------------------------------------------------------------

cor <- readRDS("cor.rds")
my_table <- readRDS("margeS_matadata.rds")

T_barcode <- my_table %>% filter(`4_position` == "T_S") %>% row.names()
N_barcode <- my_table %>% filter(`4_position` == "N_S") %>% row.names()

dis_T <- (1-cor[T_barcode,T_barcode])/2
dis_N <- (1-cor[N_barcode,N_barcode])/2


# 随机取100点,相互距离，500次
distri_T <- c()
for (i in 1:500){
  a <- sample(1:nrow(dis_T),100)
  dis_mean <- mean( dis_T[a,a])
  distri_T <- c(distri_T,dis_mean)
}
summary(distri_T)

distri_N <- c()
for (i in 1:500){
  a <- sample(1:nrow(dis_N),100)
  dis_mean <- median( dis_N[a,a]) # 对角线元素
  distri_N <- c(distri_N,dis_mean)
}
summary(distri_N)

plot_df <- data.frame(y=c(distri_T,distri_N),
                      x=rep(c("T","N"),each=500))
saveRDS(plot_df,file="boxplot_plot_df.rds")

ggplot(plot_df,aes(x=x,y=y,fill=x))+
  geom_boxplot()+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  geom_signif(aes(x=x,y=y),
              comparisons = list(c("T","N")),
              y_position = c(0.5),
              map_signif_level = T,
              test = "wilcox.test"
  )+
  theme_classic(base_line_size = 1)+
  labs(x="",y="Mean Pairwise Pearson Diatance")

pdf("C:/Users/DELL/Desktop/distance_box.pdf")
dev.off()
