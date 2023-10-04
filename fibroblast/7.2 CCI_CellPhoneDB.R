## CellPhoneDB——下游分析即可视化
## jingsiyu
## 2022.05.05


# 下游分析1：count数目-------------------------------------------------------------------------

setwd("F:/Ajingsiyu/PROJECT/HeYuFei-TME/RESULT/2.ST/neigborhood/2.cellphoneDB/")
# 加载文件1
significant_means <- read.table("neigbor/output/significant_means.txt", fill=T,header = T,sep="\t")
significant_means <- significant_means[c(2,13:ncol(significant_means))]

significant_means <- significant_means %>%
  mutate(missing =rowSums(is.na(significant_means[-1]))) %>%
  filter(missing < ncol(significant_means)-1) # 删除全部为NA的

colnames(significant_means) <- colnames(significant_means) %>%
  gsub("Bile.duct.cells","Bile_duct_cells",.) %>%
  gsub("cancer.Hepatocyte","cancer_Hepatocyte",.)
count_df <- data.frame(id = colnames(significant_means)[2:65],
                       num = colSums(!is.na(significant_means[2:65])))

count_df$from <- apply(count_df,1,function(x){strsplit(x[['id']],'\\.')[[1]][1]})
count_df$to <- apply(count_df,1,function(x){strsplit(x[['id']],'\\.')[[1]][2]})
count_df <- count_df[-1]

count_mat <- unstack(count_df[c(1,3)])
rownames(count_mat) <- count_df$to[1:8]
count_mat <- as.matrix(count_mat)
  
#画图
color_high <- colorRampPalette(c("#66CDAA","#FFC125"))(50)
color_low <- colorRampPalette(c("#000093","#66CDAA"))(100)
my_color <- c(color_low,color_high)
count_mat2 <- count_mat[-1,-1]
pheatmap(count_mat2,
         color = my_color,
         display_numbers=T,
         number_format = "%.0f")


# 下游分析：输出TSF比fib多的显著l-r表格-------------------------------------------------------------------------

library(stringr)
setwd("F:/Ajingsiyu/PROJECT/HeYuFei-TME/RESULT/2.ST/neigborhood/2.cellphoneDB/")

# 1. sig_mean 表格
means <- read.table("neigbor/output/significant_means.txt", fill=T, header=T, sep="\t", row.names = 2)

#2. fib 的相互作用
id_fib <- str_detect(colnames(means), "Fibroblast")
means_fib <- means[id_fib]  # 所有TSF 相关的l-r pairs

pair_fib <- means_fib %>%  # 删去全为0
  mutate(missing=rowSums(is.na(means_fib))) %>%
  filter( missing<15 )%>%
  row.names()

means_specific <- means_TSF %>% filter(!(row.names(.) %in% pair_fib))
write.csv(means_specific,file="TSF_specific.csv")

# TSF-cancer pairs (all and specific) -------------------------------------------------------------------------
significant_means <- read.table("ALL/output/significant_means.txt",
                                fill=T,header = T,sep="\t")
significant_means <- significant_means[c(2,13:ncol(significant_means))]

# 总共的
TSF_cancer <- significant_means %>%  # 43个
  filter( TSF_high.cancer.Hepatocyte >0 )
write.csv(TSF_cancer,file="C:/Users/54356/Desktop/TSF_cancer_all_means.csv",row.names = F) 

cancer_TSF <- significant_means %>%  #137个
  filter( cancer.Hepatocyte.TSF_high >0)
write.csv(cancer_TSF,file="C:/Users/54356/Desktop/cancer_TSF_all_means.csv",row.names = F) 

# 特异的
TSF_cancer <- significant_means %>%  # 43个
  filter( TSF_high.cancer.Hepatocyte >0 , is.na(Fibroblast.cancer.Hepatocyte))

cancer_TSF <- significant_means %>%  #137个
  filter( cancer.Hepatocyte.TSF_high >0 , is.na(Hepatocyte.TSF_high))
write.csv(TSF_cancer,file="C:/Users/54356/Desktop/TSF_cancer_means.csv",row.names = F) 