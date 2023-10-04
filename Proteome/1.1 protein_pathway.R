## 蛋白质组 概况、通路
## jingsiyu 2022.01

library(ggplot2)
library(dplyr)

#1. pca -------------------------------------------------------------------------

# 加载数据
df <- read.csv("/3.蛋白/1.pca/all_sample.csv")
rownames(df) <- df$Protein
df <- data.frame(t(df[c(-1:-3)]))
# 过滤
drop_info <- apply(df, 2, function(x){sum(x ==0)})
summary(drop_info)
df0 <- df[-which(drop_info > 3)] # 删去在3个样本中都是0表达的蛋白

# 删去T-S7
df <- df0[-c(22,24),]
# way1
data.pca <- prcomp(df,scale = T)
sum <- summary(data.pca)

pca_sample <- data.frame(data.pca[["x"]][ ,1:2])
colnames(pca_sample) <-c ('Dim.1','Dim.2')

pca_eig1 <- round(sum[["importance"]][1,1],digits = 2)
pca_eig2 <- round(sum[["importance"]][1,2],digits = 2)

# 整理数据画图
pca_sample$type <- c(rep(c("N-P","T-P","N-S","T-S"),times=5),"N-P","T-P")
pca_sample$samples <- rownames(pca_sample)
head(pca_sample) 

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = type), size = 3) +  
  scale_color_manual(values = c('DeepPink','RosyBrown2','DarkTurquoise','SkyBlue')) +  
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  

p + stat_ellipse(aes(fill = type), geom = 'polygon', level = 0.95, alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c('DeepPink','RosyBrown2','Turquoise','SkyBlue'))

pdf("C:/Users/DELL/Desktop/pca.pdf")
dev.off()

#2. 蛋白质组-间质 通路富集 -------------------------------------------------------------------------

df_ST <- read.csv("C:/Users/DELL/Desktop/S_T高.csv")
df_SN <- read.csv("C:/Users/DELL/Desktop/S_N高.csv")
df_HT <- read.csv("C:/Users/DELL/Desktop/H_T高.csv")
df_HN <- read.csv("C:/Users/DELL/Desktop/H_N高.csv")
df <- df_HT

# GO
df <- df %>% 
  filter(GO_Class == "BP") %>% 
  filter(AdjustedPv < 0.05) %>% 
  arrange(AdjustedPv)
df_HT <- df

save(df_HN,df_HT,df_SN,df_ST,file="C:/Users/DELL/Desktop/protome_GO_BP.RData")
setdiff(df_ST$GO_Term,df_HT$GO_Term) # ST中独有的

#3. S-T_N 通路气泡图 -------------------------------------------------------------------------
# CSV文件是公司的分析结果，D vs C上调下调的

df1 <- read.csv("GO_down.csv")
df1$type <-"N-S"
df2 <- read.csv("GO_up.csv")
df2$type <-"T-S"
df <- rbind(df1,df2)

GO <- df %>%
  mutate(gene_ratio = x/n) %>%
  filter(GO_Class == "BP", AdjustedPv < 0.05) %>%
  group_by(type) %>%
  top_n(10, wt=gene_ratio) %>%
  select(GO_ID,GO_Term,gene_ratio,x,AdjustedPv,type)
colnames(GO) <- c("ID","Term","gene_ratio","counts","q_value","type")

ggplot(plot_df,aes(x=type,y=Term,size=gene_ratio,colour=q_value2)) +
  geom_point()+
  scale_color_gradient(low="#FFE4C4",high="#A52A2A",
                      breaks=seq(1,15,by=3))+
  theme_bw()

pdf("S_T&N.pdf",height = 5,width = 7)
dev.off()

