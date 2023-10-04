## 建立分类器选择更关键的marker gene
## 2022.04.13

library(randomForest)
library(pROC)
library(dplyr)

set.seed(2022)

#1. 准备数据 -------------------------------------------------------------------------

load(expr_SCT.RData) # expr
my_table <- readRDS("D:/资料/metadata_0.5.rds") # metadata
genes <- c("COL1A1","CTGF","FSTL1","COL4A1","COL4A2")

expr <- expr[genes,]
score_df <- expr %>% t() %>% as.data.frame() %>%
  mutate(spot=rownames(.))
score_df <- merge(score_df,my_table[,c("spot","integrated_snn_res.0.5")],by="spot")
score_df <- score_df %>%
  mutate(group=as.factor(ifelse(score_df$integrated_snn_res.0.5 ==5,1,0))) %>%
  group_by(group) %>% 
  sample_n(376)

train_sub <- sample(nrow(score_df),0.8*nrow(score_df))
train_data <- score_df[train_sub,]
test_data <- score_df[-train_sub,]

#2. 随机森林（全部基因） -------------------------------------------------------------------------
gene_rf <- randomForest(group ~ COL4A2+CTGF+FSTL1+COL1A1+COL4A1, 
                        data=train_data, 
                        importance=TRUE,
                        proximity=TRUE)
print(gene_rf)
varImpPlot(gene_rf, main="variable importance") #变量贡献度
prob <- as.data.frame(gene_rf[["votes"]]) #概率值

# 训练集上的ROC
ROC <- roc(train_data$group, prob$`1`,ci=T,auc=T)
plot(ROC, legacy.axes = TRUE,
     print.auc=T,
     max.auc.polygon=T,
     print.thres=T,
     auc.polygon=T,
     auc.polygon.col = "#fca311")

#预测：ROC
pred <- predict(gene_rf,newdata = test_data)
table(test_data$group, pred) #混淆矩阵

pred <- predict(gene_rf,newdata = test_data,type="prob")
ROC <- roc(test_data$group, as.data.frame(pred)[,2],ci=T,auc=T)
plot(ROC, legacy.axes = TRUE,
     print.auc=T,
     max.auc.polygon=T,
     print.thres=T,
     auc.polygon=T,
     auc.polygon.col = "#fca311"
)

#3. 随机森林（最终基因） -------------------------------------------------------------------------

gene_rf1 <- randomForest(group ~ CTGF, 
                         data=train_data,proximity=TRUE)
gene_rf2 <- randomForest(group ~ CTGF+COL4A2, 
                         data=train_data,proximity=TRUE)
gene_rf3 <- randomForest(group ~ CTGF+COL4A2+FSTL1, 
                         data=train_data,proximity=TRUE)

# 训练集ROC
ROC1 <- roc(train_data$group, 
            as.data.frame(gene_rf1[["votes"]])[,2],
            ci=T,auc=T)
ROC2 <- roc(train_data$group, 
            as.data.frame(gene_rf2[["votes"]])[,2],
            ci=T,auc=T)
ROC3 <- roc(train_data$group, 
            as.data.frame(gene_rf3[["votes"]])[,2],
            ci=T,auc=T)

# 画图 -------------------------------------------------------------------------

roc_obj3 <- plot(ROC3,col ="#ff9f1c",print.auc=T)
roc_obj1 <- lines.roc(ROC1, col ="#2EC4B6") 
roc_obj2 <- lines.roc(ROC2,col ="#E71D36")
roc_obj4 <- lines.roc(ROC,col ="#0d3b66")

legend("bottomright",legend = c("CTGF","CTGF+COL4A2","CTGF+COL4A2+FSTL1","all genes"),
       col = c("#2EC4B6","#E71D36","#ff9f1c","#0d3b66"),
       lwd=2)