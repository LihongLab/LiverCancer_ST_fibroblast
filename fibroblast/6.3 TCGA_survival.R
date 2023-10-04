## TCGA 生存分析
## 2022.1.4

library(dplyr)
library(tidyverse)
library(survival)
library(survminer)

load("F:/Ajingsiyu/PROJECT/HeYuFei-TME/RESULT/ST/202111/survival/TCGA-LIHC_meta_clinical.RData")

# 整理临床变量
meta_clinical <- meta_clinical %>% 
  filter(ajcc_pathologic_stage %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IV","Stage IVA","Stage IVB")) 

meta_clinical$ajcc_stage <- as.numeric(factor(meta_clinical$ajcc_pathologic_stage))
meta_clinical[which(meta_clinical$ajcc_stage %in% c(3,4,5,6)),"ajcc_stage"] <- 3
meta_clinical[which(meta_clinical$ajcc_stage %in% c(7,8,9)),"ajcc_stage"] <- 4
meta_clinical %>% drop_na(ajcc_stage) -> meta_clinical
summary(as.factor(meta_clinical$ajcc_stage))

genes <- c("COL1A2","COL4A2","CTGF","FSTL1")


#1.  单基因cox -------------------------------------------------------------------------
genes <- c("COL1A2","COL4A2","CTGF","COL3A1", "COL1A1","COL4A1","EFEMP1","FSTL1")
meta_stage <- meta_clinical

for (gene in genes){
  meta_stage$gene <- meta_stage[,gene]
  coxfit <- coxph(Surv(time,event)~gene,data=meta_stage)
  print(gene)
  print(summary(coxfit))
}


#2. 组合基因-中位数阈值-------------------------------------------------------------------------

load("F:/Ajingsiyu/PROJECT/HeYuFei-TME/RESULT/ST/202111/survival/TCGA-LIHC_meta_clinical.RData")

meta_stage <- meta_clinical %>% # 加权gene
  mutate(risk_score = 0.0011754*CTGF + 0.005424*COL4A2 + 0.0004091*FSTL1)

(threshould <- median(meta_stage$risk_score)) # 中位数阈值
meta_stage$risk_score2 <- ifelse(meta_stage$risk_score > threshould,1,0) #分组

#OS
sfit <- survfit(Surv(time, event)~risk_score2, data=meta_stage)
print(sfit)
ggsurvplot(sfit,
           conf.int=F, pval=TRUE, risk.table = TRUE,
           #title = "TCGA-LIHC: risk_score= weighted sum",
           font.main=18, font.x=16, font.y=16,
           palette = c("#fca311", "#1a759f"),
           legend.labs = c("risk score=0", "risk score=1"),
           legend.title = "",
           legend = c(0.8,0.75)
           )
