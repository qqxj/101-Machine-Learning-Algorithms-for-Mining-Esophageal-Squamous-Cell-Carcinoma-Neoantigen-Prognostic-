#高低风险与ICD和ICP调节之间的关联
rm(list = ls())
options(stringsAsFactors = F)

setwd("/home/datahup/syj/NSCC/step6_1TCGA_GEO/")
load(file = "./TCGA_output.Rdata")

setwd('/home/datahup/syj/ESCC/ICP、ICD-9.3/')
# load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
# load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

#ICD:
# table(rownames(dat_2) %in% c("TSPAN6","HGF","EIF2AK2","EIF2AK1","P2RX7","MET","PANX1","LRP1","EIF2AK4",
#                              "ANXA1","TLR4","IFNAR2","EIF2A","IFNAR2","TLR3","CXCL10","FPR1","EIF2AK3",
#                              "P2RY2","IFNW1","CALR","IFNE","HMGB1","IFNA1"))

ICD <- c("CALR","CXCL10","FPR1","HGF","IFNAR1","IFNAR2","IFNB1","IFNE","IFNK","IFNW1","MET",
         "P2RX7","TLR4","LRP1","EIF2A","EIF2AL3","EIF2AK2","EIF2AK4","EIF2AK1","HMGB1","ANXA1","PANX1",
         "P2RY2","IFNA1","IFNA2","TLR3")
table(rownames(dat_2) %in% ICD)

#ICP
# table(rownames(dat_2) %in% c("TSPAN6","CD44","TNFRSF9","LAG3","CD200","NRP1","CD40","CD40LG",
#                              "CD276","CD86","HHLA2","CD48","CD160","TNFSF4","CD274","TNFSF18",
#                              "TNFRSF8","CD80","CD244","TNFSF9","CD70","TNFSF14","ADORA2A","VTCN1","IDO1",
#                              "HAVCR2","TNFRSF14","CD27","ICOSLG","CTLA4","ICOS","CD200R1","LAIR1","KIR3DL1",
#                              "TMIGD2","LGALS9","CD28","TNFSF15","TIGIT","BTLA","TNFRSF4","TNFRSF18","PDCD1","PDCD1LG2",
#                              "BTNL2","IDO2","TNFRSF25"))
ICP <-  c("ADORA2A","BTLA","BTNL2","CD160","CD200","CD200R1","CD244","CD27","CD274","CD276",
                             "CD28","CD40","CD40LG","CD44","CD48","CD70","CD80","CD86","CTLA4","HAVCR2","HHLA2",
                             "ICOS","ICOSLG","IDO1","IDO2","KIR3DL1","LAG3","LAIR1","LGALS9","NRP1","PDCD1","PDCD1LG2",
                             "TIGIT","TMIGD2","TNFRSF14","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14","TNFSF15",
                             "TNFSF18","TNFSF4","TNFSF9","VSIR","VTCN1")
table(rownames(dat_2) %in% ICP)

phe53625 <- phe6
dat53625 <- dat_2

phe53625 <- subset(phe53625,SubType == "ESCC")
dat53625 <- as.data.frame(t(dat53625[c(rownames(gene),ICD),]))
colnames(dat53625)
dat53625 <- dat53625[,-21]

dat53625 <- as.data.frame(t(dat53625[c(rownames(gene),ICP),]))
dat53625 <- dat53625[phe53625$Id,]
colnames(dat53625)
dat53625 <- dat53625[,-51]

exprSet <- dat53625
exprSet$RiskScore <- exprSet[,rownames(gene)[1]]*gene[1,2]+
  exprSet[,rownames(gene)[2]]*gene[2,2]+
  exprSet[,rownames(gene)[3]]*gene[3,2]+
  exprSet[,rownames(gene)[4]]*gene[4,2]+
  exprSet[,rownames(gene)[5]]*gene[5,2]
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore), "Low", "High")
table(exprSet$RiskGroup)
exprSet <- exprSet[, !(colnames(exprSet) %in% gene$sig_gene_multi_cox)]

ImmSet <- exprSet
colnames(ImmSet)
ImmSet2 <- ImmSet[,1:25]

colnames(ImmSet)
ImmSet2 <- ImmSet[,1:46]


library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == rownames(ImmSet) ,ImmSet$RiskGroup,0)
table(ImmSet3$RiskGroup)
table(ImmSet$RiskGroup)

library(ggpubr)
#ICD:16*10
#ICP:
ggplot(ImmSet3, aes(x = reorder(`Immune Cell Type`, -`Estimating Score`), 
                    y = `Estimating Score`, fill = RiskGroup)) +
  # 小提琴图层
  geom_violin(position = position_dodge(0.9), alpha = 0.7, width = 1.25, trim = TRUE) +
  # 箱线图图层
  geom_boxplot(width = 0.5, show.legend = FALSE, 
               position = position_dodge(0.9), 
               color = 'black', alpha = 0.8, 
               outlier.shape = 21, outlier.colour = "black") +
  # 主题调整
  theme_bw(base_size = 16) +
  labs(x = " ", y = 'Gene Expression') +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4,
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic")) +  # 增加 x 轴标题的上边距
  # 颜色设置
  scale_fill_manual(values = c("#FF0000", 'green'), name = 'Risk Group') +
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup, label = ..p.signif..), method = "wilcox.test") +
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.5),
        panel.grid.minor = element_blank()) 

#制箱线图
# # 加载必要的包
# library(ggplot2)
# library(ggsignif)
# 
# # 使用 ggplot 绘制箱线图并添加 p 值
# ggplot(exprSet, aes(x = RiskGroup, y = PDCD1, fill = RiskGroup)) +
#   geom_boxplot(outlier.color = "red", outlier.size = 2) +  # 绘制箱线图，标记异常值
#   geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # 添加散点
#   labs(
#     title = "Expression by RiskGroup",
#     x = "Risk Group",
#     y = "Expression"
#   ) +
#   theme_minimal() +  # 使用简洁主题
#   scale_fill_manual(values = c("Low" = "blue", "High" = "red")) +  # 自定义颜色
#   geom_signif(comparisons = list(c("Low", "High")), 
#               map_signif_level = TRUE,  # 自动显示 p 值
#               textsize = 5)  # 调整 p 值的文本大小






























































