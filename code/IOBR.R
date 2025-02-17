
#b站：IOBR（功能最强的免疫细胞浸润分析R包？）
#https://blog.csdn.net/zfyyzhys/article/details/141304607
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/ssGSEA-9.1/')
library(IOBR)
library(tidyverse)

#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
my_data <- dat53625

#查看参数和算法------
tme_deconvolution_methods
signature_score_calculation_methods
#肿瘤微环境相关特征基因集（√）
names(signature_tme)
#代谢相关基因集
names(signature_metabolism)
#与生物医学基础相关的基因集
names(signature_tumor)
#所有免疫细胞相关的基因集（√）
names(signature_collection)
#查看参考基因集的文献出处
signature_collection_citation

signature_tme$CD_8_T_effector

#cibersort分析-----
cibersort <- deconvo_tme(eset = my_data,
                         method = "cibersort",
                         arrays = F,
                         perm = 100)

# cell_bar_plot(input = cibersort,
#               title = "cell FACTION",
#               legend.position = "bottom",
#               palette = 3,#调色
#               coord_filp = T,#是否反转柱状图
#               show_col = F,
#               pattern = TG_phe$group)#c("cell_type1", "cell_type2")
ciber_plot <- cell_bar_plot(input = cibersort[1:20,], 
                            features = colnames(cibersort)[2:23], 
                            title = "CIBERSORT Cell Fraction",)

ggsave("cibersort.pdf",width = 10,height = 8)

#xcell分析-------
xcell <- deconvo_tme(eset = my_data,
                     method = "xcell",
                     arrays = F)

#epic分析-------
epic <- deconvo_tme(eset = my_data,
                     method = "epic",
                     arrays = F)

#timer分析-------
  ###报错往下看###
timer <- deconvo_tme(eset = my_data,
                    method = "timer",
                    arrays = F,
                    indications = id)

#这里导入的group_list必须是TIMER能够识别的，比如TCGA中33肿瘤类型
phe53625$group <- ifelse(phe53625$group == "cancer","ESCA","ESCA")
TIMER <- deconvo_tme(eset = my_data,
                     method = "timer",
                     group_list = phe53625$group)

load(file = "/home/datahup/syj/ESCC/100种机器学习算法-8.1/GEO训练集.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
phe <- subset(phe53625,group == "cancer")
rownames(phe) #<- phe$Id
rownames(TIMER) <- TIMER$ID
TIMER1 <- TIMER[rownames(phe),]

save(TIMER1,TIMER,file = "TIMER.Rdata")


#准备数据
load(file = "TIMER.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
df <- as.data.frame(scale(dat53625))
exprSet <- as.data.frame(t(df[rownames(gene),]))
exprSet$RiskScore <- exprSet[,rownames(gene)[1]]*gene[1,2]+
  exprSet[,rownames(gene)[2]]*gene[2,2]+
  exprSet[,rownames(gene)[3]]*gene[3,2]+
  exprSet[,rownames(gene)[4]]*gene[4,2]+
  exprSet[,rownames(gene)[5]]*gene[5,2]#+
  # exprSet[,rownames(gene)[6]]*gene[6,2]+
  # exprSet[,rownames(gene)[7]]*gene[7,2]+
  # exprSet[,rownames(gene)[8]]*gene[8,2]+
  # exprSet[,rownames(gene)[9]]*gene[9,2]+
  # exprSet[,rownames(gene)[10]]*gene[10,2]+
  # exprSet[,rownames(gene)[11]]*gene[11,2]
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore), "Low", "High")
table(exprSet$RiskGroup)
exprSet$id <- rownames(exprSet)
#table(TG_phe$group)
#TG_phe <- subset(TG_phe,TG_phe$group == 'cancer')
#exprSet <- exprSet[TG_phe$Id,]
#table(exprSet$RiskGroup)

names(TIMER1)[1] <- "id"
#TIMER1 <- TIMER[rownames(phe53625),]
ImmSet <- merge(TIMER1,exprSet,by = "id")
colnames(ImmSet)
names(ImmSet)[2:7] <- c( "B cell", "CD4 T cell", "CD8 T cell", "Neutrophil", "Macrophage", "Dendritic Cell")
rownames(ImmSet) <- ImmSet$id
colnames(ImmSet)
ImmSet2 <- ImmSet[2:7]

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == ImmSet$id ,ImmSet$RiskGroup,0)
table(ImmSet3$RiskGroup)
table(ImmSet$RiskGroup)
save(ImmSet,ImmSet3,ImmSet2,file = 'ssGSEA_output_IOBR.Rdata')

# 小提琴图----------
library(ggplot2)
library(tidyverse)
library(ggpubr)
# plot
colnames(ImmSet3)
ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  theme_bw(base_size = 16) +
  labs(x = "Immune Cell Type", y = 'Immune infiltration') +
  theme(axis.text.x = element_text(angle = 65,hjust = 1,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  scale_fill_manual(values = c("#FF0000",'#00CCFF'))+ 
  stat_compare_means(aes(group = RiskGroup,label = ..p.format..),#..p.format..  ..p.signif..
                     method = "wilcox.test")#kruskal.test  p.signif wilcox.test

#ggsave(filename = '免疫浸润2.pdf',width = 8,height = 5)

#美化-----
ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  theme_bw(base_size = 16) +
  labs(x = " ", y = 'Immune infiltration') + #x = "Immune Cell Type"
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  scale_fill_manual(values = c("#FF5733",'#33FF57'))+ 
  stat_compare_means(aes(group = RiskGroup,label = ..p.signif..),#..p.format..  ..p.signif..
                     method = "wilcox.test")#kruskal.test  p.signif wilcox.test


ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  # 小提琴图层
  geom_violin(position = position_dodge(0.9),
              alpha = 1.2,
              width = 1,trim = F,
              color = NA) +
  # 箱线图图层
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  theme_bw(base_size = 16) +
  labs(x = " ", y = 'Immune infiltration') + #x = "Immune Cell Type"
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  scale_fill_manual(values = c("#FF5733",'#33FF57'))+ 
  stat_compare_means(aes(group = RiskGroup,label = ..p.signif..),
                     method = "wilcox.test")#kruskal.test  p.signif wilcox.test

#美化
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
  labs(x = "Immune Cell Type", y = 'ssGSEA Estimating Score') +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4,
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
        axis.title.x = element_text(margin = margin(t = 10))) +  # 增加 x 轴标题的上边距
  # 颜色设置
  scale_fill_manual(values = c("#FF0000", 'green'), name = 'Risk Group') +
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup, label = ..p.signif..), method = "wilcox.test") +
  theme(panel.grid.major = element_line(color = "lightgrey", size = 0.5),
        panel.grid.minor = element_blank()) #+
# # 调整 x 轴的宽度
# scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))  # 适当调整间距
