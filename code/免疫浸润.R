
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/ssGSEA-9.1/')

# 获取28种免疫细胞基因集---------
library(dplyr)
library(tidyverse)
geneSet <- read.csv("CellReports.txt",header = F,sep = "\t",) # 用EXCEL打开删除NA列
class(geneSet)
geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

# Imu = as.list(c(l[["Central memory CD8 T cell"]],l[["Effector memeory CD8 T cell"]],l[["Activated CD4 T cell"]],
#                 l[["Central memory CD4 T cell"]],l[["Effector memeory CD4 T cell"]],l[["T follicular helper cell"]],
#                 l[["Gamma delta T cell"]],l[["Type 1 T helper cell"]],l[["Type 17 T helper cell"]],l[["Type 2 T helper cell"]],
#                 l[["Regulatory T cell"]],l[["Activated B cell"]],l[["Immature  B cell"]],l[["Memory B cell"]],l[["Natural killer cell"]],
#                 l[["CD56bright natural killer cell"]],l[["CD56dim natural killer cell"]],l[["Myeloid derived suppressor cell"]],
#                 l[["Natural killer T cell"]],l[["Activated dendritic cell"]],l[["Plasmacytoid dendritic cell"]],l[["Immature dendritic cell"]],
#                 l[["Macrophage"]],l[["Eosinophil"]],l[["Mast cell"]],l[["Monocyte"]],l[["Neutrophil"]]))

#Imu = c(l[["Central memory CD8 T cell"]],l[["Effector memeory CD8 T cell"]],l[["Activated CD4 T cell"]])

#准备表达矩阵-----
load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
df <- scale(dat53625)
df <- as.data.frame(df)
# 导出数据后，需要去除第一行第一列的空值，并改为.txt格式
#write.csv(df,file = "df_scale.csv")

# 运行ssGSEA--------
library(GSVA)
library(limma)
exprSet <- as.matrix(dat53625)
exprSet[1:6,1:6]
ssgsea<- gsva(exprSet, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#ssgsea<- gsva(exprSet, l,method='ssgsea',mx.diff=F,verbose = F)

Immune2 = as.data.frame(t(ssgsea))
Immune2$id <- rownames(Immune2)
save(Immune2,file = 'FPKM_ssGSEA.Rdata')

#load(file = 'FPKM_ssGSEA.Rdata')
# 加载训练集-------------
#load(file = "../机器学习和验证然后取跑100种机器学习-8.2/Train_Result.Rdata") 
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
#df <- read.csv(file = 'df_scale.csv')
#rownames(df) <- df$X
#df <- df[,-1]
df[1:6,1:6]
exprSet <- as.data.frame(t(df[rownames(gene),]))
dim(exprSet)
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
#exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore) , "Low","High")


exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore), "Low", "High")
table(exprSet$RiskGroup)
exprSet$id <- rownames(exprSet)
table(phe53625$group)
phe53625 <- subset(phe53625,phe53625$group == 'cancer')
#exprSet <- subset(exprSet,exprSet$Type == 'Tumor')
exprSet <- exprSet[rownames(phe53625),]

ImmSet <- merge(Immune2,exprSet,by = "id")
rownames(ImmSet) <- ImmSet$id
colnames(ImmSet)
ImmSet2 <- ImmSet[,2:29]
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
dim(ImmSet2)
dim(ImmSet3)
ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == ImmSet$id ,ImmSet$RiskGroup,0)
table(ImmSet3$RiskGroup)
table(ImmSet$RiskGroup)
save(ImmSet,ImmSet3,Immune2,file = 'ssGSEA_output.Rdata')

# 小提琴图----------
library(ggplot2)
library(tidyverse)
library(ggpubr)
# plot
colnames(ImmSet3)
ggplot(ImmSet3,aes(x = reorder(`Immune Cell Type`,-`Estimating Score`) , y = `Estimating Score`, fill = RiskGroup)) +
  # # 小提琴图层
  # geom_violin(position = position_dodge(0.9),alpha = 1.2,
  #             width = 1.8,trim = T,
  #             color = NA) +
  # 箱线图图层
  geom_boxplot(width = 0.35,show.legend = F,
               position = position_dodge(0.9),
               color = 'black',alpha = 1.2,
               outlier.shape = 21) +
  # outlier.shape = 21,
  # outlier.shape = NA
  # 主题调整
  theme_bw(base_size = 16) +
  labs(x = "Immune Cell Type", y = 'ssGSEA Estimating Score') +
  theme(axis.text.x = element_text(angle = 65,hjust = 1,color = 'black'),
        legend.position = 'top',
        aspect.ratio = 0.4) +
  # 颜色设置
  # scale_fill_manual(values = c('Low'='#398AB9','High'='red'),
  #                   name = '') +
  scale_fill_manual(values = c("#FF0000",'#00CCFF'))+ 
  # 添加显著性标记
  stat_compare_means(aes(group = RiskGroup,label = ..p.signif..),#..p.format..   ..p.signif..
                     method = "wilcox.test")#kruskal.test  p.signif wilcox.test
# #添加显著性标记
# stat_compare_means(aes(group=RiskGroup),
#                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
# 
#                                     symbols = c("***", "**", "*", "NS")),label = "p.signif",
#                    label.y = 1.7,size = 4) + ylim(0.3,1.7)
ggsave(filename = './11.14结果/ssGSEA.pdf',width = 16,height = 10)

#美化(16*10)-----
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
  theme(axis.text.x = element_text(angle = 65, hjust = 1, color = 'black'),
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



































