
#预后模型的临床应用
rm(list = ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/预后模型的临床应用-8.3/")

load(file = "../合并数据集(临床样本)-1.4/合并表达矩阵.Rdata")
load(file = "../合并数据集(临床样本)-1.4/合并的临床数据.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

exprSet <- as.data.frame(t(merge_dat[gene$sig_gene_multi_cox,]))
exprSet$RiskScore <- exprSet[,rownames(gene)[1]]*gene[1,2]+
  exprSet[,rownames(gene)[2]]*gene[2,2]+
  exprSet[,rownames(gene)[3]]*gene[3,2]+
  exprSet[,rownames(gene)[4]]*gene[4,2]+
  exprSet[,rownames(gene)[5]]*gene[5,2]+
  exprSet[,rownames(gene)[6]]*gene[6,2]+
  exprSet[,rownames(gene)[7]]*gene[7,2]+
  exprSet[,rownames(gene)[8]]*gene[8,2]+
  exprSet[,rownames(gene)[9]]*gene[9,2]+
  exprSet[,rownames(gene)[10]]*gene[10,2]+
  exprSet[,rownames(gene)[11]]*gene[11,2]
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore) , "Low","High")
table(exprSet$RiskGroup)
rownames(merge_phe) <- merge_phe$Id
exprSet <- exprSet[rownames(merge_phe),]
merge_phe$RiskGroup <- exprSet$RiskGroup
merge_phe$RiskScore <- exprSet$RiskScore

phe <- merge_phe
phe$OS <- ifelse(phe$OS == "Alive",0 ,1)
table(is.na(phe$OS.time))
library(ggstatsplot)
library(tidyverse) 
library(ggpubr)
datSet <- phe

datSet$group <- ifelse(datSet$group == "normal", 'Normal','Tumor')
table(datSet$gender)
datSet$gender <- ifelse(datSet$gender == "female" , "Female", "Male")
datSet$Age <- ifelse(datSet$Age <= 60 , '<=60','>60')
table(datSet$Age)

#年龄---------
table(is.na(datSet$Age))
# datSet$Age <- ifelse(datSet$Age <= 60 , '<=60','>60')
table(datSet$Age)
p2 = ggbetweenstats(data = datSet, x = Age,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric",#选择非参数检验
                    package = "ggsci",#调用调色板
                    palette = "lanonc_lancet",
                    title = 'Age')+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
p2
ggsave(filename = 'Age.pdf',plot = p2, width = 8,height = 6,path = '../../Fig/Step5-2/')

# Stage------
datSet$Stage <- datSet$tnm_stage
table(datSet$Stage)
table(datSet$Stage != 'not reported')
datSet5 <- subset(datSet,datSet$Stage != 'not reported')
table(datSet5$Stage)
# datSet5$Stage <- ifelse(datSet5$Stage == ' stage i','I' ,
#                         ifelse(datSet5$Stage == ' stage ii' ,'II',
#                                ifelse(datSet5$Stage == ' stage iii', 'III', 'IV')))
#table(datSet5$Stage)

p6 = ggbetweenstats(data = datSet5, x = Stage,y = RiskScore,
                    bf.message = FALSE,#去除贝叶斯相关的统计值
                    type = "nonparamatric", #选择非参数检验
                    grouping.var = Stage, # 分组变量
                    package = "ggsci", #调用调色板
                    palette = "lanonc_lancet",
                    title = 'Stage',
                    pairwise.comparisons=TRUE,
                    pairwise.display ="all")+
  theme(plot.title = element_text(size = 15,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))
# outlier.tagging = TRUE,#标记异常值
# p.adjust.method = "fdr",
# pairwise.comparisons = TRUE)
# ggpubr::stat_compare_means(comparisons=my_comparisons3, label.y = c(0.08, 0.12, 0.18, 0.22, 0.26,
#                                                                     0.30, 0.34, 0.38, 0.42, 0.52),
#                            method = "wilcox.test",
#                            label = 'p.signif') # p.signif p.format
p6
ggsave(filename = 'Tumor Stage.pdf',plot = p6, width = 8,height = 6,
       path = './')


