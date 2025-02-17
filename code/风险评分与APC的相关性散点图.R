#风险评分与抗原呈递细胞的相关性散点图
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/ssGSEA-9.2/')

# 获取免疫细胞基因集---------
library(dplyr)
library(tidyverse)
geneSet <- read.csv("整理.txt", header = F, sep = "\t")

geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")

df <- scale(dat53625)
df <- as.data.frame(df)

# 运行ssGSEA--------
library(GSVA)
library(limma)
exprSet <- as.matrix(dat53625)
exprSet[1:6,1:6]
ssgsea<- gsva(exprSet, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
#ssgsea<- gsva(exprSet, l,method='ssgsea',mx.diff=F,verbose = F)

Immune2 = as.data.frame(t(ssgsea))
Immune2$id <- rownames(Immune2)

#分组--------
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

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

exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore), "Low", "High")
table(exprSet$RiskGroup)
exprSet$id <- rownames(exprSet)
table(phe53625$group)
phe53625 <- subset(phe53625,phe53625$group == 'cancer')
exprSet <- exprSet[rownames(phe53625),]

ImmSet <- merge(Immune2,exprSet,by = "id")
rownames(ImmSet) <- ImmSet$id
colnames(ImmSet)
ImmSet2 <- ImmSet[,2:13]

# library(dplyr)
# library(tidyr)
# library(tibble)
# library(ggplot2)
# ImmSet3 <- ImmSet2 %>% as.data.frame() %>%
#   rownames_to_column("Sample") %>% 
#   gather(key = 'Immune Cell Type',value = 'Estimating Score',-Sample)
# dim(ImmSet2)
# dim(ImmSet3)
# ImmSet3$RiskGroup <- ifelse(ImmSet3$Sample == ImmSet$id ,ImmSet$RiskGroup,0)
# table(ImmSet3$RiskGroup)
# table(ImmSet$RiskGroup)

#save(ImmSet,file = "./ImmSet.Rdata")
# 相关性分析---------
library(ggplot2)
library(ggstatsplot)
library(cowplot)
colnames(ImmSet)
Immcell <- colnames(ImmSet)[2:13]
# 检查是否符合正态分布
for (i in 1:length(Immcell)) {
  x=i+1
  a <- shapiro.test(ImmSet[,x])
  ifelse(a[["p.value"]] > 0.05,
         print(paste0(Immcell[i],' is Yes')),
         print(paste0(Immcell[i],' is No')))
}
plotlist <- list()
library(rlang)
for (i in 1:length(Immcell)) {
  print(paste0('Now is ',Immcell[i]))
  xlab1 <- sym(Immcell[i])
  p1 <- ggscatterstats(data = ImmSet,
                       y = RiskScore,
                       x = !!xlab1,
                       # centrality.para = "mean",
                       # margins = "both",
                       bf.message = FALSE,#去除贝叶斯相关的统计值
                       type = "nonparamatric",#选择非参数检验
                       xfill = "#CC79A7",
                       yfill = "#009E73",
                       # marginal.type = "density", # 其中marginal.type可选 histograms，boxplots，density，violin，densigram (density + histogram)
                       title = paste0("Relationship between RiskScore and ",Immcell[i]))
  # print(p1)
  # ggsave(filename = paste0(Immcell[i],'.png'), plot = p1, width = 8, height = 6, path = '../Figer/Relationship')
  plotlist[[i]] = p1
}
#保存600*1500
plot_grid(plotlist[[8]],plotlist[[6]],plotlist[[7]]
         ,nrow = 1)#,plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],

































































