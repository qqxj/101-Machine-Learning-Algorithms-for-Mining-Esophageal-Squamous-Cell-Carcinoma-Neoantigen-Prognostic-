#单因素cox
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/机器学习和验证然后取跑100种机器学习-8.2/")

load(file = "../100种机器学习算法-8.1/验证集数据格式.Rdata")
load(file = "/home/datahup/syj/ESCC/取交集韦恩图-7/交集结果.Rdata")
datset2$Id
tcga <- subset(datset2, grepl("^TCGA", Id))
rownames(datset2) <- datset2$Id
geo <- datset2[c(1:119),]
save(geo,tcga,file = "拆分数据.Rdata")







































































