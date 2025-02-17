#后台
rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/ESCC/copykat-5/")
library(Seurat)
library(infercnv)
#### CopyKAT（自动化识别肿瘤）#####
load(file = "/home/datahup/syj/NSCC/step10_1copykat_CNV/Epi_output.Rdata")
table(Epi@active.ident)
counts <- as.matrix(Epi@assays$RNA@counts)
##运行时间较长 
#时间可能要两天   Time difference of 2.45014 days
library(copykat)
cnv <- copykat(rawmat=counts,
               ngene.chr=5,
               sam.name="ESCC",
               n.cores=8)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV(拷贝数变异）预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
save(cnv, "/home/datahup/syj/ESCC/copykat-5/cnv.rdata")



