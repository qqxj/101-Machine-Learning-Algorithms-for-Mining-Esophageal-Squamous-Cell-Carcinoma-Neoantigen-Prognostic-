#备注：tcga的数据+gse53624数据作为测试集

rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/合并数据集(临床样本)-1.4/")
load(file = "合并表达矩阵.Rdata")
load(file = "合并的临床数据.Rdata")

#先把tcga的数据拆出来
merge_phe$Id

TCGA_phe <- subset(merge_phe, grepl("^TCGA", Id))
TCGA_dat <- merge_dat[,TCGA_phe$Id]

# GSE53624--------
GSE53625_phe <- merge_phe[86:443,]
GSE53625_dat <- merge_dat[,GSE53625_phe$Id]

setwd("/home/datahup/syj/ESCC/验证集15/")
####下载bulk数据 GSE53625 
library(GEOquery)
dir.create("./GSE53624",recursive = T)
eSet <- getGEO('GSE53624', destdir="./GSE53624",
               AnnotGPL = F,
               getGPL = F)
#save(eSet,file='./GSE53625/GSE53625_eSet.Rdata')

load(file='./GSE53625/GSE53625_eSet.Rdata')
b = eSet[[1]]
raw_exprSet= exprs(b)
raw_exprSet[1:4,1:4]
phe=pData(b)
save(phe,file = "GSEGSE53624_phe.rdata")
table(phe$geo_accession %in% GSE53625_phe$Id)

GSE53624_phe <- GSE53625_phe[GSE53625_phe$Id %in% phe$geo_accession,]
GSE53624_dat <- GSE53625_dat[,GSE53624_phe$Id]

#合并验证集-----
TG_phe <- rbind(TCGA_phe,GSE53624_phe)
TG_dat <- cbind(TCGA_dat,GSE53624_dat)
save(TG_dat,TG_phe,file = "./TG.rdata")















































