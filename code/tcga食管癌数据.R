
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/NSCC/step6_1TCGA_GEO/")

load(file = "./TCGA_output.Rdata")

library(dplyr)
library(DESeq2)

dat <- dat_2
# dat <- na.omit(dat)
dat <- 2^dat -1 # FPKM转位Count
dat <- round(dat,digits = 0) # 取整
All_dat <- dat

#构建dds矩阵
countData <- All_dat
#condition <- factor(Group$Sample)
phe6$group <- ifelse(phe6$SubType == "normal","normal","cancer")
condition <- factor(phe6$group)
head(condition)
countData[1:4,1:4]

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
head(dds)
dim(dds)
dds <- DESeq(dds)#时间较长
resultsNames(dds)
res <- results(dds)
summary(res)

DEG_INFO <- as.data.frame(res)

# 提取差异分析结果
# 获取padj（p值经过多重校验校正后的值）小于0.01，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
table(res$padj<0.05) #取P值小于0.01的结果
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)# 所有差异基因

Up_gene_deseq2 <- subset(res,padj < 0.05 & (log2FoldChange > 1)) 
Up_gene_deseq2 <- row.names(Up_gene_deseq2) # 上调差异基因

Down_gene_deseq2 <-  subset(res,padj < 0.05 & (log2FoldChange <  -1))
Down_gene_deseq2 <- row.names(Down_gene_deseq2) # 下调差异基因

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)

setwd("/home/datahup/syj/ESCC/下载tcga数据-1.1/")
save(resdata,diff_gene_deseq2,Up_gene_deseq2,Down_gene_deseq2,
     file = "DEseq2.Rdata")

load(file = "DEseq2.Rdata")

#查看交集基因集
setwd("/home/datahup/syj/ESCC/取交集韦恩图-7/")
#突变基因
lusc.snv = read.delim(gzfile('/home/datahup/syj/NSCC/step8_1TCGA_snv/TCGA-ESCA.mutect2_snv.tsv.gz'))
load(file = "../../ESCC/拟时序-6/癌细胞_Statemarker.Rdata")
table(unique(Statemarker$gene) %in% unique(lusc.snv$gene))#用这个
# FALSE  TRUE 
# 342   743 
Mut_genes <- unique(Statemarker$gene)[unique(Statemarker$gene) %in% unique(lusc.snv$gene)]

table(Up_gene_deseq2 %in% Mut_genes)


