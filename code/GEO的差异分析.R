rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/DEG-2.3/")
#load(file = "DEseq2.Rdata")#结果

load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")

library(dplyr)
library(DESeq2)

dat <- dat53625
# dat <- na.omit(dat)
dat <- 2^dat -1 # FPKM转位Count
dat <- round(dat,digits = 0) # 取整

All_dat <- dat

#构建dds矩阵
countData <- All_dat
#condition <- factor(Group$Sample)
condition <- factor(phe53625$group)
condition#注意一定是癌症在后、癌旁和正常组在前#Levels:normal cancer
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

save(resdata,diff_gene_deseq2,Up_gene_deseq2,Down_gene_deseq2,
     file = "DEseq2.Rdata")


res = resdata
rownames(res) = res$Row.names
res = res[,-1]
library(ggplot2)
for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                          'padj' = res$padj,
                          'State' = rep('No', length(res$log2FoldChange)),
                          row.names = rownames(res))
up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 1), which(for_volcano$padj < 0.05))
down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -1), which(for_volcano$padj < 0.05))
for_volcano[up_sig_indices,'State'] <- 'Up'
for_volcano[down_sig_indices,'State'] <- 'Down'
for_volcano$State <- as.factor(for_volcano$State)
for_volcano$padj <- -log10(for_volcano$padj)

this_tile <- paste0('Cutoff for logFC is 2',
                    '\nThe number of Up gene is ',nrow(for_volcano[for_volcano$State =='Up',]) ,
                    '\nThe number of Down gene is ',nrow(for_volcano[for_volcano$State =='Down',]))

p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = State))+
  geom_point(size = I(1))+
  scale_color_manual(values = c('No'='black', 'Up' = 'red', 'Down' = 'blue'))+
  geom_vline(xintercept = c(1, -1), lty=2, size=I(0.4), colour = 'grey11')+
  geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2, size=I(0.1),colour = 'grey11')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'),
        panel.grid = element_blank())+
  labs(x='log2FoldChange', y = '-log10Pvalue')+
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5)) 

p

#美化（选择这个）----
library(ggplot2)
data <- resdata
data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))

this_tile <- paste0('Cutoff for logFC is 1',
                    '\nThe number of Up gene is ',nrow(for_volcano[for_volcano$State =='Down',]) ,
                    '\nThe number of Down gene is ',nrow(for_volcano[for_volcano$State =='Up',]))

p <- ggplot(data,aes(log2FoldChange, -log10(padj)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,1.5))+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5)) 
p#pdf为8*6


#美化2------
data <- resdata
data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))

ggplot(data,aes(log2FoldChange, -log10(padj)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(padj), color= -log10(padj)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  # 调整主题和图例位置：
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.7),
        legend.justification = c(0,1)
  )+
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none")+
  # 添加标签：
  geom_text(aes(label=label, color = -log10(padj)), size = 3, vjust = 1.5, hjust=1)+
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")+
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5)) 












