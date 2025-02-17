#9.30
#过表达基因和突变基因取交集
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/取交集韦恩图-7/")

#过表达基因
load(file = "../DEG-2.3/DEseq2.Rdata")
Up_gene_deseq2
#突变基因
lusc.snv = read.delim(gzfile('/home/datahup/syj/NSCC/step8_1TCGA_snv/TCGA-ESCA.mutect2_snv.tsv.gz'))
load(file = "../../ESCC/拟时序-6/癌细胞_Statemarker.Rdata")
table(unique(Statemarker$gene) %in% unique(lusc.snv$gene))#用这个
# FALSE  TRUE 
# 342   743 
Mut_genes <- unique(Statemarker$gene)[unique(Statemarker$gene) %in% unique(lusc.snv$gene)]

table(Up_gene_deseq2 %in% Mut_genes)
# FALSE  TRUE 
# 1167    36 
a <- unique(Up_gene_deseq2)[unique(Up_gene_deseq2) %in% unique(Mut_genes)]
save(a,file = "交集结果.Rdata")

library(VennDiagram)
venn_list <- list(Mut_genes = Mut_genes,
                  Over_genes= Up_gene_deseq2)
# venn.diagram(venn_list, filename = 'IntersectGene.png', imagetype = 'png', 
#              fill = c('red', 'green','purple'), alpha = 0.50, 
#              cat.col = c('red', 'green','purple'), cat.cex = 1.2, cat.fontfamily = 'serif',cat.default.pos = "outer",cat.pos = c(-5, 0, 0),  # 位置，用圆的度数
#              cat.dist = c(0.055, 0.055, 0.035),  # 位置，离圆的距离
#              col = c('red', 'green','purple'), cex = 1.5, fontfamily = 'serif')

venn.diagram(venn_list, filename = 'IntersectGene.png', imagetype = 'png', 
             fill = c("cyan", 'red'), alpha = 0.50, #,'purple'
             cat.col = c("cyan", 'red'), #,'purple'
             cat.cex = 1.2, cat.fontfamily = 'serif',
             cat.default.pos = "outer",cat.pos = c(-5, 0),  # 位置，用圆的度数#, 0
             cat.dist = c(0.02, 0.015),  #, 0.035# 位置，离圆的距离
             col = c("cyan", 'red'), cex = 1.5, fontfamily = 'serif')#,'purple'


#install.packages("ggvenn") 
library(ggvenn)
names(venn_list)
names(venn_list) <- c("Mutant genes", "Overexpression of genes")
names(venn_list)
p1 <- ggvenn(venn_list,
             columns = c("Mutant genes", "Overexpression of genes"),
             fill_color = c("cyan", "red")) 
#+ggtitle("Intersection of Genes Set")  # 添加标题    

p1#pdf为6*6









































