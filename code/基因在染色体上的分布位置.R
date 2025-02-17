#9.30
#基因在染色体上的分布位置
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/DEG-2.3/")

load(file = "../ESCC/DEG-2.3/DEseq2.Rdata")

library(karyoploteR)
library(biomaRt)#获取基因再染色体位置信息
library(AnnotationHub)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#gene_names <- c("BRCA1", "TP53", "EGFR")  # 加入想要查的基因
gene_names <- Up_gene_deseq2
gene_locations <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                        filters = "external_gene_name",
                        values = gene_names,
                        mart = ensembl)
save(gene_locations,file = "基因在染色体的位置信息.Rdata")

#创建GRanges格式
 #所有基因的GRanges格式
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# all.genes <- genes(txdb)
# head(all.genes)
table(gene_locations$chromosome_name)
gene_locations <- gene_locations[gene_locations$chromosome_name %in% c(1:22, "X", "Y"), ]
table(duplicated(gene_locations$external_gene_name))
library(dplyr)
gene_locations %>% distinct(external_gene_name, .keep_all = TRUE)
gene_locations$chromosome_name <- paste0("chr", gene_locations$chromosome_name)

gene_locations_gr <- GRanges(seqnames = Rle(gene_locations$chromosome_name),
                             ranges = IRanges(start = gene_locations$start_position, 
                                              end = gene_locations$end_position))
kp <- plotKaryotype(genome="hg19")
#kp <- kpPlotDensity(kp, all.genes)
kp <- kpPlotDensity(kp, gene_locations_gr)

#美化
kp <- plotKaryotype()
kpPlotCoverage(kp, data=gene_locations_gr)
kpPlotCoverage(kp, data=gene_locations_gr, col="red")  












































