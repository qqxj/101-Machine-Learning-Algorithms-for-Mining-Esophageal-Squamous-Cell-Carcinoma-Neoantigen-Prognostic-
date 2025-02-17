rm(list = ls())
options(stringsAsFactors = F)
gc()
setwd("/home/datahup/syj/ESCC/通路评分-5.1/")
load(file = "/home/datahup/syj/ESCC/copykat-5/copykat_Epi.Rdata")
WNT <- read.table("KEGG_WNT_SIGNALING_PATHWAY.v2023.2.Hs.tsv", 
                   header = TRUE, 
                   sep = "\t")
notch <- read.table("KEGG_NOTCH_SIGNALING_PATHWAY.v2023.2.Hs.tsv", 
                  header = TRUE, 
                  sep = "\t")
p13k <- read.table("HALLMARK_PI3K_AKT_MTOR_SIGNALING.v2023.2.Hs.tsv", 
                  header = TRUE, 
                  sep = "\t")
#gsva-----------
library(Seurat)
library(AUCell)
library(GSVA)
table(scRNA@active.ident)
Idents(scRNA) <- scRNA@meta.data[["copykat.pred"]]
geneSet <- list( `WNT Score` = c("FRAT1","APC2","NFAT5","WIF1","FZD10","CHP1","CSNK1A1L","CREBBP","PRICKLE1","CSNK1A1","CSNK1E","CSNK2A1","CSNK2A2","CSNK2B","CTBP1","CTBP2","CTNNB1","PRICKLE2","DVL1","DVL2","DVL3","EP300","DKK1","DAAM1","PLCB1","FBXW11","FRAT2","DAAM2","FZD2","CACYBP","DKK4","DKK2","GSK3B","APC","JUN","RHOA","LRP6","LRP5","SMAD2","SMAD3","SMAD4","MMP7","MYC","NFATC1","NFATC2","NFATC3","NFATC4","LEF1","WNT16","NLK","PLCB2","PLCB3","PLCB4","WNT4","PPARD","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B","PPP2R5A","PPP2R5B","PPP2R5C","PPP2R5D","PPP2R5E","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP3R2","PRKACA","PRKACB","PRKACG","PRKCA","PRKCB","PRKCG","MAPK8","MAPK9","MAPK10","PRKX","PSEN1","CTNNBIP1","VANGL2","CHD8","RAC1","RAC2","RAC3","SENP2","CCND1","ROCK1","CHP2","SFRP1","SFRP2","SFRP4","SFRP5","SOX17","SIAH1","PORCN","SKP1","MAP3K7","TBL1X","TCF7","TCF7L2","TP53","SKP1P2","WNT1","WNT2","WNT3","WNT5A","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT10B","WNT11","WNT2B","WNT9A","WNT9B","FZD5","TBL1XR1","FZD3","CXXC4","WNT10A","FOSL1","WNT5B","CAMK2A","CAMK2B","CAMK2D","CAMK2G","VANGL1","AXIN1","AXIN2","FZD1","FZD4","FZD6","FZD7","FZD8","FZD9","TCF7L1","CUL1","NKD1","NKD2","RUVBL1","CCND2","BTRC","CCND3","WNT3A","TBL1Y","CER1","ROCK2","RBX1"),
                `NOTCH Score` = c("DLL3","RBPJL","DTX2","CREBBP","CTBP1","CTBP2","DTX3L","PTCRA","JAG1","DTX1","DVL1","DVL2","DVL3","DTX3","EP300","SNW1","DTX4","NCSTN","KAT2A","DLL1","HDAC1","HDAC2","HES1","RBPJ","JAG2","HES5","LFNG","MFNG","NOTCH1","NOTCH2","NOTCH3","NOTCH4","APH1A","DLL4","MAML3","PSENEN","PSEN1","PSEN2","RFNG","ADAM17","MAML2","NUMB","KAT2B","NUMBL","CIR1","NCOR2","MAML1"),
                `PI3K/AKT/mTOR Score` = c("ACACA","ACTR2","ACTR3","ADCY2","GRK2","AKT1","AKT1S1","AP2M1","ARF1","ARHGDIA","ARPC3","ATF1","CAB39","CAB39L","CALR","CAMK4","CDK1","CDK2","CDK4","CDKN1A","CDKN1B","CFL1","CLTC","CSNK2B","CXCR4","DAPP1","DDIT3","DUSP3","E2F1","ECSIT","EGFR","EIF4E","FASLG","FGF17","FGF22","FGF6","GNA14","GNGT1","GRB2","GSK3B","HRAS","HSP90B1","IL2RG","IL4","IRAK4","ITPR2","LCK","MAP2K3","MAP2K6","MAP3K7","MAPK1","MAPK10","MAPK8","MAPK9","MAPKAP1","MKNK1","MKNK2","MYD88","NCK1","NFKBIB","NGF","NOD1","PAK4","PDK1","PFN1","PIK3R3","PIKFYVE","PIN1","PITX2","PLA2G12A","PLCB1","PLCG1","PPP1CA","PPP2R1B","PRKAA2","PRKAG1","PRKAR2A","PRKCB","PTEN","PTPN11","RAC1","RAF1","RALB","RIPK1","RIT1","RPS6KA1","RPS6KA3","RPTOR","SFN","SLA","SLC2A1","SMAD2","SQSTM1","STAT2","TBK1","THEM4","TIAM1","TNFRSF1A","TRAF2","TRIB3","TSC2","UBE2D3","UBE2N","VAV3","YWHAB")
)

cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA@data, splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(geneSet, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.1)
head(cells_AUC)
#提取打分添加至metadata中
AUCell_auc <- as.data.frame(t(getAUC(cells_AUC)[names(geneSet), ]))
scRNA <- AddMetaData(scRNA, AUCell_auc)
head(scRNA@meta.data)
save(scRNA,
     file="/home/datahup/syj/ESCC/通路评分-5.1/AUCcell.Rdata")

####载入颜色####
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

####-2.8 箱线图 ####
load(file="/home/datahup/syj/ESCC/通路评分-5.1/AUCcell.Rdata")
table(scRNA@meta.data[["copykat.pred"]])
library(dplyr)

scRNA@meta.data <- scRNA@meta.data %>%
  mutate(copykat.pred = recode(copykat.pred, 
                               `cancer cell` = "Cancer cell"))
table(scRNA@meta.data[["copykat.pred"]])

# WNT Score
p <- ggboxplot(scRNA@meta.data, x="copykat.pred", y="WNT.Score", width = 0.6, 
               #color = "black",#轮廓颜色
               fill="copykat.pred",#填充
               palette =c("#98FB98","#6A5ACD"),#分组着色
               xlab = F, #不显示x轴的标签
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=1, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right") #图例
comparisons <- list(c("Epithelial cell", "Cancer cell")) # 根据具体的组名进行调整
p1 <- p + stat_compare_means(comparisons = scRNA@active.ident, method = "t.test", label = "p.signif", p.adjust.method = "bonferroni")
p1



p <- ggboxplot(scRNA@meta.data, x="copykat.pred", y="WNT.Score", width = 0.6, 
               fill="copykat.pred", palette =c("#98FB98","#6A5ACD"),
               xlab = F, bxp.errorbar=T, bxp.errorbar.width=0.5, size=1, 
               outlier.shape=NA, legend = "right")

comparisons <- list(c("Epithelial cell", "cancer cell")) 
p1 <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                             label = "p.signif", p.adjust.method = "bonferroni")
print(p1)


#合并1800*600-------
library(ggpubr)

# WNT Score
# p1 <- ggboxplot(scRNA@meta.data, x="copykat.pred", y="WNT.Score", width = 0.6, 
#                 fill="copykat.pred", palette = c("#40E0D0","#FF6347"),
#                 xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5, size = 1, 
#                 outlier.shape = NA, legend = "none") +
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
#                      label = "p.signif", p.adjust.method = "bonferroni") +
#   ggtitle("WNT Score")
p1 <- ggboxplot(scRNA@meta.data, x = "copykat.pred", y = "WNT.Score", width = 0.6, 
                fill = "copykat.pred", palette = c("#40E0D0","#FF6347"),
                xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5, size = 1, 
                outlier.shape = NA, legend = "none") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", p.adjust.method = "bonferroni") +
  ggtitle("WNT Score") +
  ylab("Score") +  # 设置 y 轴标签
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中

# NOTCH Score
# p2 <- ggboxplot(scRNA@meta.data, x="copykat.pred", y="NOTCH.Score", width = 0.6, 
#                 fill="copykat.pred", palette = c("#40E0D0","#FF6347"),
#                 xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5, size = 1, 
#                 outlier.shape = NA, legend = "none") +
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
#                      label = "p.signif", p.adjust.method = "bonferroni") +
#   ggtitle("NOTCH Score")
p2 <- ggboxplot(scRNA@meta.data, x = "copykat.pred", y = "NOTCH.Score", width = 0.6, 
                fill = "copykat.pred", palette = c("#40E0D0","#FF6347"),
                xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5, size = 1, 
                outlier.shape = NA, legend = "none") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", p.adjust.method = "bonferroni") +
  ggtitle("NOTCH Score") +
  ylab("Score") +  # 设置 y 轴标签
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中

# PI3K.AKT.mTOR Score
# p3 <- ggboxplot(scRNA@meta.data, x="copykat.pred", y="PI3K.AKT.mTOR.Score", width = 0.6, 
#                 fill="copykat.pred", palette = c("#40E0D0","#FF6347"),
#                 xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5, size = 1, 
#                 outlier.shape = NA, legend = "none") +
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
#                      label = "p.signif", p.adjust.method = "bonferroni") +
#   ggtitle("PI3K.AKT.mTOR Score")
p3 <- ggboxplot(scRNA@meta.data, x = "copykat.pred", y = "PI3K.AKT.mTOR.Score", width = 0.6, 
                fill = "copykat.pred", palette = c("#40E0D0","#FF6347"),
                xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5, size = 1, 
                outlier.shape = NA, legend = "none") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", p.adjust.method = "bonferroni") +
  ggtitle("PI3K.AKT.mTOR Score") +
  ylab("Score") +  # 设置 y 轴标签
  theme(plot.title = element_text(hjust = 0.5))  # 设置标题居中
# Combine the plots into one
combined_plot <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1)

# Print the combined plot
print(combined_plot)

#aucell可视化600*600------
load(file="/home/datahup/syj/ESCC/通路评分-5.1/AUCcell.Rdata")

library(ggraph)
library(ggrepel)
library(dplyr)
umap <- data.frame(scRNA@meta.data, scRNA@reductions$umap@cell.embeddings)
head(umap)
cell_type_med <- umap %>%
  group_by(copykat.pred) %>%
  summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))

# WNT.Score NOTCH.Score PI3K.AKT.mTOR.Score
p_uamp1 <- ggplot(umap, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = `WNT.Score`)) + 
  ggrepel::geom_label_repel(aes(label = copykat.pred),fontface="bold",data = cell_type_med,
                            point.padding=unit(0.5, "lines")) +
  scale_color_viridis(option="B")  + #改颜色的参数调option="B",C,F,E,H都还不错
  theme_light(base_size = 15)+
  labs(title = "WNT.Score", colour = "Score")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))

p_uamp1

#NOTCH.Score
p_uamp2 <- ggplot(umap, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = `NOTCH.Score`)) + 
  ggrepel::geom_label_repel(aes(label = copykat.pred),fontface="bold",data = cell_type_med,
                            point.padding=unit(0.5, "lines")) +
  scale_color_viridis(option="B")  + #改颜色的参数调option="B",C,F,E,H都还不错
  theme_light(base_size = 15)+
  labs(title = "NOTCH.Score", colour = "Score")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))

p_uamp2
#PI3K.AKT.mTOR.Score
p_uamp3 <- ggplot(umap, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = `PI3K.AKT.mTOR.Score`)) + 
  ggrepel::geom_label_repel(aes(label = copykat.pred),fontface="bold",data = cell_type_med,
                            point.padding=unit(0.5, "lines")) +
  scale_color_viridis(option="B")  + #改颜色的参数调option="B",C,F,E,H都还不错
  theme_light(base_size = 15)+
  labs(title = "PI3K.AKT.mTOR.Score", colour = "Score")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))

p_uamp3

library(ggplot2)
library(patchwork)
#combined_plot_1 <- ggarrange(p_uamp1, p_uamp2, p_uamp3, ncol = 3, nrow = 1)
combined_plot_1 <- p_uamp1 + p_uamp2+p_uamp3
print(combined_plot_1)
#保存pdf为18*6
ggsave(
  filename = "./combined_plot_1.pdf", 
  plot = combined_plot_1, 
  width = 18, # 图宽
  height = 6, # 图高
  dpi = 300   # 分辨率
)









