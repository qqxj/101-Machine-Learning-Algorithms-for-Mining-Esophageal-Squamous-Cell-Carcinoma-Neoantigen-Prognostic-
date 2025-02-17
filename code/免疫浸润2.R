# rm(list = ls())
 options(stringsAsFactors = F)
 setwd('/home/datahup/syj/ESCC/ssGSEA-9.1/')
# load(file = 'FPKM_ssGSEA.Rdata')
# load(file = 'ssGSEA_output.Rdata')
# 相关性分析
library(ggplot2)
library(ggstatsplot)
library(cowplot)

 load(file = './ssGSEA_output.Rdata')
 
 colnames(ImmSet)
Immcell <- colnames(ImmSet)[2:10]
for (i in 1:length(Immcell)) {
  x=i+1
  a <- shapiro.test(ImmSet[,x])
  ifelse(a[["p.value"]] > 0.05,
         print(paste0(Immcell[i],' is Yes')),
         print(paste0(Immcell[i],' is No')))
}

Immcell <- colnames(ImmSet)[2:10]
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

plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],
          plotlist[[4]],plotlist[[5]],plotlist[[6]],
          nrow = 2)#,plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]]
# p1 <- ggscatterstats(data = ImmSet,
#                      y = RiskScore,
#                      x = 'Activated CD8 T cell',
#                      # centrality.para = "mean",
#                      # margins = "both",
#                      bf.message = FALSE,#去除贝叶斯相关的统计值
#                      type = "nonparamatric",#选择非参数检验
#                      xfill = "#CC79A7",
#                      yfill = "#009E73",
#                      # marginal.type = "density", # 其中marginal.type可选 histograms，boxplots，density，violin，densigram (density + histogram)
#                      title = "Relationship between RiskScore and Activated CD8 T cell")
# p1

# 计算相关性
colnames(ImmSet)
Cor_imm = ImmSet[,c(2:10,16)]
colnames(Cor_imm)
Cor_imm = Cor_imm[,c(10,1:9)]
Cor_imm1 = cor(Cor_imm, method = "pearson")#pearson spearman
head(Cor_imm1)
# 计算显著性
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(Cor_imm1)
library(corrplot)
pdf(file = '/home/datahup/syj/ESCC/ssGSEA-9.1/Imm_risk.pdf',width = 14,height = 12)
corrplot(Cor_imm1, type="upper",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between Risk Scores and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
#p.mat = p.mat, sig.level = 0.05 , insig = "blank")  #添加显著性

dev.off()

#美化-------
pdf(file = './11.14结果/MHC_Imm_risk.pdf',width = 14,height = 12)
corrplot(Cor_imm1, type="upper",
         title = 'Spearman Correlation between Risk Scores and Immune Cells',
         mar=c(0, 0,3, 0),
         tl.col="black", tl.srt=45, tl.cex = 0.9,
         addCoef.col="black", number.cex=0.7,
         method="ellipse",
         diag=F,
         col=colorRampPalette(c("green","#FFFFFF", "#FF0000"))(100))  
dev.off()

"#FF0000", 'green'"#191970", "#FFFFFF", "#C72C41"
# 11gene immune(跳过)--------
#load(file = 'ssGSEA_output.Rdata')
cor_df = ImmSet[,1:40]
rownames(cor_df) = cor_df$id
cor_df = cor_df[,-1]
cor_df = cor_df[,c(29:39,1:28)]
# 计算相关性
cor_imm = cor(cor_df, method = "spearman")
head(cor_imm)
# 计算显著性
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(cor_imm)

# DUSP2
cor_imm1 = cor_imm[-c(2:11),-c(2:11)]
p.mat1 = p.mat[-c(2:11),-c(2:11)]
library(corrplot)
pdf(file = '28Imm_DUSP2.pdf',width = 14,height = 12)
corrplot(cor_imm1, type="lower",
         # addCoef.col = "black", #添加相关系数
         title = 'Spearman Correlation between MS4A7 Expression and Immune Cells',mar=c(0, 0, 1.2, 0),
         tl.col="black", tl.srt=45,tl.cex = 0.9, #修改字体
         addCoef.col="black",number.cex=0.7,
         method="ellipse",
         # addCoefasPercent = T,
         diag = F)
# p.mat = p.mat1, sig.level = 0.05 , insig = "blank")  #添加显著性
dev.off()




























