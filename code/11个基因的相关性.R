library(corrplot)
library(ggplot2)
library(ggpubr)
#基因与基因之间相关性-----------
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的功能-9.2/')

load(file = "../100种机器学习算法-8.1/GEO训练集.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

sv2 <- datset1[,c("Id","OS","OS.time",gene$sig_gene_multi_cox)]
td <- datset1[,gene$sig_gene_multi_cox]
cor (td, method="pearson")

tdc <- cor (td, method="pearson")#"pearson", "kendall",
corrplot(tdc)
write.csv(tdc, file = "correlation_matrix.csv", row.names = TRUE)

#添加p值
library(corrplot)
testRes = cor.mtest(td, method="pearson",conf.level = 0.95)
my_colors <- colorRampPalette(colors = c("blue", "white", "red"))(100)
corrplot(tdc, method = "color", col = my_colors, 
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",
         p.mat = testRes$p, diag = T, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')
corrplot(tdc, method = "number", type = "lower",col = my_colors, 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE',
         add = T)

#饼图
addcol <- colorRampPalette(c("blue","white", "red"))
corrplot(tdc, method = "pie", type = "upper",col = addcol(100), 
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,
         tl.pos = "lt")
corrplot(tdc, method = "number", type = "lower",col = addcol(100), 
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",
         add = T)

#直接绘制
library(PerformanceAnalytics)
chart.Correlation(td)

#环形相关性热图--------
library(psych)
library(igraph)
library(tidygraph)
library(ggraph)
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的功能-9.2/')
#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

dat53625[1:4,1:4]
dat53625 <- dat53625[rownames(gene),]
phe <- subset(phe53625,group == "cancer")
dat53625 <- dat53625[,rownames(phe)]

data <- dat53625
pvalue <- 0.05
rvalue <- 0.6#rvalue表示相关性，绝对值大于0.3才有相关性
#pearson相关性
corr_data1 <- corr.test(t(data),method = "spearman",adjust = "fdr")
p <- corr_data1$p
p[p > 0.05] <- -1
p[p <= 0.05 & p>=0] <- 1
p[p == -1] <- 0
corr_data1$r[abs(corr_data1$r) < 0.2] <- 0
corr_data1$r <- corr_data1$r * p
s <- corr_data1$r
m <- graph.adjacency(s,weighted = TRUE,mode = "undirected",diag =FALSE)
m <- simplify(m)
m <- delete.vertices(m,names(degree(m)[degree(m) == 0]))
E(m)$correlation <- E(m)$weight
E(m)$weight <- abs(E(m)$weight)
E(m)$Cor <- ifelse(E(m)$correlation > 0, "Positive","Negative")

taxonomy <- data.frame(name = rownames(dat53625),
                   group = rownames(dat53625))

rownames(taxonomy) <- taxonomy$name
taxonomy <- data.frame(taxonomy[as.character(V(m)$name),])
library(psych)
library(igraph)
library(tidygraph)
library(ggraph)
taxonomy$group <- 1
taxonomy$group <- as.factor(taxonomy$group)
V(m)$Group <- taxonomy$group
V(m)$Name <- taxonomy$name
V(m)$Degree <- degree(m) 

library(ggraph)
p1 <- ggraph(m,layout = "linear",circular = TRUE)+
  geom_edge_arc(aes(width = weight,color = Cor))+
  scale_edge_color_manual(values = c("#F8766D","#00BFC4"))+#"#A9A9A9","#808080","#00BFC4"
  scale_edge_width_continuous(range = c(0.5,2))+
  geom_node_point(aes(fill = Group,size = Degree),shape=21,stroke=0)+
  scale_fill_manual(values = c("#4A90BD","#4CA85F"), guide = "none")+
  scale_size_continuous(range = c(5,8))+
  geom_node_text(aes(x = x*1.1,
                     y = y*1.1,
                     label = Name,
                     angle =-((-node_angle(x,y)+90) %% 180)+90),
                 size = 5,
                 fontface = 2,
                 family = "Arial",
                 hjust = "outward")+
  coord_cartesian(xlim = c(-1.5,1.5),ylim = c(-1.5,1.5))+
  theme_graph()+
  theme(text = element_text(size = 15,face = 2,family = "Arial"))+
  guides(size = "none")

p1                































