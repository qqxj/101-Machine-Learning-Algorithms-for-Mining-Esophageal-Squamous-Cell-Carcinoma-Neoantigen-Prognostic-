library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(cowplot)

rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的功能-9.2/')
#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

Gene <- rownames(gene)
# 检查分组情况
#rownames(dat) <- dat$Id
# a <- dat[dat$type == "Normal",]
# table(substr(rownames(a), 14,15))
# expr1 <- as.data.frame(t(NSCLCcount[Gene,rownames(sv1)]))
# # expr1 <- as.data.frame(t(NSCLCcount[Gene,]))
# table(substr(rownames(expr1), 14,15))
# expr1$Group <- ifelse(substr(rownames(expr1), 14,15) == 11 ,"Normal","Tumor")
# table(expr1$Group)
# Gene

expr1 <- as.data.frame(t(dat53625[rownames(gene),]))
expr1$Group <- phe53625$group
expr1$Group <- ifelse(expr1$Group == "cancer","Tumor","Normal")
rownames(expr1)
#expr1 <- expr1[grep("^GSM129", rownames(expr1)), ]
# 箱图------
{#CSRP1
p1 <- ggplot(data=expr1)+ 
  geom_boxplot(mapping=aes(x=Group,y=CSRP1,colour = Group ), 
               alpha = 0.5,
               size=1.2,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=Group,y=CSRP1,colour = Group), #散点
              alpha = 0.3,size=2)+
  scale_color_manual(limits=c("Tumor","Normal"), 
                     values=c("red","blue"))+ #颜色
  geom_signif(mapping=aes(x=Group,y=CSRP1), # 不同组别的显著性
              comparisons = list(c("Tumor", "Normal")), # 哪些组进行比较
              # c("A", "C"),
              # c("A", "D"),
              # c("B", "C"),
              # c("B", "D"),
              # c("C", "D")),
              map_signif_level=T, # T显示显著性，F显示p value
              # tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              # y_position = c(40,41,42,39,38,40), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 6, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(colour="none")+ # 删除图注
  labs(title='CSRP1',y="Expression(log2(Counts+1))")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 18,
                                  colour = "black",
                                  face = "bold",
                                  vjust = 1.9,
                                  hjust = 0.5),
        axis.title.y = element_text(size = 12, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold.italic", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        axis.title.x = element_blank(),#删除X轴标签
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold.italic", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p1
}

# 分半小提琴图------
library(devtools)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)
{
  a <- rownames(gene)
  gene1 <- as.data.frame(t(dat53625[a[1],]))
  names(gene1) <- "Expression"
  gene1$Gene <- rep(a[1],length(rownames(gene1)))
  #gene1$Type <- ifelse(substr(rownames(gene1), 14,15) == 11 ,"Normal","Tumor")
  gene1$Type <- phe53625$group
  
  gene2 <- as.data.frame(t(dat53625[a[2],]))
  names(gene2) <- "Expression"
  gene2$Gene <- rep(a[2],length(rownames(gene2)))
  #gene2$Type <- ifelse(substr(rownames(gene2), 14,15) == 11 ,"Normal","Tumor")
  gene2$Type <- phe53625$group
  
  
}

a <- rownames(gene)
for (i in 1:5) {
  # 创建临时数据框
  gene_temp <- as.data.frame(t(dat53625[a[i], ]))
  names(gene_temp) <- "Expression"
  gene_temp$Gene <- rep(a[i], nrow(gene_temp))  # 修正为 a[i]
  gene_temp$Type <- phe53625$group
  
  # 使用 assign 创建动态变量名
  assign(paste0("gene", i), gene_temp)
}

df <- rbind(gene1,gene2,gene3,gene4,gene5)#,gene6,gene7,gene8,gene9,gene10,gene11
table(df$Gene)
df$ID <- rownames(df)
rownames(df) <- NULL
#可视化
{
  mypalette
  [1] "#FF0000FF" "#FF9900FF" "#FFCC00FF" "#00FF00FF"
  [5] "#6699FFFF" "#CC33FFFF" "#99991EFF" "#999999FF"
  [9] "#FF00CCFF" "#CC0000FF" "#FFCCCCFF" "#FFFF00FF"
  [13] "#CCFF00FF" "#358000FF" "#0000CCFF" "#99CCFFFF"
  [17] "#00FFFFFF" "#CCFFFFFF" "#9900CCFF" "#CC99FFFF"
}
{
  mypalette = pal_ucscgb()(20)
  ggplot(df,aes(x = Gene, y = Expression, fill = Type)) +
    # split violin
    geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 0,color = 'black',hjust = 0.5),
          legend.position = 'top') + # 图例的位置
    # scale_fill_brewer(palette = 'Set1') +
    # scale_fill_jco(name = '') +
    labs(y="Expression(log2(Count+1))",x=NULL)+ # 添加标题，x轴，y轴内容
    # scale_color_lancet(name = '')+
    #scale_fill_manual(values = mypalette[c(6,1)]) + 
    scale_fill_manual(values = c("#FF0000FF","#99CCFFFF")) + 
    # scale_fill_manual(values = c(pal[1],pal[2]),name = '') +
    # ylim(1,8) + #限定y轴范围
    # 添加显著性标记
    stat_compare_means(aes(group=Type),
                       # symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "NS")),
                       label = "p.signif",
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "NS")),
                       label.y = 22,size = 5)
}
#美化癌症和癌旁可视化--------
ggplot(df, aes(x = Gene, y = Expression, fill = Type)) +
  # 箱线图
  geom_boxplot(alpha = 0.5, color = "black", width = 0.6, outlier.shape = NA, 
               position = position_dodge(0.8)) +
  # 散点图（jittered散点）
  geom_jitter(aes(color = Type), width = 0.2, height = 0, size = 2, alpha = 0.7) +
  # mean point
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.8), size = 3, color = "red") +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = 0.2, 
               size = 0.3, position = position_dodge(0.8)) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 0, color = 'black', hjust = 0.5),
        legend.position = 'top') + # 图例的位置
  labs(y = "Expression(log2(Count+1))", x = NULL) + # 添加标题，x轴，y轴内容
  scale_fill_manual(values = c("#FF0000FF", "#999999FF")) + 
  # 添加显著性标记
  stat_compare_means(aes(group = Type),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "NS")),
                     label.y = 22, size = 5)



#km-----
colnames(sv2)
range(sv2$OS.time)
#sv2 <- subset(sv2,OS.time > "")
dat3 <- sv2[,c("ID","OS","OS.time",rownames(gene))]
dat3$OS.time=dat3$OS.time/365
dat4 = data.frame(OS = dat3$OS, OS.time = dat3$OS.time,row.names = dat3$ID)

# 最佳生存节点
library(survival)
library(survminer)
library(ggsci)
#load(file = "../100种机器学习算法-8.1/geo数据格式.Rdata")
#dat3 <- datset1
# PMEPA1 
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "PMEPA1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  PMEPA1 = ifelse(dat3$PMEPA1<res.cut$cutpoint[,1],"Low","High")
  PMEPA1 = factor(PMEPA1)
  sfit <- survfit(Surv(OS.time, OS)~PMEPA1, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="PMEPA1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for PMEPA1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# MAGEA4  
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "MAGEA4") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  MAGEA4 = ifelse(dat3$MAGEA4<res.cut$cutpoint[,1],"Low","High")
  MAGEA4 = factor(MAGEA4)
  sfit <- survfit(Surv(OS.time, OS)~MAGEA4, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="MAGEA4",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for MAGEA4", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# RCN1  
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "RCN1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  RCN1 = ifelse(dat3$RCN1<res.cut$cutpoint[,1],"Low","High")
  RCN1 = factor(RCN1)
  sfit <- survfit(Surv(OS.time, OS)~RCN1, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="RCN1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for RCN1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# DLX5  
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "DLX5") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  DLX5 = ifelse(dat3$DLX5<res.cut$cutpoint[,1],"Low","High")
  DLX5 = factor(DLX5)
  sfit <- survfit(Surv(OS.time, OS)~DLX5, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="DLX5",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for DLX5", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# TIMP1  
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "TIMP1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  TIMP1 = ifelse(dat3$TIMP1<res.cut$cutpoint[,1],"Low","High")
  TIMP1 = factor(TIMP1)
  sfit <- survfit(Surv(OS.time, OS)~TIMP1, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="TIMP1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for TIMP1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# CSRP1  这些都不是
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "CSRP1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  CSRP1 = ifelse(dat3$CSRP1<res.cut$cutpoint[,1],"Low","High")
  CSRP1 = factor(CSRP1)
  sfit <- survfit(Surv(OS.time, OS)~CSRP1, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="CSRP1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for CSRP1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# CLDN7  这些都不是
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "CLDN7") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  CLDN7 = ifelse(dat3$CLDN7<res.cut$cutpoint[,1],"Low","High")
  CLDN7 = factor(CLDN7)
  sfit <- survfit(Surv(OS.time, OS)~CLDN7, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="CLDN7",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for CLDN7", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# SLC26A2  这些都不是
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "SLC26A2") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  SLC26A2 = ifelse(dat3$SLC26A2<res.cut$cutpoint[,1],"Low","High")
  SLC26A2 = factor(SLC26A2)
  sfit <- survfit(Surv(OS.time, OS)~SLC26A2, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="SLC26A2",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for SLC26A2", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# FOS  没有 这些都不是
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "FOS") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  FOS = ifelse(dat3$FOS<res.cut$cutpoint[,1],"Low","High")
  FOS = factor(FOS)
  sfit <- survfit(Surv(OS.time, OS)~FOS, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'jco', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="FOS",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for FOS", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# EGR1  没有 这些都不是
{
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "EGR1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  EGR1 = ifelse(dat3$EGR1<res.cut$cutpoint[,1],"Low","High")
  EGR1 = factor(EGR1)
  sfit <- survfit(Surv(OS.time, OS)~EGR1, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'jco', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="EGR1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for EGR1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
# GFRA1  没有  这些都不是
{
  #datset1 <- subset(datset1, OS.time > 100)
  res.cut <- surv_cutpoint(datset1, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = "GFRA1") #需要计算的数据列名
  summary(res.cut) #查看数据最佳截断点及统计量
  GFRA1 = ifelse(datset1$GFRA1<10.67593 ,"Low","High")
  #不行GFRA1 = ifelse(datset1$GFRA1<median(datset1$GFRA1),"Low","High")
  GFRA1 = factor(GFRA1)
  dat4 <- datset1[,c(2,3)]
  rownames(dat4) <- datset1$Id
  sfit <- survfit(Surv(OS.time, OS)~GFRA1, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = 'npg', 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title="GFRA1",
                  legend.labs=c("High","Low"),
                  title="Survival Curve for GFRA1", 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  #ggsave(filename = 'MS4A7_KM.pdf', width = 6.5,height = 6, path = '../../Fig/Step4-4/')
}
library(survival)
library(survminer)

# 假设 datset1 包含数据，包含 OS.time、OS 和 GFRA1 列
results <- data.frame(Cutpoint = numeric(), PValue = numeric())

# 计算百分位数
percentiles <- quantile(datset1$GFRA1, probs = seq(0.01, 0.99, by = 0.01), na.rm = TRUE)

# 循环遍历各个百分位数截断点
for (cutpoint in percentiles) {
  
  # 创建分组变量
  GFRA1_group <- ifelse(datset1$GFRA1 < cutpoint, "Low", "High")
  GFRA1_group <- factor(GFRA1_group)
  
  # 准备数据集
  dat4 <- datset1[, c("OS.time", "OS")]
  rownames(dat4) <- datset1$Id
  
  # 进行生存分析
  sfit <- survfit(Surv(OS.time, OS) ~ GFRA1_group, data = dat4)
  
  # 进行 Log-rank 检验
  logrank_test <- survdiff(Surv(OS.time, OS) ~ GFRA1_group, data = dat4)
  p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  
  # 存储结果
  results <- rbind(results, data.frame(Cutpoint = cutpoint, PValue = p_value))
}

# 输出 p 值小于 0.05 的结果
significant_results <- results[results$PValue < 0.05, ]
print(significant_results)


library(ggplot2)

# 创建一个空的列表来保存图形
plots <- list()

# 循环绘制每个显著的截断点
for (i in seq_along(significant_results$Cutpoint)) {
  cutpoint <- significant_results$Cutpoint[i]
  
  # 创建分组变量
  GFRA1_group <- ifelse(datset1$GFRA1 < cutpoint, "Low", "High")
  GFRA1_group <- factor(GFRA1_group)
  
  # 准备数据集
  dat4 <- datset1[, c("OS.time", "OS")]
  rownames(dat4) <- datset1$Id
  
  # 进行生存分析
  sfit <- survfit(Surv(OS.time, OS) ~ GFRA1_group, data = dat4)
  
  # 绘制生存曲线
  p <- ggsurvplot(sfit,
                  palette = 'npg',
                  conf.int = TRUE,
                  pval = TRUE,
                  legend.title = "GFRA1",
                  legend.labs = c("Low", "High"),
                  title = paste("Survival Curve for Cutpoint:", round(cutpoint, 2)),
                  xlab = "Time (Years)",
                  ylab = "Survival Probability",
                  ggtheme = theme_bw(base_size = 12))
  
  # 将图形添加到列表
  plots[[i]] <- p
}

# 打印每个图形
for (plot in plots) {
  print(plot)
}

# 创建一个空的数据框来存储结果
risk_counts <- data.frame(Cutpoint = numeric(), Low_Risk_Count = numeric(), High_Risk_Count = numeric())

# 循环绘制每个显著的截断点
for (i in seq_along(significant_results$Cutpoint)) {
  cutpoint <- significant_results$Cutpoint[i]
  
  # 创建分组变量
  GFRA1_group <- ifelse(datset1$GFRA1 < cutpoint, "Low", "High")
  GFRA1_group <- factor(GFRA1_group)
  
  # 统计低风险组和高风险组的样本数量
  low_risk_count <- sum(GFRA1_group == "Low")
  high_risk_count <- sum(GFRA1_group == "High")
  
  # 将结果存储到数据框中
  risk_counts <- rbind(risk_counts, data.frame(Cutpoint = cutpoint, 
                                               Low_Risk_Count = low_risk_count, 
                                               High_Risk_Count = high_risk_count))
}

# 打印每个截断点的样本数量
print(risk_counts)


#循环(用不了)
for (i in 1:11) {
  res.cut <- surv_cutpoint(dat3, #数据集
                           time = "OS.time", #生存时间
                           event = "OS", #生存状态
                           variables = a[i]) #需要计算的数据列名
  b <- a[i]
  b = ifelse(dat3[,c(i+3)]<res.cut$cutpoint[,1],"Low","High")
  b = factor(b)
  #KM
  sfit <- survfit(Surv(OS.time, OS)~b, data=dat4)
  p1 = ggsurvplot(sfit,
                  palette = c("#FF0000FF","#00FF00FF"), 
                  conf.int = T,conf.int.style='step', 
                  pval = T,pval.method = T,
                  # risk.table = T,risk.table.pos='in',
                  legend=c(0.85,0.85),
                  legend.title=a[i],
                  legend.labs=c("High","Low"),
                  title=paste0("Survival Curve for ", a[i]), 
                  xlab ="Time(Years)",
                  surv.median.line = "hv",
                  ggtheme = theme_bw(base_size = 12))
  p1
  assign(paste0("p", i), p1)
}

p1
p2
p3
p4
p5
p6
p7
p8
p9
p10
p11



















