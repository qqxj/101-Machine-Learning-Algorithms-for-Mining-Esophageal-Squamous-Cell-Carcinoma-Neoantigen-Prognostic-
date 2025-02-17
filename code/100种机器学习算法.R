#100种机器学习算法
#https://blog.csdn.net/zfyyzhys/article/details/141048663#:~:text=Golang%20%E4%BD%9C%E4%B8%BA%E4%B8%80

# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
# if (!requireNamespace("BiocManager", quietly = TRUE)) 
#   install.packages("BiocManager")

# depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
# for(i in 1:length(depens)){
#   depen<-depens[i]
#   if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
# }

# if (!requireNamespace("CoxBoost", quietly = TRUE))
#   devtools::install_github("binderh/CoxBoost")

# if (!requireNamespace("fastAdaboost", quietly = TRUE))
#   devtools::install_github("souravc83/fastAdaboost")
# 
# if (!requireNamespace("Mime", quietly = TRUE))
#   devtools::install_github("l-magnificence/Mime")
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/100种机器学习算法-8.1/")

#数据准备
{
  load(file = "/home/datahup/syj/ESCC/取交集韦恩图-7/交集结果.Rdata")
  #geo训练集数据------
  load(file = "/home/datahup/syj/ESCC/下载geo数据-1.2/dat53625.Rdata")
  load(file = "/home/datahup/syj/ESCC/下载geo数据-1.2/phe53625.Rdata")
  table(a %in% rownames(dat53625))
  datset_geo <- as.data.frame(t(dat53625[a,]))
  datset_geo[1:4,1:4]
  datset_geo$Id <- rownames(datset_geo)
  phe1 <- subset(phe53625,group =="cancer")
  phe1$OS <- ifelse(phe1$OS == "Dead",1,0)
  colnames(phe1)
  phe1 <- phe1[,c(7,8)]
  phe1$Id <- rownames(phe1)
  table(phe1$Id %in% datset_geo$Id)
  dat1 <- datset_geo[rownames(phe1),]
  datset1 <- merge(phe1, datset_geo, by = "Id")
  datset1[1:4,1:4]
save(datset1,file = "GEO训练集.Rdata")

 #验证集------
rm(list=ls())
load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
load(file = "/home/datahup/syj/ESCC/取交集韦恩图-7/交集结果.Rdata")
dat2 <- as.data.frame(t(TG_dat[a,]))
dat2[1:4,1:4]
dat2$Id <- rownames(dat2)
phe2 <- subset(TG_phe,group =="cancer")
rownames(phe2) <- phe2$Id
phe2$OS <- ifelse(phe2$OS == "Dead",1,0)
colnames(phe2)
phe2 <- phe2[,c(6,9)]
phe2$Id <- rownames(phe2)
table(phe2$Id %in% dat2$Id)
dat2 <- dat2[rownames(phe2),]
datset2 <- merge(phe2, dat2, by = "Id")
datset2[1:4,1:4]
save(datset2,file = "TCGA验证集.Rdata")
  
}

library(CoxBoost)
library(fastAdaboost)
library(Mime1)

rm(list=ls())
load(file = "GEO训练集.Rdata")
load(file = "TCGA验证集.Rdata")
load(file = "/home/datahup/syj/ESCC/取交集韦恩图-7/交集结果.Rdata")

genelist <- data.frame(gene = a)
genelist <- as.character(genelist$gene)

library(dplyr)
head(colnames(datset1))
datset1 <- datset1[, c(1, 3, 2, 4:ncol(datset1))]
names(datset1)[1] <- "ID"
datset2 <- datset2[, c(1, 3, 2, 4:ncol(datset2))]
names(datset2)[1] <- "ID"

datset2$ID
is_tcga <- grepl("^TCGA", datset2$ID)
datset2_sorted <- datset2[order(!is_tcga), ]
range(datset2_sorted$OS.time)
range(datset1$OS.time)
datset2_sorted <- subset(datset2_sorted,datset2_sorted$OS.time > 0)
datset1$OS.time <- datset1$OS.time/365
datset2_sorted$OS.time <- datset2_sorted$OS.time/365
# datset2 <- subset(datset2,datset2$OS.time > 0)
# datset2 <- subset(datset2,datset2$OS.time < 4000)

list_train_vali_Data <- list(Training_set = datset1, 
                             Validation_set = datset2_sorted)
# names(list_train_vali_Data) <- c("Training_set",
#                                  "Validation_set")
library(rms)
dd <- datadist(list_train_vali_Data$Training_set)  # 使用 Training_set 定义数据分布对象
options(datadist = "dd")  # 设置全局选项
library(RSpectra)
setwd("/home/datahup/syj/ESCC/100种机器学习算法-8.1/")
# library(snowfall)
# sfInit(parallel = FALSE) 
#跑通这一步就跳过
{
  res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Training_set,
                        list_train_vali_Data = list_train_vali_Data,
                        unicox.filter.for.candi = T,
                        unicox_p_cutoff = 0.05,
                        candidate_genes = genelist,
                        mode = 'all',nodesize =5,seed = 123 )
  
  save(res,file = "100种机器学习.Rdata")
  }

load(file = "100种机器学习.Rdata")

# 设置 PDF 保存路径和文件名
pdf("./cindex_plot.pdf", width = 6, height = 12)  # 设置 PDF 文件名和尺寸
#结果1(10*15)
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order =names(list_train_vali_Data),
               width = 0.5)

# 关闭设备，保存文件
dev.off()


#指定模型(可以自定义)
cindex_dis_select(res,
                  model="StepCox[forward] + RSF",#自定义选择最高的
                  order= names(list_train_vali_Data))


#根据不同数据集中特定模型计算的风险评分绘制患者的生存曲线：
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + RSF",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)

#最后，可以看看指定模型的森林图(X)
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = "StepCox[forward] + RSF",
                                   type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
























































































