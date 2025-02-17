rm(list=ls())
options(stringsAsFactors = F)
#加载
setwd("/home/datahup/syj/NSCC/step6_1TCGA_GEO/")
load(file = "./TCGA_output.Rdata")
setwd("/home/datahup/syj/ESCC/合并数据集-1.3/")
load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")

#检查数据---
library(sva)
library(limma)
tcga_exp<-dat_2
boxplot(tcga_exp[,1:20],las="2")
tcga_group<-data.frame(row.names = colnames(tcga_exp),
                       Sample=colnames(tcga_exp),
                       DataSet="TCGA")

gtex_exp<-dat53625
boxplot(gtex_exp[,1:20],las="2")
gtex_group<-data.frame(row.names = colnames(gtex_exp),
                       Sample=colnames(gtex_exp),
                       DataSet="GTEx")

dat_group<-rbind(tcga_group,gtex_group)
com_ensg<-intersect(rownames(gtex_exp),rownames(tcga_exp))

dat_exp_before<-cbind(tcga_exp[com_ensg,],
                      gtex_exp[com_ensg,])
#查看两组数据差别-----
boxplot(dat_exp_before[,70:120],las="2")
library(factoextra)
library(FactoMineR)
dat_group$DataSet[dat_group$DataSet == "GTEx"] <- "GSE53624"
before_exp_pca<-PCA(t(dat_exp_before[,rownames(dat_group)]),
                    scale.unit=T,ncp=5,graph=F)
before_exp_pca.plot<-fviz_pca_ind(before_exp_pca,
                                  axes=c(1,2),
                                  label="none",
                                  addEllipses = T,
                                  ellipse.level=0.9,
                                  habillage = factor(dat_group$DataSet),
                                  palette = "aaas",
                                  mean.point=F,
                                  title="")
before_exp_pca.plot#保存pdf为8*6

# before_exp_pca.plot<-before_exp_pca.plot+theme_bw()+
#   theme(legend.direction = "horizontal",legend.position = "top")+
#   xlim(-150,150)+ylim(-150,150)+
#   xlab("Dim1")+ylab("Dim2")
# before_exp_pca.plot

#去除批次效应-----
library(limma)
dat_exp<-removeBatchEffect(dat_exp_before[,rownames(dat_group)],
                           batch = dat_group$DataSet)
#查看
dim(dat_exp)
boxplot(dat_exp[,70:120],las="2")
#dat_exp<-normalizeBetweenArrays(dat_exp)
#boxplot(dat_exp[,70:120],las="2")

afterexp_pca<-PCA(t(dat_exp[,rownames(dat_group)]),
                  scale.unit=T,ncp=5,graph=F)
afterexp_pca.plot<-fviz_pca_ind(afterexp_pca,
                                axes=c(1,2),
                                label="none",
                                addEllipses = T,
                                ellipse.level=0.9,
                                habillage = factor(dat_group$DataSet),
                                palette = "aaas",
                                mean.point=F,
                                title="")
afterexp_pca.plot#保存pdf为8*6
# afterexp_pca.plot<-afterexp_pca.plot+theme_bw()+
#   theme(legend.direction = "horizontal",legend.position = "top")+
#   xlim(-150,150)+ylim(-150,150)+
#   xlab("Dim1")+ylab("Dim2")
# afterexp_pca.plot

#保存-----
merge_dat <- dat_exp
save(merge_dat,file = "合并数据.Rdata")
write.csv(merge_dat, file = "合并数据.csv")





















