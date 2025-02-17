
rm(list = ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/预后模型的临床应用-8.3/")

load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

exprSet <- as.data.frame(t(dat53625[gene$sig_gene_multi_cox,]))
exprSet$RiskScore <- exprSet[,rownames(gene)[1]]*gene[1,2]+
  exprSet[,rownames(gene)[2]]*gene[2,2]+
  exprSet[,rownames(gene)[3]]*gene[3,2]+
  exprSet[,rownames(gene)[4]]*gene[4,2]+
  exprSet[,rownames(gene)[5]]*gene[5,2]#+
  # exprSet[,rownames(gene)[6]]*gene[6,2]+
  # exprSet[,rownames(gene)[7]]*gene[7,2]+
  # exprSet[,rownames(gene)[8]]*gene[8,2]+
  # exprSet[,rownames(gene)[9]]*gene[9,2]+
  # exprSet[,rownames(gene)[10]]*gene[10,2]+
  # exprSet[,rownames(gene)[11]]*gene[11,2]
exprSet$RiskGroup <- ifelse(exprSet$RiskScore < median(exprSet$RiskScore) , "Low","High")
table(exprSet$RiskGroup)
rownames(phe53625)# <- phe53625$Id
exprSet <- exprSet[rownames(phe53625),]
phe53625$RiskGroup <- exprSet$RiskGroup
phe53625$RiskScore <- exprSet$RiskScore

phe <- phe53625
phe$OS <- ifelse(phe$OS == "Alive",0 ,1)
table(is.na(phe$OS.time))
library(ggstatsplot)
library(tidyverse) 
library(ggpubr)
datSet <- phe

datSet$group <- ifelse(datSet$group == "normal", 'Normal','Tumor')
table(datSet$gender)
datSet$gender <- ifelse(datSet$gender == "female" , "Female", "Male")
datSet$Age <- ifelse(datSet$Age <= 60 , '<=60','>60')
table(datSet$Age)

save(datSet,file = "./datSet.Rdata")
# 临床信息 单因素Cox----------
datSet0 <- datSet
table(datSet0$gender)
datSet0$gender <- ifelse(datSet0$gender == 'Female',1,2) #女性为1，男性为2
table(datSet0$gender)

table(datSet0$Age)
datSet0$Age <- ifelse(datSet0$Age == '<=60',1,2) # <60为1， >60为2
table(datSet0$Age)

table(datSet0$n_stage)
datSet0$n_stage <- ifelse(datSet0$n_stage == "N0",1,
                            ifelse(datSet0$n_stage == "N1",2,
                                   ifelse(datSet0$n_stage == "N2",3,4)))
table(datSet0$n_stage)

table(datSet0$t_stage)
#datSet0 <- subset(datSet0,datSet0$t_stage != 'TX')
datSet0$t_stage <- ifelse(datSet0$t_stage == "T1",1,
                            ifelse(datSet0$t_stage== "T2",2,
                                   ifelse(datSet0$t_stage == "T3",3,
                                          ifelse(datSet0$t_stage == "T4",4,5))))
table(datSet0$t_stage)

table(datSet0$tnm_stage)
#datSet0 <- subset(datSet0,datSet0$Stage != 'not reported')
datSet0$tnm_stage <- ifelse(datSet0$tnm_stage == "I",1,
                        ifelse(datSet0$tnm_stage == "II",2,
                               ifelse(datSet0$tnm_stage == "III",3,4)))
table(datSet0$tnm_stage)

table(datSet0$RiskGroup)
datSet0$RiskGroup <- ifelse(datSet0$RiskGroup == 'Low',1,2)
table(datSet0$RiskGroup)

colnames(datSet0)
colnames(datSet0)[3:6] <- c("Gender","TStage","NStage","MStage")
covariates <- c("Age","Gender","TStage","NStage","MStage","RiskGroup")

datSet0 <- subset(datSet0,group == "Tumor")

library(survival)
library(survminer)
#分别对每一个变量，构建生存分析的公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time,OS)~', x)))
univ_formulas
#对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = datSet0)})
univ_models
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)",           "wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
str(univ_results)
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)

rownames(res)
res$Names <- c("Age","Gender","T Stage","N Stage","M Stage","Risk Score")
colnames(res)
res<-res[, c("Names","HR (95% CI for HR)","wald.test","p.value","beta")]

#############################################################
#对HR (95% CI for HR)做处理，得到HR和low .95和high .95
#当然也可以改计算univ_results这一步的代码，不要将HR和CI贴起来
############################################################
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)

HR1 <- HR
HR1 <- as.data.frame(t(HR1))
res$HR <- as.numeric(HR1$V1)
res$upper <- as.numeric(HR1$V3)
res$lower <- as.numeric(HR1$V2)
res$p.value <- as.numeric(res$p.value)
res$beta <- as.numeric(res$beta)
rownames(res) <- res$Names
str(res)

res$'P value' = ifelse(
  res$p.value < 0.001,
  "<0.001 ***",
  ifelse(
    res$p.value < 0.01,
    "<0.01  **",
    ifelse(
      res$p.value < 0.05,
      paste(round(res$p.value, 2), " *"),
      round(res$p.value, 2)
    )
  )
)

# 作图1(自定义画法)-----
# #左边和右边边距稍微留多一点来写变量名称，pvalue和HR
# par(mar=c(5,6,4,13))
# #先用小方块画出HR
# plot(as.numeric(HR[1,]),1:dim(HR)[2],
#      pch=15,cex=2,col="blue",bty='n',yaxt='n',ylab=NA,xlab="Hazard Ratio",
#      xlim=range(as.numeric(unlist(HR)))
# )
# #添加中线
# abline(v=1,col="blue",lwd=2,lty=2)
# 
# for(i in 1:ncol(HR)){
#   x=as.numeric(HR[2:3,i])
#   #循环画出CI
#   lines(x,c(i,i),col="blue")
#   #添加变量名
#   text(0.2,i,rownames(res)[i],xpd=T,adj = c(0,0))
#   #添加HR和CI
#   text(1.9,i,as.character(res[i,2]),xpd=T,adj = c(0,0))
#   #添加p值
#   text(2.8,i,as.numeric(res[i,9]),xpd=T,adj = c(0,0))
# }
# #添加标题
# text(1.9,ncol(HR)+1,"HR(95% CI)",xpd=T,adj = c(0,0),font=2)
# text(2.8,ncol(HR)+1,"P value",xpd=T,adj = c(0,0),font=2)
# text(0.2,ncol(HR)+1,"Names",xpd=T,adj = c(0,0),font=2)
# text(0.8,ncol(HR)+2,"Univariate Cox Regression Analysis",xpd=T,adj = c(0,0),font = 4)
# 
# lines(x = c(0.2,2.8),y=c(7,7),col="blue")

# 作图2--------
tabletext <- cbind(c("Variable",res$Names),
                   c("HR(95% CI)",res$`HR (95% CI for HR)`),
                   c("P Value",res$`P value`))
str(tabletext)
class(tabletext)
tabletext0 <- data.frame(mean = res$HR,lower =res$lower, upper = res$upper,row.names = rownames(res))


# jpeg(file = "results_Value_2.jpg",width =2000,height = 1800,units = "px",res =300) #结果保存
library(forestplot)
forestplot(tabletext,  #显示的文本
           mean=c(NA,tabletext0$mean),
           lower=c(NA,tabletext0$lower), upper=c(NA,tabletext0$upper),
           # c(NA,tabletext0$mean), #误差条的均值(此处为差值的中值)
           # c(NA,tabletext0$lower), #误差条的下界(此处为差值的25%分位数)
           # c(NA,tabletext0$upper), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           graph.pos=2,
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.3, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 3,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.8), cex = 1), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Hazard Ratio(HR)",#x轴的标题
           #xticks = F,
           is.summary = c(T, rep(F, 7)),
           align = "l",
           hrzl_lines = list(
             "1" = gpar(lty=1),
             "2" = gpar(lty=1),
             "8"= gpar(lty=1)),
           colgap = unit(6, 'mm'),
           title="Univariate Cox Regression Analysis",
)
# dev.off()
"#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f","#b2db87","#7ee7bb","#64cccf","#a9dce6","#a48cbe","#e4b7d6"
#美化单因素cox------
#保存1000*800
forestplot(tabletext,  
           mean=c(NA, tabletext0$mean),
           lower=c(NA, tabletext0$lower), 
           upper=c(NA, tabletext0$upper),
           zero = 1, 
           xlog=FALSE, 
           graph.pos=2,
           fn.ci_norm = fpDrawCircleCI, 
           boxsize = 0.3, 
           col=fpColors(line = "#ea5c6f", 
                        box="#64cccf"), 
           lty.ci = 7,   
           lwd.ci = 3,   
           ci.vertices.height = 0.15, 
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.8), cex = 1), 
           lineheight = "auto", 
           xlab="Hazard Ratio(HR)",
           is.summary = c(TRUE, rep(FALSE, 7)),
           align = "l",
           hrzl_lines = list(
             "1" = gpar(lty=1),
             "2" = gpar(lty=1),
             "3" = gpar(lty=1),  # 添加有效编号
             "4" = gpar(lty=1),  # 添加有效编号
             "5" = gpar(lty=1),  # 添加有效编号
             "6" = gpar(lty=1),  # 添加有效编号
             "7" = gpar(lty=1),  # 添加有效编号
             "8" = gpar(lty=1)),  # 添加有效编号
           colgap = unit(6, 'mm'),
           title="Univariate Cox Regression Analysis")


# 多因素Cox----------
library(glmnet)
library(survival)
library(survminer)

formula_for_multivarirate <- as.formula(paste0('Surv(OS.time, OS)~',paste(covariates,sep = '',collapse = "+")))
multi_varirate_cox <- coxph(formula_for_multivarirate, data = datSet0)
ph_hypo_multi <- cox.zph(multi_varirate_cox)
ph_hypo_table <-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),] ######
#formula_for_multivarirate <- as.formula(paste0('Surv(os.time, os)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = " ",collapse = "+")))
#multi_varirate_cox <- coxph(formula_for_multivarirate,data = sv1)
multiCoxSum <- summary(multi_varirate_cox)
multi_cox<- as.data.frame(multi_varirate_cox$coefficients)
multi_cox

out_multi <- cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

out_multi
class(out_multi)
out_multi <- as.data.frame(out_multi)

#HR和它的置信区间
dat2 = as.data.frame(round(multiCoxSum$conf.int[, c(1, 3, 4)], 2))
dat2 = tibble::rownames_to_column(dat2, var = "Variable")
colnames(dat2)[2:4] = c("HR", "lower", "upper")
#需要在图上显示的HR文字和p值
dat2$'HR(95% CI)' = paste0(dat2[, 2], "(", dat2[, 3], "-", dat2[, 4], ")")
dat2$'P Value' = out_multi$pvalue

dat2$`P Value` = ifelse(
  dat2$`P Value` < 0.001,
  "<0.001 ***",
  ifelse(
    dat2$`P Value` < 0.01,
    "<0.01  **",
    ifelse(
      dat2$`P Value` < 0.05,
      paste(round(dat2$`P Value`, 2), " *"),
      round(dat2$`P Value`, 2)
    )
  )
)
dat2$Variable <- c("Gender","Age","N Stage","T Stage","M Stage","Risk Score")
str(dat2)

dat3 <- as.data.frame(t(dat2))
dat3$V8 <- rownames(dat3)
dat3 <- dat3[,c(7,1:6)]
dat3 <- as.data.frame(t(dat3))
rownames(dat3) <- NULL
colnames(dat3) <- NULL
dat4 <- dat3[,c(1,5,6)]


#可视化--------
library(forestplot)
forestplot(dat4,  #显示的文本
           mean=c(NA,dat2$HR),
           lower=c(NA,dat2$lower), upper=c(NA,dat2$upper),
           # c(NA,tabletext0$mean), #误差条的均值(此处为差值的中值)
           # c(NA,tabletext0$lower), #误差条的下界(此处为差值的25%分位数)
           # c(NA,tabletext0$upper), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           graph.pos=2,
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.3, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 3,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.8), cex = 1), #文本大小设置
           lineheight = "auto", #线的高度
           xlab="Hazard Ratio(HR)",#x轴的标题
           #xticks = F,
           is.summary = c(T, rep(F, 7)),
           align = "l",
           hrzl_lines = list(
             "1" = gpar(lty=1),
             "2" = gpar(lty=1),
             "8"= gpar(lty=1)),
           colgap = unit(6, 'mm'),
           title="Multivariate Cox Regression Analysis"
)
 

forestplot(dat4,  
           mean=c(NA, dat2$HR),
           lower=c(NA, dat2$lower), 
           upper=c(NA, dat2$upper),
           zero = 1, 
           xlog=FALSE, 
           graph.pos=2,
           fn.ci_norm = fpDrawCircleCI, 
           boxsize = 0.3, 
           col=fpColors(line = "#ea5c6f", 
                        box="#64cccf"), 
           lty.ci = 7,   
           lwd.ci = 3,   
           ci.vertices.height = 0.15, 
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.8), cex = 1), 
           lineheight = "auto", 
           xlab="Hazard Ratio(HR)",
           is.summary = c(TRUE, rep(FALSE, 7)),  # 确保长度一致
           align = "l",
           hrzl_lines = list(
             "1" = gpar(lty=1),
             "2" = gpar(lty=1),
             "3" = gpar(lty=1),  # 添加有效编号
             "4" = gpar(lty=1),  # 添加有效编号
             "5" = gpar(lty=1),  # 添加有效编号
             "6" = gpar(lty=1),  # 添加有效编号
             "7" = gpar(lty=1),  # 添加有效编号
             "8" = gpar(lty=1)),  # 添加有效编号
           colgap = unit(6, 'mm'),
           title="Multivariate Cox Regression Analysis"
)










