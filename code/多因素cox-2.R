#多因素cox
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/机器学习和验证然后取跑100种机器学习-8.2/")

# load(file = "../100种机器学习算法-8.1/geo数据格式.Rdata")
# load(file = "../100种机器学习算法-8.1/验证集数据格式.Rdata")
# load(file = "/home/datahup/syj/ESCC/取交集韦恩图-7/交集结果.Rdata")
# load(file = "lasso结果.Rdata")

# 多因素Cox
{
  library(glmnet)
  library(survival)
  library(survminer)
  covariates3 <- rownames(gene)
  formula_for_multivarirate <- as.formula(paste0('Surv(OS.time, OS)~',paste(covariates3,sep = '',collapse = "+")))
  multi_varirate_cox <- coxph(formula_for_multivarirate, data = sv2)
  ph_hypo_multi <- cox.zph(multi_varirate_cox)
  ph_hypo_table <-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),] ######
  
  multiCoxSum <- summary(multi_varirate_cox)
  multi_cox<- as.data.frame(multi_varirate_cox$coefficients)
  multi_cox
  
  correlation <- cor(sv2[,rownames(ph_hypo_table)],method = 'pearson')
  
  library('GGally')
  ggpairs(sv2[,rownames(ph_hypo_table)],
          axisLabels = "show")+
    theme_bw()+
    theme(panel.background = element_rect(colour = 'red',size = 1,fill = "white"),
          panel.grid = element_blank())
  
  
  library("rms")
  vif <- rms::vif(multi_varirate_cox)
  sqrt(vif) < 2
  library(survival)
  library(survminer)
  ggforest(model = multi_varirate_cox,data = sv2, main =  "Hazard",fontsize = 1)
  
  
  C_index <- multi_varirate_cox$concordance["concordance"]
  if(C_index >= 0.9){print("high accuracy")
  }else{
    if(C_index <0.9 & C_index >= 0.7){
      print("Medium accuracy")
    }else{print('low accuracy')
    }
  }
  
  out_multi <- cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  
  out_multi
  class(out_multi)
  out_multi <- as.data.frame(out_multi)
  gene_multi <- subset(out_multi,out_multi$pvalue < 0.05)
  
  gene_multi
  rownames(gene_multi)
  
}

sv2$exp <- sv2[,rownames(gene_multi)[1]]*gene_multi[1,1]+
  sv2[,rownames(gene_multi)[2]]*gene_multi[2,1]+
  sv2[,rownames(gene_multi)[3]]*gene_multi[3,1]#+
  # sv2[,rownames(gene_multi)[4]]*gene_multi[4,1]+
  # sv2[,rownames(gene_multi)[5]]*gene_multi[5,1]+
  # sv2[,rownames(gene_multi)[6]]*gene_multi[6,1]
save(sv2,gene_multi,file = "Train_Result.Rdata")

dat1<- sv2
# ROC 1,3,5
{library(survivalROC)
  cutoff <- 1*365
  ROC1= survivalROC(Stime=dat1$OS.time,##生存时间
                    status=dat1$OS,## 终止事件    
                    marker =dat1$exp, ## marker value    
                    predict.time = cutoff## 预测时间截点
                    ,method="KM")##span,NNE法的namda
  str(ROC1)## list结构
  
  cutoff <- 3*365
  ROC3= survivalROC(Stime=dat1$OS.time,##生存时间
                    status=dat1$OS,## 终止事件    
                    marker =dat1$exp, ## marker value    
                    predict.time = cutoff## 预测时间截点
                    ,method="KM")##span,NNE法的namda
  str(ROC3)## list结构
  
  cutoff <- 5*365
  ROC5= survivalROC(Stime=dat1$OS.time,##生存时间
                    status=dat1$OS,## 终止事件    
                    marker = dat1$exp, ## marker value    
                    predict.time = cutoff## 预测时间截点
                    ,method="KM")##span,NNE法的namda
  str(ROC5)## list结构
  
  
  plot(ROC5$FP, ROC5$TP, ## x=FP,y=TP
       type="l",col="blue",lwd=2,##线条设置
       xlim=c(0,1), ylim=c(0,1),   
       xlab=paste( "False positive rate"), ##连接
       ylab="True positive rate",
       main="ROC of Training Set"
  )## \n换行符
  # lines(ROC12$FP,ROC12$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC10$FP,ROC10$TP, type="l",col="yellowgreen",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC9$FP,ROC9$TP, type="l",col="pink",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC7$FP,ROC7$TP, type="l",col="red",xlim=c(0,1), ylim=c(0,1),lwd=2)
  # lines(ROC5$FP,ROC5$TP, type="l",col="salmon",xlim=c(0,1), ylim=c(0,1),lwd=2)
  lines(ROC3$FP,ROC3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1),lwd=2)
  lines(ROC1$FP,ROC1$TP, type="l",col="red",xlim=c(0,1), ylim=c(0,1),lwd=2)
  
  legend(0.5,0.4,c(paste("AUC of 1 years =",round(ROC1$AUC,3)),
                   paste("AUC of 3 years =",round(ROC3$AUC,3)),
                   paste("AUC of 5 years =",round(ROC5$AUC,3))),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,col=c("red","green","blue"),
         bty = "n",# bty框的类型
         seg.len=1,cex=0.8)# 
  abline(0,1,col="black",lty=1,lwd=2)##线条颜色
  box(lwd=2)
}



















































