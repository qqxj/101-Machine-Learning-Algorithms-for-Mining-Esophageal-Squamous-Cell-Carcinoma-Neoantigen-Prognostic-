#可视化
#训练集auc-----
sv2 = sv1
sv2$exp <- sv2[,gene[1,1]]*gene[1,2]+
  sv2[,gene[2,1]]*gene[2,2]+
  sv2[,gene[3,1]]*gene[3,2]+
  sv2[,gene[4,1]]*gene[4,2]+
  sv2[,gene[5,1]]*gene[5,2]+
  sv2[,gene[6,1]]*gene[6,2]+
  sv2[,gene[7,1]]*gene[7,2]+
  sv2[,gene[8,1]]*gene[8,2]+
  sv2[,gene[9,1]]*gene[9,2]+
  sv2[,gene[10,1]]*gene[10,2]+
  sv2[,gene[11,1]]*gene[11,2]
dat1<- sv2
# ROC曲线图(受试者工作特征曲线也叫感受性曲线) 1,3,5
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
       type="l",col="blue",lwd=4,##线条设置
       xlim=c(0,1), ylim=c(0,1),   
       xlab=paste( "1-Sensitivity"), ##连接
       ylab="Sensitivity",
       main="Training set AUC")## \n换行符
  #lines(ROC5$FP,ROC5$TP, type="l",col="salmon",xlim=c(0,1), ylim=c(0,1),lwd=4)
  lines(ROC3$FP,ROC3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1),lwd=4)
  lines(ROC1$FP,ROC1$TP, type="l",col="red",xlim=c(0,1), ylim=c(0,1),lwd=4)
  
  legend(0.6,0.4,#图例放置位置
         c(paste("AUC of 1 =",round(ROC1$AUC,3)),
           paste("AUC of 3 =",round(ROC3$AUC,3)),
           paste("AUC of 5 =",round(ROC5$AUC,3))#,
         ),#paste("AUC of 7 =",round(ROC7$AUC,3))
         x.intersp=1, y.intersp=0.8,#控制文本行之间和列之间的间距
         lty= 1 ,lwd= 4,col=c("red","green","blue"),#线的类型、宽度和颜色
         bty = "n",# bty框的类型
         seg.len=1,cex=0.8)# cex控制文本的大小
  abline(0,1,col="black",lty=2,lwd=3)##线条颜色
  box(lwd=2)
}


#训练集km-----
##### 生存分析
load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
table(TG_phe$Id %in% dat1$Id)
phe <- TG_phe[TG_phe$Id %in% dat1$Id,]
{
  # factGene <- c("CYGB","CLVS1","CBY3","CPNE6")
  
  fac_mix <- dat1
  library(survival)
  library(survminer)
  library(ggplotify)
  library(cowplot)
  library(Hmisc)
  library(pheatmap)
  library(gridExtra)
  #s = as.formula(paste('Surv(OS.time, OS)~', noquote(paste(factGene,collapse = ' + '))))
  #model <- coxph(s, data = fac_mix )
  #summary(model,data=fac_mix)
  # RiskScore <- predict(model,type = "risk")
  RiskScore <- fac_mix$exp
  #names(RiskScore) = rownames(fac_mix)
  fp <- RiskScore
  phe <- fac_mix
  phe$id <- phe$Id
  
  # 生存图 
  {
    #最佳节点
    # res.cut <- surv_cutpoint(dat1, #数据集
    #                          time = "OS.time", #生存时间
    #                          event = "OS", #生存状态
    #                          variables = "exp") #需要计算的数据列名
    # 
    # summary(res.cut) #查看数据最佳截断点及统计量
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    dat1$OS.time=dat1$OS.time/365
    sfit <- survfit(Surv(OS.time, OS)~RiskGroup, data=dat1)
    # ggsurvplot(sfit, pval=TRUE,xlab ="Time(Years)",surv.median.line = "hv",pval.method = T)
    # ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
    #            risk.table =TRUE,pval =TRUE,
    #            conf.int =TRUE,xlab ="Time in years",
    #            ggtheme =theme_light(),
    #            ncensor.plot = TRUE)
    ggsurvplot(sfit,
               palette = c("#FF5733", "#33FF57"),
               conf.int = T,conf.int.style='step', 
               pval = T,pval.method = T,
                risk.table = T,risk.table.pos='in',
               legend=c(0.85,0.85),
               legend.title="Risk Group",
               legend.labs=c("High","Low"),
               title="Survival Curve for Training Set", 
               xlab ="Time(Years)",
               surv.median.line = "hv",
               ggtheme = theme_bw(base_size = 20))
  }
  
  # 散点+热图 中位数
  {
    fp_dat=data.frame(patientid=phe$id,fp=phe$exp)
    fp_dat$RiskGroup= ifelse(fp_dat$fp>= median(fp_dat$fp),'High','Low')
    library(dplyr)
    fp_dat=arrange(fp_dat,fp)
    
    sur_dat=data.frame(patientid=phe$id,time=phe$OS.time/365,Status=phe$OS)
    sur_dat$Status=ifelse(sur_dat$Status==0,'Alive','Dead')
    sur_dat$Status=factor(sur_dat$Status,levels = c("Dead","Alive"))
    sur_dat$time <- sur_dat$time
    rownames(sur_dat)=sur_dat$patientid
    sur_dat=sur_dat[fp_dat$patientid,]

    exp_dat=dat1[,rownames(gene)]
    rownames(exp_dat)=phe$id
    exp_dat=exp_dat[fp_dat$patientid,]
    fp_dat$patientid=1:length(fp)
    sur_dat$patientid=1:length(fp)
    rownames(exp_dat)=1:length(fp)
    rownames(sur_dat)=1:length(fp)
    
    # RiskGroup = ifelse(fp<res.cut$cutpoint[,1],"Low","High")
    RiskGroup = ifelse(fp<median(fp_dat$fp),"Low","High")
    RiskGroup = factor(RiskGroup)
    
    
    # p1 <- ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=RiskGroup))+
    #   scale_colour_manual(values = c("#FF6666","#00CCCC"))+
    #   theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
    #   geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
    #   # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
    #   geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
    #              colour="black", linetype="dotted",size=0.8) +
    #   theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
    #         axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
    #         legend.text=element_text(size=14),legend.title =element_text(size=14))
    # p1
    
    library(ggsci)
    library(scales)
    
    fp_dat$catpo <- rep(median(fp_dat$fp),length(fp_dat$patientid))
    p1 <- ggplot(fp_dat,aes(patientid))+
      geom_ribbon(aes(ymin = catpo, ymax = fp, fill = RiskGroup), alpha = 1)+
      scale_fill_manual(values = c("#FF5733", "#33FF57")) +
      #scale_fill_jco()+
      #scale_colour_manual(values = c("#FF6666","#00CCCC"))+#"#FF5733", "#33FF57"
      theme_bw()+labs(x="Patient ID(Increasing Risk Score)",y="Risk Score")+
      geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
      # geom_hline(yintercept=res.cut$cutpoint[,1],colour="black", linetype="dotted",size=0.8)+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),
                 colour="black", linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p1
    
    pal <- pal_jco('default')(10)
    p2 <- ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=Status),size = 3)+theme_bw()+
      scale_colour_manual(values = c("#FF5733", "#33FF57"))+
      labs(x="Patient ID(Increasing Risk Score)",y="Survival Time(Years)")+
      geom_vline(xintercept=sum(fp_dat$RiskGroup=="Low"),colour="black", 
                 linetype="dotted",size=0.8) +
      theme(axis.title.x=element_text(size = 14),axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
            legend.text=element_text(size=14),legend.title =element_text(size=14))
    p2
    
    # mycolors <- colorRampPalette(c("#64b5f6", "#fffde7", "#ff5252"), bias = 1.2)(100)
    # mycolors <- colorRampPalette(c("#FF6666", "White", "#00CCCC"), bias = 1.2)(100)
    mycolors <- colorRampPalette(c("#FF5733", "White","#33FF57"), bias = 1.2)(100)
    tmp=t(scale(exp_dat))
    tmp[tmp > 1] = 1
    tmp[tmp < -1] = -1
    p3=pheatmap(tmp,col= mycolors,show_colnames = F,cluster_cols = F,
                show_rownames = T,cluster_rows = F,
                fontsize_row = 14,fontsize = 14)
    p3
    plots = list(p1,p2,as.ggplot(as.grob(p3)))
    lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7)))
    grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))
    # plots = list(p1,as.ggplot(as.grob(p3)))
    # lay1 = rbind(c(rep(1,7)),c(rep(2,7)))
    # grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10))
  }
}

#列线图----------
# nomaogram(诺莫图)
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2) 
library(survival)
library(rms)

load(file = "/home/datahup/syj/ESCC/机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
#train_phe = sv1[,c(2,6,7,9)]
colnames(TG_phe)
train_phe=TG_phe[,c(1,2,5,6,8,9)]
train_phe$id = rownames(train_phe)
risk = data.frame(row.names = sv2$Id, 
                  RiskScore = sv2$exp, 
                  Id = sv2$Id)
rownames(train_phe) <- train_phe$Id
train_phe <-train_phe[rownames(risk),]
train_phe = merge(train_phe, risk , by = 'Id')
rownames(train_phe) = train_phe$Id
train_phe = train_phe[,-c(1,7)]
colnames(train_phe)
# 自定义排序顺序
custom_order <- c("Age", "gender", "tnm_stage", "OS", "OS.time", "RiskScore")
train_phe <- train_phe[, custom_order]
train_phe$OS <- ifelse(train_phe$OS == "Alive",0,1)
save(train_phe,file = "train_phe.RData")

{
  Nmdat <- train_phe
  # Nmdat <- Nmdat[,c('RiskScore','OS','OS.time')]
  # Nmdat$OS.time <- as.numeric(Nmdat$OS.time)*30
  dd <- datadist(Nmdat)
  options(datadist="dd")
  multivarl <- as.formula(paste0('Surv(OS.time,OS)~', 
                                 paste('RiskScore', sep = '', collapse = '+')))
  
  coxm_1 <- cph(formula = multivarl,data=Nmdat,surv=T,x=T,y=T,time.inc = 365)
  surv <- Survival(coxm_1)
  surv1 <- function(x) surv(1*365,x)
  surv3 <- function(x) surv(3*365,x)
  surv5 <- function(x) surv(5*365,x)
  nomo <- nomogram(coxm_1,
                   fun = list(surv1,surv3,surv5),
                   lp = T,
                   funlabel = c('1-year survival Probability',
                                '3-year survival Probability',
                                '5-year survival Probability'),
                   maxscale = 100,
                   fun.at = c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
 
  #pdf(file = '../../Fig/Step4-1/nomaogram.pdf',width = 6,height = 5)
  plot(nomo,
       lplabel = 'Linear Preadictor',
       xfrac = .35,
       varname.label = T,
       varname.label.sep = '=',
       ia.space = .2,
       tck = NA,
       tcl = 0.2,
       lmgp = 0.3,
       points.label = 'Points',
       total.points.label = 'Total Points',
       total.sep.page = F,
       cap.labels = F,cex.var = 0.53,cex.axis = 0.53,lwd = 0.53,
       label.every = 1,col.grid = gray(c(0.8,0.95)))
  dev.off()
  
  ## 校准曲线
  ## 参数说明：
  ## 1、绘制校正曲线前需要在模型函数中添加参数x=T, y=T，详细内容参考帮助
  ## 2、u需要与之前模型中定义好的time.inc一致，即365或730；
  ## 3、m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）
  ## 而m代表每组的样本量数，因此m*3应该等于或近似等于样本量；
  ## 4、b代表最大再抽样的样本量
  f1 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=365)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1*365, m=64, B=1000)
  f3 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=3*365)
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3*365, m=64, B=1000)
  f5 <- cph(formula = multivarl, x=T, y=T, surv=T, data=Nmdat, time.inc=5*365)
  cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5*365, m=64, B=1000)
  
  #pdf(file = '../Figer/nomaogram_calibrate.pdf',width = 6,height = 5)
  plot(cal1,lwd=1,lty=1, cex.axis = 1,cex.lab=1,
       errbar.col = 'blue',
       xlab='Nomogram-Predicted Probability',
       ylab='Actual',
       col = 'blue',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal3,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = 'green',
       col = 'green',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  plot(cal5,lwd=1,lty=1, cex.axis = 1,cex.lab=1,add=T,
       errbar.col = '#FF0033',
       col = '#FF0033',subtitles = F,
       xlim = c(0,1),ylim = c(0,1))
  abline(0,1,lty=1,lwd=1)
  legend("bottomright",legend=c("1 - year","3 - year","5 - year"), 
         col=c("blue","green","#FF0033"),
         lty= 1 ,lwd= 4,
         bty = "n",
         seg.len=1,cex=1)
  dev.off()
}


#美化列线图
{
  library(survival)
  colnames(train_phe)
  Coxfit<-coxph(Surv(OS.time,OS==1)~Age+gender+tnm_stage+RiskScore,data=train_phe)  
  library(regplot)
  regplot(Coxfit, plots=c("density","boxes") ,
          observation=FALSE,points=TRUE,
          title="Survival Nomogram",
          failtime=c(366,731,1827),
          dencol="#ADD8E6",boxcol="#FFCCCB",droplines=TRUE)
  regplot(Coxfit, plots=c("bean","boxes"), 
          observation=FALSE,points=TRUE,#observation=survdata[20,],
          title="Survival Nomogram", failtime=c(365,730,1826),
          prfail=F, clickable=TRUE,  
          dencol="green", boxcol="yellow",droplines=TRUE)
}








































































