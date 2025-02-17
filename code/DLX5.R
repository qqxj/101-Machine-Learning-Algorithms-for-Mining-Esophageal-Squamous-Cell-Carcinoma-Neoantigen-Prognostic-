#DLX5
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

load(file = "../../ESCC/下载geo数据-1.2/dat53625.Rdata")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

#基因
gene_multi <- gene
gene = rownames(gene_multi)
#gene <- c("SAP18","CSRP1","GFRA1")
#基因矩阵
NSCLCcount <- dat53625
NSCLCcount[1:6,1:6]
dat = as.data.frame(t(NSCLCcount))
dat[1:6,1:6]
phe <- subset(phe53625,group == "cancer")
dat0 <- dat[rownames(phe),]
df = data.frame(gene = dat0[,gene],
                row.names = rownames(dat0)) 

#Step1 高低分组做差异分析----------
library(DESeq2)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')
for (i in 4) {
  df = df[order(df[,i]),]
  a = quantile(df[,i],c(0.5,0.5))
  df$group = ifelse(df[,i] <= a[1],'Low',
                    ifelse(df[,i] >= a[2],'High', 'no'))
  table(df$group)  
  df1 = subset(df,df$group == 'High' | df$group == 'Low')
  
  #构建dds矩阵
  dat1 = as.data.frame(t(dat0[rownames(df1),]))
  range(dat1)
  dat1 <- 2^dat1 -1 #转Count
  dat1 <- round(dat1,digits = 0) # 取整
  dat1[dat1 < 0 ] = 0
  countData <- dat1
  countData[1:6,1:6]
  Group <- data.frame(group = df1$group, row.names = rownames(df1))
  condition <- factor(Group$group)
  table(condition)
  head(condition)
  # 差异分析
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition )
  head(dds)
  dim(dds)
  dds <- DESeq(dds) 
  resultsNames(dds)
  res <- results(dds)
  summary(res)
  DEG_INFO <- as.data.frame(res)
  # 提取差异分析结果
  table(res$padj<0.05) #取P值小于0.05的结果
  res <- res[order(res$padj),]
  resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
  rownames(resdata) <- resdata$Row.names
  resdata <- resdata[,-1]
  resdata <- na.omit(resdata)
  diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  diff_gene_deseq2 <- row.names(diff_gene_deseq2)# 所有差异基因
  diff_gene_deseq2_df <- resdata[diff_gene_deseq2,1:6]
  diff_gene_deseq2_df$state <- ifelse(diff_gene_deseq2_df$log2FoldChange > 1, 'Up','Down')
  table(diff_gene_deseq2_df$state)
  
  write.csv(diff_gene_deseq2_df,file = paste0(gene[i],'_DEG_DLX5.csv'))
  
  # 火山图
  library(ggplot2)
  for_volcano <- data.frame('log2FoldChange' = res$log2FoldChange,
                            'padj' = res$padj,
                            'State' = rep('No', length(res$log2FoldChange)))
  up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 1), which(for_volcano$padj < 0.05))
  down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -1), which(for_volcano$padj < 0.05))
  for_volcano[up_sig_indices,'State'] <- 'Up'
  for_volcano[down_sig_indices,'State'] <- 'Down'
  for_volcano$State <- as.factor(for_volcano$State)
  for_volcano$padj <- -log10(for_volcano$padj)
  
  this_tile <- paste0('Cutoff for logFC is 2',
                      '\nThe number of Up gene is ',nrow(diff_gene_deseq2_df[diff_gene_deseq2_df$state =='Up',]) ,
                      '\nThe number of Down gene is ',nrow(diff_gene_deseq2_df[diff_gene_deseq2_df$state =='Down',]))
  
  
  p <- ggplot(for_volcano,aes(x = log2FoldChange, y = padj, colour = State))+
    geom_point(size = I(0.7))+
    scale_color_manual(values = c('No'='black', 'Up' = 'red', 'Down' = 'blue'))+
    geom_vline(xintercept = c(1, -1), lty=2, size=I(0.4), colour = 'grey11')+
    geom_hline(yintercept = c(-log(x=0.05,base = 10)),lty=2, size=I(0.1),colour = 'grey11')+
    theme_bw()+
    theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'),
          panel.grid = element_blank())+
    labs(x='log2FoldChange', y = '-log10Pvalue')+
    ggtitle( this_tile ) +
    theme(plot.title = element_text(size=15,hjust = 0.5)) 
  
  p
  
  ggsave(filename = paste0(gene[i],'_vol.png'),plot = p,width = 8,height = 6,
         path = "./")
  
}

#  go|KEGG-GSVA--------
#rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

#kegg
kk_list = list()
kegg_gsva_list = list()
#go
go_list = list()
go_gsva_list = list()

# 加载GO通路数据库
library(clusterProfiler)
{
  source("getGoTerm.R")
  GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
  save(GO_DATA, file = "GO_DATA.RData")
  findGO <- function(pattern, method = "key"){
    
    if(!exists("GO_DATA"))
      load("GO_DATA.RData")
    if(method == "key"){
      pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
    } else if(method == "gene"){
      pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
    }
    
    colnames(pathways) = "pathway"
    
    if(length(pathways) == 0){
      cat("No results!\n")
    } else{
      return(pathways)
    }
  } # 用于寻找 GO ID
  getGO <- function(ID){
    
    if(!exists("GO_DATA"))
      load("GO_DATA.RData")
    allNAME = names(GO_DATA$PATHID2EXTID)
    if(ID %in% allNAME){
      geneSet = GO_DATA$PATHID2EXTID[ID]
      names(geneSet) = GO_DATA$PATHID2NAME[ID]
      return(geneSet)     
    } else{
      cat("No results!\n")
    }
  } # 获取 GO geneSet
  load("GO_DATA.RData") # 载入数据 GO_DATA
}

#DLX5----
colnames(df)
{
  df = df[order(df[,4]),]
  a = quantile(df[,4],c(0.5,0.5))
  df[,paste0(gene[4],'_group')] = ifelse(df[,4] <= a[1],'Low',
                                         ifelse(df[,4] >= a[2],'High', 'no'))
  table(df[,paste0(gene[4],'_group')])  
  df1 = subset(df,df[,paste0(gene[4],'_group')] == 'High' | df[,paste0(gene[4],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[4],'_DEG_DLX5.csv'))
  rownames(df2) <- df2$X
  df2 <- df2[,-1]
  
  ## KEGG 富集
  ##  ID转换
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  # degene 
  deg_entr <- bitr(rownames(df2), fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)
  # KEGG
  kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.9,
                       qvalueCutoff =0.9)
  kk_list[[gene[4]]] <- kk_deg
  
  
  # 提取通路里的基因集
  library(KEGGREST) 
  keggpathway <- kk_deg@result[["ID"]]
  kegg_genelist <- list()
  for (j in keggpathway) {
    gs <- keggGet(j)  # 使用 keggGet 函数获取人类基因信号通路 的信息，并缓存
    #获取通路中gene信息 
    # gs[[1]]$GENE 
    #查找所有基因 
    # print(paste0('Now is ',i))
    
    genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
    genelist <- genes[1:length(genes)%%3 ==2] 
    gs[[1]]$NAME # 通路名称
    
    kegg_genelist[[gs[[1]]$NAME]] <- genelist # 保存为list，方便下一步gsva
    
  }
  library(GSVA)
  meta <- data.frame(group = df1[,paste0(gene[4],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[4]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[4]]] <- go_deg
  gopathway <- go_deg@result[["ID"]]
  
  go_genelist <- list()
  # 批量获取通路基因集
  # for (x in gopathway) {
  #   go_genelist <- getGO(x)
  # }
  for (i in gopathway) {
    needline <- subset(go_deg, go_deg$ID == i)
    go_genelist[[needline$Description]] <- as.character(str_split(needline$geneID, '/', simplify = T))
  }
  #go-gsva
  go_gsva <- gsva(df3, go_genelist, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[4]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
}

# 6.28 GSVA得分差异分析 找到关键通路--------
#rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

#load(file = 'Step6-1 Result_1.Rdata')
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
gene = rownames(gene)

KEGG_logFC_list = list()
GO_logFC_list = list()

library(stringr)

#DLX5-------
gene
{
  kk_re = kk_list[[1]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[1]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[1]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[1]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("DLX5",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), 
                     group = group$group)
  table(Group$group)
  Group = as.matrix(Group)
  class(Group)
  
  design = model.matrix(~0+factor(Group))
  colnames(design) = levels(factor(Group))
  rownames(design) = rownames(Group)
  head(design)
  
  # 比较矩阵
  contrast.matrix<-makeContrasts("High-Low",
                                 levels = design)
  contrast.matrix ##这个矩阵声明，我们要把 1组跟-1进行差异分析比较
  # Contrasts
  # Levels High-Low
  # High        1
  # Low        -1
  
  # 差异分析函数
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    ##这一步很重要，大家可以自行看看效果
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  # KEGG 差异分析
  KEGG_logFC = deg(kk_gsva_df,design,contrast.matrix)
  KEGG_logFC_list[[4]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[4]] = GO_logFC
  
}

# 作图--------
#rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

#load(file = 'Step6-1 Result_MAGEA4.Rdata')
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")

library(dplyr)
library(tidyr)
library(tibble)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(stringr)

# DLX5 没有--------
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[4]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["DLX5"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("DLX5",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in GFRA1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "t.test",label.y = 1.1,label = 'p.signif')#$wilcox.test
  #ggsave(filename = 'GFRA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  plot5 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','DLX5',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.5,label = 'p.signif') + # wilcox.test t.test
    coord_flip() 
  ggsave(filename = './KEGG_DLX5.pdf', plot = plot5, width = 8, height = 12)
  # GO
  # go_re = go_list[["DLX5"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[4]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.12)
  
  go_gsva_df = as.data.frame(go_gsva_list[["DLX5"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  go_gsva_df <- go_gsva_df[,-8]
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in GFRA1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'GFRA1 go gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  plot10 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','DLX5',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.5,label = 'p.signif') + # wilcox.test t.test
    coord_flip() 
  
}




















