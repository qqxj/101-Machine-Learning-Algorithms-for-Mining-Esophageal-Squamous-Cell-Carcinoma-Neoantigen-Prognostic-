#5个基因的GSVA通路评分分析
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

#load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
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
for (i in 1:length(gene)) {
  df = df[order(df[,i]),]
  a = quantile(df[,i],c(0.3,0.7))
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
  
  write.csv(diff_gene_deseq2_df,file = paste0(gene[i],'_DEG.csv'))
  
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
  
  this_tile <- paste0('Cutoff for logFC is 1',
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

#PMEPA1----
colnames(df)
{
df = df[order(df[,1]),]
a = quantile(df[,1],c(0.3,0.7))
df[,paste0(gene[1],'_group')] = ifelse(df[,1] <= a[1],'Low',
                                       ifelse(df[,1] >= a[2],'High', 'no'))
table(df$PMEPA1_group)
df1 = subset(df,df[,paste0(gene[1],'_group')] == 'High' | df[,paste0(gene[1],'_group')] == 'Low')

# deg
df2 = read.csv(file = paste0(gene[1],'_DEG.csv'))
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
#deg_entr$ENTREZID <- as.numeric(deg_entr$ENTREZID)
# KEGG
kk_deg <- enrichKEGG(gene = deg_entr$ENTREZID,
                     organism = 'hsa',
                     pvalueCutoff = 0.5,
                     qvalueCutoff =0.5)
# kk_deg <- enrichKEGG(gene         =  deg_entr$ENTREZID,
#                       organism     = 'hsa', 
#                       pvalueCutoff = 0.05,
#                       qvalueCutoff =0.05)
kk_list[[gene[1]]] <- kk_deg


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
#kegg-GSVA
library(GSVA)
meta <- data.frame(group = df1[,paste0(gene[1],'_group')], row.names = rownames(df1))
df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
kegg_gsva <- gsva(df3, kegg_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
kegg_gsva_list[[gene[1]]] <- kegg_gsva


# GO
go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                   OrgDb          = org.Hs.eg.db,
                   ont            = 'ALL', 
                   pAdjustMethod  = "BH",
                   pvalueCutoff   = 0.05, 
                    qvalueCutoff   = 0.2, 
                   readable       = TRUE)
go_list[[gene[1]]] <- go_deg
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
                parallel.sz=10) #gsva
go_gsva_list[[gene[1]]] <- go_gsva

#print(paste0(gene[1],' is over'))
}

#MAGEA4----
colnames(df)
{
  df = df[order(df[,2]),]
  a = quantile(df[,2],c(0.3,0.7))
  df[,paste0(gene[2],'_group')] = ifelse(df[,2] <= a[1],'Low',
                                         ifelse(df[,2] >= a[2],'High', 'no'))
  #table(df$paste0(gene[6],'_group'))
  table(df[,paste0(gene[2],'_group')])  
  df1 = subset(df,df[,paste0(gene[2],'_group')] == 'High' | df[,paste0(gene[2],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[2],'_DEG.csv'))
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
                       pvalueCutoff = 0.7,
                       qvalueCutoff =0.7)
  kk_list[[gene[2]]] <- kk_deg
  
  
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
  meta <- data.frame(group = df1[,paste0(gene[2],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[2]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.5, 
                     qvalueCutoff   = 1, 
                     readable       = TRUE)
  go_list[[gene[2]]] <- go_deg
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
  go_gsva <- gsva(df3, go_genelist, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) #gsva
  go_gsva_list[[gene[2]]] <- go_gsva
  
  #print(paste0(gene[6],' is over'))
}

#RCN1----
colnames(df)
{
  df = df[order(df[,3]),]
  a = quantile(df[,3],c(0.3,0.7))
  df[,paste0(gene[3],'_group')] = ifelse(df[,3] <= a[1],'Low',
                                         ifelse(df[,3] >= a[2],'High', 'no'))
  #table(df$paste0(gene[6],'_group'))
  table(df[,paste0(gene[3],'_group')])  
  df1 = subset(df,df[,paste0(gene[3],'_group')] == 'High' | df[,paste0(gene[3],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[3],'_DEG.csv'))
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
                       pvalueCutoff = 0.6,
                       qvalueCutoff =0.6)
  kk_list[[gene[3]]] <- kk_deg
  
  
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
  meta <- data.frame(group = df1[,paste0(gene[3],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist, 
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[3]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[3]]] <- go_deg
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
  #gsva
  go_gsva <- gsva(df3, go_genelist,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=10) 
  go_gsva_list[[gene[3]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
}

#DLX5----
colnames(df)
{
  df = df[order(df[,4]),]
  a = quantile(df[,4],c(0.3,0.7))
  df[,paste0(gene[4],'_group')] = ifelse(df[,4] <= a[1],'Low',
                                         ifelse(df[,4] >= a[2],'High', 'no'))
  table(df[,paste0(gene[4],'_group')])  
  df1 = subset(df,df[,paste0(gene[4],'_group')] == 'High' | df[,paste0(gene[4],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[4],'_DEG.csv'))
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
                       pvalueCutoff = 0.4,
                       qvalueCutoff =0.6)
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

#TIMP1----
colnames(df)
{
  df = df[order(df[,5]),]
  a = quantile(df[,5],c(0.3,0.7))
  df[,paste0(gene[5],'_group')] = ifelse(df[,5] <= a[1],'Low',
                                         ifelse(df[,5] >= a[2],'High', 'no'))
  table(df[,paste0(gene[5],'_group')])  
  df1 = subset(df,df[,paste0(gene[5],'_group')] == 'High' | df[,paste0(gene[5],'_group')] == 'Low')
  
  # deg
  df2 = read.csv(file = paste0(gene[5],'_DEG.csv'))
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
                       pvalueCutoff = 0.05)
  kk_list[[gene[5]]] <- kk_deg
  
  
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
  meta <- data.frame(group = df1[,paste0(gene[5],'_group')], row.names = rownames(df1))
  df3 <- as.matrix(t(dat0[rownames(df1),rownames(df2)]))
  #gsva
  kegg_gsva <- gsva(df3, kegg_genelist,
                    kcdf="Gaussian",
                    method = "gsva",
                    parallel.sz=10) 
  kegg_gsva_list[[gene[5]]] <- kegg_gsva
  
  
  # GO
  go_deg <- enrichGO(gene           = deg_entr$ENTREZID, 
                     OrgDb          = org.Hs.eg.db,
                     ont            = 'ALL', 
                     pAdjustMethod  = "BH",
                     pvalueCutoff   = 0.05, 
                     # qvalueCutoff   = 0.2, 
                     readable       = TRUE)
  go_list[[gene[5]]] <- go_deg
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
  go_gsva_list[[gene[5]]] <- go_gsva
  
  #print(paste0(gene[11],' is over'))
}

save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,
     file = 'Step6-1 Result_1.Rdata')


# 6.28 GSVA得分差异分析 找到关键通路--------
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

load(file = 'Step6-1 Result_1.Rdata')
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
gene = rownames(gene)

KEGG_logFC_list = list()
GO_logFC_list = list()

#go和kegg的gsva差异分析(循环,跳过挨个跑)-------
for (i in gene) {
  
  kk_re = kk_list[[i]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[i]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  
  go_re = go_list[[i]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[i]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  
  group = data.frame(row.names = rownames(df), group = df[,paste0(i,'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  
  # gsva矩阵差异分析
  library(limma)
  library(stringr)
  # 分组矩阵
  Group = data.frame(row.names = rownames(group), group = group$group)
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
  KEGG_logFC_list[[i]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[i]] = GO_logFC
  
  
}
# save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,
#      KEGG_logFC_list,GO_logFC_list,
#      file = 'Step6-1 Result.Rdata')

library(stringr)
#PMEPA1------
gene
{
  kk_re = kk_list[[1]]@result
  # kk_re_1 = as.data.frame(kk_list[["PMEPA1"]])
  # table(kk_re$Description %in% kk_re_1$Description)
  
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[1]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[1]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[1]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("PMEPA1",'_group')], Sample = rownames(df))
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
  KEGG_logFC_list[[1]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[1]] = GO_logFC
  
}

#MAGEA4-------
gene
{
  kk_re = kk_list[[2]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[2]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[2]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[2]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  go_gsva_df <- na.omit(go_gsva_df)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("MAGEA4",'_group')], Sample = rownames(df))
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
  KEGG_logFC_list[[2]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[2]] = GO_logFC
  
}

#RCN1------
gene
{
  kk_re = kk_list[[3]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[3]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[3]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[3]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("RCN1",'_group')], Sample = rownames(df))
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
  KEGG_logFC_list[[3]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[3]] = GO_logFC
  
}

#DLX5-------
gene
{
  kk_re = kk_list[[4]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[4]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[4]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[4]])
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

#TIMP1------
gene
{
  kk_re = kk_list[[5]]@result
  kk_re = kk_re %>% filter(pvalue < 0.05)
  kk_gsva_df = as.data.frame(kegg_gsva_list[[5]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = kk_gsva_df[kk_re$Description,]
  kk_gsva_df <- na.omit(kk_gsva_df)
  
  go_re = go_list[[5]]@result
  go_re = go_re %>% filter(pvalue < 0.05)
  go_gsva_df = as.data.frame(go_gsva_list[[5]])
  go_gsva_df = go_gsva_df[go_re$Description,]
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("TIMP1",'_group')], Sample = rownames(df))
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
  KEGG_logFC_list[[5]] = KEGG_logFC
  
  # GO 差异分析
  GO_logFC = deg(go_gsva_df,design,contrast.matrix)
  GO_logFC_list[[5]] = GO_logFC
  
}

save(df,kk_list,kegg_gsva_list,go_list,go_gsva_list,
     KEGG_logFC_list,GO_logFC_list,
     file = 'Step6-1 Result.Rdata')

# 作图--------
rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/11个基因的GSVA分析-9.4/')

load(file = 'Step6-1 Result.Rdata')
load(file = "../机器学习和验证然后取跑100种机器学习-8.2/lasso结果.Rdata")
library(dplyr)
library(tidyr)
library(tibble)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(stringr)

# PMEPA1-----
gene
{
# KEGG
KEGG_logFC = KEGG_logFC_list[[1]]
kk_gsva_df = as.data.frame(kegg_gsva_list[["PMEPA1"]])
rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))

colnames(kk_gsva_df)
kk_gsva_df <- kk_gsva_df[,c("ABC transporters",#
                            "Arachidonic acid metabolism",#
                            #"Axon guidance",
                            "Basal cell carcinoma",#
                            "cAMP signaling pathway",#
                            "Cytokine-cytokine receptor interaction",#
                            #"Cytoskeleton in muscle cells",
                            #"GABAergic synapse",
                            "Hepatocellular carcinoma",#
                            "Linoleic acid metabolism",#
                            #"Melanogenesis",
                            "Neuroactive ligand-receptor interaction",#
                            "Phospholipase D signaling pathway",#
                            "Proteoglycans in cancer",#
                            "Signaling pathways regulating pluripotency of stem cells",
                            "Wnt signaling pathway")]#

kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

group = data.frame(row.names = rownames(df), group = df[,paste0("PMEPA1",'_group')], Sample = rownames(df))
group = group[group$group != 'no',]
table(group$group)

kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(kk_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in PMEPA1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
#ggsave(filename = 'PMEPA1 kegg gsva.pdf',width = 10,height = 8,path = './')

#格式800*1200
plot1 <-ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
  # split violin
  geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw() + 
  labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','PMEPA1',' Grouping')) +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.3,label = 'p.signif') + # wilcox.test t.test
  coord_flip() 
# 保存为PDF文件
ggsave(filename = './KEGG_PMEPA1_Grouping.pdf', 
       plot = plot1, width = 8, height = 12)
# GO
# go_re = go_list[["PMEPA1"]]@result
# go_re = go_re %>% filter(pvalue < 0.05)

GO_logFC = GO_logFC_list[[1]]
GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.25)

go_gsva_df = as.data.frame(go_gsva_list[["PMEPA1"]])
go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))

colnames(go_gsva_df)
go_gsva_df <- go_gsva_df[,c(
  "extracellular matrix organization",
  "collagen-containing extracellular matrix",
  "Wnt-protein binding",
  "transmembrane receptor protein kinase activity",
  "transmembrane receptor protein tyrosine kinase activity",
  "xenobiotic metabolic process",
  "retinol binding",
  "xenobiotic metabolic process",
  "aromatase activity"
)]

go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
table(go_gsva_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in PMEPA1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
#ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)

#格式800*1200
plot6 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  # split violin
  geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw() + 
  labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','PMEPA1',' Grouping')) +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  stat_compare_means(aes(group = group),method = "t.test",label.y = 1.3,label = 'p.signif') + # wilcox.test
  coord_flip() 
ggsave(filename = './GO_PMEPA1.pdf',plot6,width = 8,height = 12)


#将两者合并
colnames(kk_gsva_df)
colnames(kk_gsva_df) <- paste0("KEGG_", colnames(kk_gsva_df))
colnames(go_gsva_df)
colnames(go_gsva_df) <- paste0("GO_", colnames(go_gsva_df))
kk_gsva_df <- kk_gsva_df[rownames(go_gsva_df),]
merged_df <- merge(kk_gsva_df, go_gsva_df, by = "row.names")
rownames(merged_df) <- merged_df$Row.names
merged_df <- merged_df[,-1]

merged_df2 <- merged_df %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Pathway,value = `GSVA Score`,-Sample)

merged_df2$group = ifelse(merged_df2$Sample %in% group$Sample, group$group, 'unkown')
table(merged_df2$group)

mypalette = pal_ucscgb()(10)
ggplot(merged_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in PMEPA1 Grouping') +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
#ggsave(filename = 'SAP18 go gsva.pdf',width = 10,height = 8,path = './')

library(ggplot2)
library(ggunchained) 
# library(devtools)  
# install_github("JanCoUnchained/ggunchained")

mypalette = pal_ucscgb()(26)
#格式800*1200
ggplot(merged_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
  # split violin
  geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw() + 
  labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG/GO Pathway in ','PMEPA1',' Grouping')) +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
        plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  stat_compare_means(aes(group = group),method = "t.test",label.y = 1.1,label = 'p.signif') + # wilcox.test
  coord_flip() 
}


# MAGEA4 ----
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[2]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["MAGEA4"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  #kk_gsva_df <- kk_gsva_df[, colSums(is.na(kk_gsva_df)) == NA]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("MAGEA4",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in CSRP1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'CSRP1 kegg gsva.pdf',width = 10,height = 8,path = './')
 
  # GO
  # go_re = go_list[["MAGEA4"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[2]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 )#& abs(logFC) > 0.25
  
  go_gsva_df = as.data.frame(go_gsva_list[["CSRP1"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  go_gsva_df2 <- go_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  go_gsva_df2$group = ifelse(go_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(go_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(go_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "GO Terms", y = "GSVA Score", title = 'GO Pathway GSVA Score in CSRP1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6,3)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'CSRP1 go gsva.pdf',width = 10,height = 8,path = './')
  
}

# RCN1------
gene
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[3]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["RCN1"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  kk_gsva_df <- kk_gsva_df[,c(
    "Signaling pathways regulating pluripotency of stem cells",
    "Cell adhesion molecules",
    "Ether lipid metabolism",
    "IL-17 signaling pathway",
    "Linoleic acid metabolism",
    "Renin secretion",
    "Cytokine-cytokine receptor interaction",
    "Neuroactive ligand-receptor interaction",
    "Phospholipase D signaling pathway",
    "Arachidonic acid metabolism"
)]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("RCN1",'_group')], Sample = rownames(df))
  group = group[group$group != 'no',]
  table(group$group)
  
  kk_gsva_df2$group = ifelse(kk_gsva_df2$Sample %in% group$Sample, group$group, 'unkown')
  table(kk_gsva_df2$group)
  
  mypalette = pal_ucscgb()(10)
  ggplot(kk_gsva_df2,aes(Pathway,`GSVA Score`,fill = group)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + 
    theme_bw() + 
    labs(x = "KEGG Terms", y = "GSVA Score", title = 'KEGG Pathway GSVA Score in RCN1 Grouping') +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=75,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'GFRA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  plot2 <-ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','RCN1',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    stat_compare_means(aes(group = group),method = "t.test",label.y = 1.5,label = 'p.signif') + # wilcox.test
    coord_flip()  
  ggsave(filename = './KEGG_RCN1_Grouping.pdf',
         plot = plot2, width = 8, height = 12)
  
  # GO
  # go_re = go_list[["GFRA1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[3]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.12)
  
  go_gsva_df = as.data.frame(go_gsva_list[["RCN1"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
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
  plot7 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`, fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','RCN1',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    stat_compare_means(aes(group = group),method = "t.test",label.y = 1.4,label = 'p.signif') + # wilcox.test
    coord_flip()  
  ggsave(filename = './GO_RCN1.pdf', plot = plot7, width = 8, height = 12)
}

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
  
  # GO
  # go_re = go_list[["DLX5"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[4]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.12)
  
  go_gsva_df = as.data.frame(go_gsva_list[["DLX5"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
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
  
  
}

# TIMP1---------
gene$sig_gene_multi_cox
{
  # KEGG
  KEGG_logFC = KEGG_logFC_list[[5]]
  kk_gsva_df = as.data.frame(kegg_gsva_list[["TIMP1"]])
  rownames(kk_gsva_df) = str_sub(rownames(kk_gsva_df), start = 1, end = -24)
  kk_gsva_df = as.data.frame(t(kk_gsva_df[rownames(KEGG_logFC),]))
  
  kk_gsva_df <- kk_gsva_df[,c(
    "Wnt signaling pathway",
    "Chemokine signaling pathway",
    "Viral protein interaction with cytokine and cytokine receptor",
    "Complement and coagulation cascades",
    "Cytokine-cytokine receptor interaction",
    "Thyroid hormone synthesis",
    "Aldosterone synthesis and secretion",
    "Cell adhesion molecules",
    "Linoleic acid metabolism",
    "Neuroactive ligand-receptor interaction",
    "Phospholipase D signaling pathway",
    "Signaling pathways regulating pluripotency of stem cells",
    "Arachidonic acid metabolism",
    "Inflammatory mediator regulation of TRP channels"
  )]
  
  kk_gsva_df2 <- kk_gsva_df %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Pathway,value = `GSVA Score`,-Sample)
  
  group = data.frame(row.names = rownames(df), group = df[,paste0("TIMP1",'_group')], Sample = rownames(df))
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
    scale_fill_manual(values = mypalette[c(1,6)])+ stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.1,label = 'p.signif')
  #ggsave(filename = 'GFRA1 kegg gsva.pdf',width = 10,height = 8,path = './')
  
  #格式800*1200
  plot3 <- ggplot(kk_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of KEGG Pathway in ','TIMP1',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1.5,label = 'p.signif') + # wilcox.test t.test
    coord_flip() 
  ggsave(filename = './KEGG_TIMP1_Grouping.pdf', plot = plot3, width = 8, height = 12)

  # GO
  # go_re = go_list[["TIMP1"]]@result
  # go_re = go_re %>% filter(pvalue < 0.05)
  
  GO_logFC = GO_logFC_list[[5]]
  GO_logFC1 = GO_logFC %>% filter(P.Value < 0.05 & abs(logFC) > 0.12)
  
  go_gsva_df = as.data.frame(go_gsva_list[["TIMP1"]])
  go_gsva_df = as.data.frame(t(go_gsva_df[rownames(GO_logFC1),]))
  
  colnames(go_gsva_df)
  go_gsva_df <- go_gsva_df[,c(
    "collagen-containing extracellular matrix",
    "extracellular matrix structural constituent",
    "positive regulation of angiogenesis",
    "regulation of angiogenesis",
    "regulation of vasculature development",
    "immune receptor activity",
    "cytokine receptor activity",
    "chemokine binding",
    "chemokine-mediated signaling pathway",
    "chemokine receptor activity",
    "cytokine activity",
    "cytokine receptor binding"
  )]
  
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
  plot8 <- ggplot(go_gsva_df2,aes(x = Pathway, y = `GSVA Score`,   fill = group)) +
    # split violin
    geom_split_violin(alpha = 1, trim = F,color = NA,width = 1.2) +
    # mean point
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw() + 
    labs(x = " ", y = "GSVA Score", title = paste0('GSVA Score of GO Pathway in ','TIMP1',' Grouping')) +
    theme(legend.position = "right") + 
    theme(axis.text.x = element_text(angle=0,vjust = 0.5,size = 8,color = "black", hjust = 0.5),
          plot.title = element_text(size = 12,colour = "black",face = "bold",vjust = 1.9,hjust = 0.5))+
    scale_fill_manual(values = c("#F8766D","#00BFC4"))+
    stat_compare_means(aes(group = group),method = "wilcox.test",label.y = 1,label = 'p.signif') + # wilcox.test t.test
    coord_flip() 
  ggsave(filename = './GO_TIMP1.pdf', plot = plot8, width = 8, height = 12)
}








