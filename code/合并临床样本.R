rm(list=ls())
options(stringsAsFactors = F)
#加载
setwd("/home/datahup/syj/NSCC/step6_1TCGA_GEO/")
load(file = "./TCGA_output.Rdata")
setwd("/home/datahup/syj/ESCC/合并数据集(临床样本)-1.4/")
load(file = "../../ESCC/下载geo数据-1.2/phe53625.Rdata")

#TCGA数据-------
colnames(phe6)
TCGA_phe <- phe6[,c(2:11,13)]
table(TCGA_phe$Stage)
TCGA_phe <- subset(TCGA_phe, Stage != "not reported")
TCGA_phe$Stage_Converted <- ifelse(TCGA_phe$Stage %in% c("stage i", "stage ia", "stage ib"), "I", 
                                   ifelse(TCGA_phe$Stage %in% c("stage iia", "stage iib"), "II", 
                                          ifelse(TCGA_phe$Stage %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic"), "III", "IV")))

TCGA_phe$OS.time <- ifelse(is.na(TCGA_phe$`Days to Death`),
                           TCGA_phe$`Days to Last Follow`,TCGA_phe$`Days to Death`)
TCGA_phe <- TCGA_phe[,-c(3,7,8,10)]
colnames(TCGA_phe) <- c("Id","Age","n_stage","t_stage","gender","OS","group","tnm_stage","OS.time")
TCGA_phe$group <- ifelse(TCGA_phe$group == "normal", "normal", "cancer")

phe53625$Id <- rownames(phe53625)
phe53625 <- phe53625[, match(colnames(TCGA_phe), colnames(phe53625))]

merge_phe <- rbind(TCGA_phe,phe53625)

save(merge_phe,file = "合并的临床数据.Rdata")
write.csv(merge_phe, file = "合并临床数据.csv", row.names = FALSE)


setwd("/home/datahup/syj/ESCC/合并数据集(临床样本)-1.4/")
load(file = "合并的临床数据.Rdata")
#匹配dat数据----------
load(file = "../../ESCC/合并数据集-1.3/合并数据.Rdata")
merge_dat[1:4,1:4]
merge_dat <- as.data.frame(merge_dat)
merge_dat <- merge_dat[,merge_phe$Id]
save(merge_dat,file = "合并表达矩阵.Rdata")
write.csv(merge_dat, file = "合并表达矩阵.csv", row.names = FALSE)

#拆出来TCGA的数据--------
rm(list=ls())
setwd("/home/datahup/syj/ESCC/合并数据集(临床样本)-1.4/")
load(file = "合并表达矩阵.Rdata")
load(file = "合并的临床数据.Rdata")
merge_phe$Id

TCGA_phe <- subset(merge_phe, grepl("^TCGA", Id))
TCGA_dat <- merge_dat[,TCGA_phe$Id]
save(TCGA_phe,TCGA_dat,file = "tcga_dat_phe.Rdata")




