
rm(list=ls())
options(stringsAsFactors = F)
setwd("/home/datahup/syj/ESCC/下载geo数据-1.2/")

####下载bulk数据 GSE53625####
#library(GEOquery)
#dir.create("./GSE53625",recursive = T)
# eSet <- getGEO('GSE53625', destdir="./GSE53625",
#                AnnotGPL = F,
#                getGPL = F)
#save(eSet,file='./GSE53625/GSE53625_eSet.Rdata')

load(file='/home/datahup/syj/NSCC/step6_1TCGA_GEO/GSE53625/GSE53625_eSet.Rdata')
b = eSet[[1]]
raw_exprSet= exprs(b)
raw_exprSet[1:4,1:4]
phe=pData(b)

####处理临床数据------
colnames(phe)
GSE53625_phe <- data.frame(row.names = phe$geo_accession,
                           group = phe$title,
                           Age = phe$characteristics_ch1.1, 
                           gender = phe$characteristics_ch1.2, 
                           tobacco = phe$characteristics_ch1.3,#吸烟
                           alchol = phe$characteristics_ch1.4,#酒精
                           tumor_loation = phe$characteristics_ch1.5,
                           tumor_grade = phe$characteristics_ch1.6,
                           t_stage = phe$characteristics_ch1.7,
                           n_stage = phe$characteristics_ch1.8,
                           tnm_stage = phe$characteristics_ch1.9,
                           OS = phe$characteristics_ch1.14,
                           OS.time = phe$characteristics_ch1.15)#不是按天算的（月）
head(GSE53625_phe$group)
GSE53625_phe$group_1 <- sapply(strsplit(GSE53625_phe$group, " "), `[`, 1)
head(GSE53625_phe$Age)
#GSE53625_phe$Age_1 <- sub("age: (\\d+)\\..*", "\\1", GSE53625_phe$Age)
GSE53625_phe$Age_1 <- sub("age: ", "", GSE53625_phe$Age)
GSE53625_phe$Age_1 <- as.numeric(sub("\\..*$", "", GSE53625_phe$Age_1))
head(GSE53625_phe$gender)
GSE53625_phe$gender_1 <- sub("Sex: ", "", GSE53625_phe$gender)
head(GSE53625_phe$tobacco)#烟草
GSE53625_phe$tobacco_1 <- sub("tobacco use: ", "", GSE53625_phe$tobacco)
head(GSE53625_phe$alchol)#酒精
GSE53625_phe$alchol_1 <- sub("alcohol use: ", "", GSE53625_phe$alchol)
head(GSE53625_phe$tumor_loation)
GSE53625_phe$tumor_loation_1 <- sub("tumor loation: ", "", GSE53625_phe$tumor_loation)
head(GSE53625_phe$tumor_grade)
GSE53625_phe$tumor_grade_1 <- sub("tumor grade: ", "", GSE53625_phe$tumor_grade)
head(GSE53625_phe$t_stage)
GSE53625_phe$t_stage_1 <- sub("t stage: ", "", GSE53625_phe$t_stage)
head(GSE53625_phe$n_stage)
GSE53625_phe$n_stage_1 <- sub("n stage: ", "", GSE53625_phe$n_stage)
head(GSE53625_phe$tnm_stage)
GSE53625_phe$tnm_stage_1 <- sub("tnm stage: ", "", GSE53625_phe$tnm_stage)
head(GSE53625_phe$OS)
GSE53625_phe$OS_1 <- sub("death at fu: ", "", GSE53625_phe$OS)
head(GSE53625_phe$OS.time)
#GSE53625_phe$OS.time_1 <- sub("survival time(months): ", "", GSE53625_phe$OS.time)
GSE53625_phe$OS.time_1 <- as.numeric(sub(".*: ([0-9.]+)", "\\1", GSE53625_phe$OS.time))
head(GSE53625_phe$OS.time_1)
days_per_month <- 30.44
GSE53625_phe$OS.time_2 <- round(GSE53625_phe$OS.time_1 * days_per_month)

colnames(GSE53625_phe)
GSE53625_phe_1 <- GSE53625_phe[, c(13:25)]

phe53625 <- GSE53625_phe_1[,-c(4:7,12)]
colnames(phe53625) <- gsub("_1$", "", colnames(phe53625))
colnames(phe53625)[8] <- "OS.time"
phe53625$OS <- ifelse(phe53625$OS == "yes", "Dead", "Alive")
save(phe53625,file = "phe53625.Rdata")
write.csv(phe53625, file = "phe53625.csv", row.names = FALSE)

#处理表达矩阵----
raw_exprSet[1:4,1:4]
dim(raw_exprSet)

#id转换
ids <- read.csv("ids_GPL18109.csv")
table(ids$geneType)
ids <- subset(ids, geneType == "mRNA")
dat = as.data.frame(raw_exprSet)
dat[1:4,1:4]
table(rownames(dat) %in% ids$probeId)
dat53625 <- dat[rownames(dat) %in% ids$probeId, ]

ids_1 <- ids[,c(1,2)]
#删除重复基因
ids_1$median=apply(dat53625,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids_1=ids_1[order(ids_1$geneName,ids_1$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids_1=ids_1[!duplicated(ids_1$geneName),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
table(duplicated(ids_1$geneName))

dat53625 <- dat53625[rownames(dat53625) %in% ids_1$probeId, ]
head(rownames(dat53625))
head(ids_1)

gene_names <- ids_1$geneName[match(rownames(dat53625), ids_1$probeId)]
rownames(dat53625) <- gene_names
save(dat53625,file = "dat53625.Rdata")
write.csv(dat53625, file = "dat53625.csv", row.names = FALSE)









