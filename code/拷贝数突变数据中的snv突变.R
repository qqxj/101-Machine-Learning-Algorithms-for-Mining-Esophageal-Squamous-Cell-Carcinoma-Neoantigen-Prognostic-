#9.27：SNP
rm(list = ls())
options(stringsAsFactors = F)
gc()
library(maftools)
#snv突变数据、临床数据
setwd('/home/datahup/syj/NSCC/step8_1TCGA_snv/')
lusc.snv = read.delim(gzfile('./TCGA-ESCA.mutect2_snv.tsv.gz'))
load(file = "../step6_1TCGA_GEO/TCGA_output.Rdata")
#cnv数据
setwd("/home/datahup/syj/ESCC/SNP-TCGA-6.1/")
load(file = "../../ESCC/拟时序-6/癌细胞_Statemarker.Rdata")
gene = unique(Statemarker$gene)
table(gene %in% lusc.snv$gene)
snv <- lusc.snv[lusc.snv$gene %in% gene, ]
table(duplicated(snv$gene))


#数据清洗(用不了、用NSCC文件里的就行)-----
dat <- phe6
tmp <- snv

dat$Id2 = gsub('[.]','-',dat$Id)
rownames(dat) = dat$Id2
length(unique(tmp$gene))
length(unique(tmp$Sample_ID))
phe = subset(dat,dat$Id2 %in% unique(tmp$Sample_ID))
head(tmp)   
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
table(tmp$Variant_Type )
tmp2 = subset(tmp, tmp$Tumor_Sample_Barcode  %in% dat$Id2 ) # 添加临床信息
phe$age = ifelse(phe$Age <= 60,'<=60','>60')
tmp2$age = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$age, 'Unknow')
table(tmp2$age)
tmp2$Gender = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$Gender, 'Unknow')
table(tmp2$Gender)
tmp2$`Survive State` = ifelse(tmp2$Tumor_Sample_Barcode %in% phe$Id2 , phe$OS, 'Unknow')
table(tmp2$`Survive State`)
neu.state = subset(tmp2,tmp2$Hugo_Symbol %in% tmp$Hugo_Symbol)

tcga.neu.state = read.maf(maf = tmp,
                          vc_nonSyn=names(tail(sort(table(tmp$Variant_Classification )))))

save(neu.state,tcga.neu.state,file = "./snv.Rdata")

# 添加临床信息
a = tcga.neu.state@clinical.data
a$age = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$age,'Unkown')
a$age[is.na(a$age)] = 'Unkown'
table(a$age)
a$Gender = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$Gender,'Unkown')
table(a$Gender)
a$`Survive State` = ifelse(a$Tumor_Sample_Barcode %in% tmp2$Tumor_Sample_Barcode, tmp2$`Survive State`,'Unkown')
table(a$`Survive State`)
tcga.neu.state@clinical.data = a

####瀑布图-----
oncoplot(maf = tcga.neu.state,top = 10) # 高频突变的前10个基因
plot <- oncoplot(maf = tcga.neu.state, top = 10)
ggsave(file = "./前十个基因突变——瀑布图.png", plot, width = 10, height = 6, dpi = 300)

getFields(tcga.neu.state)
getClinicalData(tcga.neu.state) #查看临床信息
getSampleSummary(tcga.neu.state)#查看每个样品发生突变的情况，此处就可以计算tumor mutation load,TML=Missense_Mutation/外显子数。

plotmafSummary(maf = tcga.neu.state, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE,
               titvRaw=FALSE)#绘制整体的突变情况

library(ggsci)
col = pal_lancet()(6)#6种
# col = RColorBrewer::brewer.pal(n = 6, name = 'Paired')
table(tcga.neu.state@data[["effect"]])
names(col) = c('frameshift_variant','synonymous_variant','missense_variant', 'intron_variant', '3_prime_UTR_variant' , 'stop_gained')

tcga.neu.state
getClinicalData(x = tcga.neu.state)

pdf(file  = './neu_oncoplot.pdf',width = 12,height = 8)
oncoplot(maf = tcga.neu.state, top = 30, colors = col,
         draw_titv = T,
         clinicalFeatures = c('age','Gender','Survive State'),
         sortByAnnotation = TRUE,
         # annotationColor = c(anocol,anocol2,anocol3)
)
dev.off()

# 转换 颠倒 
pdf(file  = './neu_titv.pdf',width = 12,height = 8)
titv(tcga.neu.state, useSyn = FALSE, plot = TRUE, file = NULL)
dev.off()

####提取snv突变数据####
snv_data <- tcga.neu.state@data[tcga.neu.state@data$Variant_Type == "SNP", ]
selected_columns <- c("Hugo_Symbol", "Chromosome", "Start_Position", 
                      "End_Position", "Reference_Allele", 
                      "Tumor_Seq_Allele2", "Variant_Classification")
snv_data <- snv_data[, selected_columns, with = FALSE]
save(snv_data,file = "./snv_data.Rdata")


# 不同State突变情况-----
rm(list = ls())
setwd("/home/datahup/syj/ESCC/SNP-TCGA-6.1/")
lusc.snv = read.delim(gzfile('/home/datahup/syj/NSCC/step8_1TCGA_snv/TCGA-ESCA.mutect2_snv.tsv.gz'))
load(file = "../../ESCC/拟时序-6/癌细胞_Statemarker.Rdata")
table(Statemarker$cluster)
table(unique(State1$gene) %in% unique(lusc.snv$gene))
table(unique(State2$gene) %in% unique(lusc.snv$gene))
table(unique(State3$gene) %in% unique(lusc.snv$gene))
# State1 
table(Statemarker$cluster)
State1 = subset(Statemarker,Statemarker$cluster == 'State1')
#备注：不知道为什么这样匹配出来的结果不能用，直接用总数居进行匹配百分比
 #table(unique(tmp2$Hugo_Symbol) %in% unique(State1$gene))
length(unique(State1$gene))
State1_pre = (420/600)*100
# neu.state1 = subset(tmp,tmp$Hugo_Symbol %in% State1$gene)
# table(neu.state1$effect)

# State2
table(Statemarker$cluster)
State2 = subset(Statemarker,Statemarker$cluster == 'State2')
 #table(unique(tmp2$Hugo_Symbol) %in% unique(State2$gene))
length(unique(State2$gene))
State2_pre = (573/857)*100
# neu.state2 = subset(tmp,tmp$Hugo_Symbol %in% State2$gene)

# State3
table(Statemarker$cluster)
State3 = subset(Statemarker,Statemarker$cluster == 'State3')
 #table(unique(tmp2$Hugo_Symbol) %in% unique(State3$gene))
length(unique(State3$gene))
State3_pre = (196/276)*100
# neu.state3 = subset(tmp,tmp$Hugo_Symbol %in% State3$gene)

df = data.frame(Type = c('State1','State1','State2','State2','State3','State3'), 
                Mut = c('Mutation','Unmutated','Mutation','Unmutated','Mutation','Unmutated'),
                Freq = c(420,600-420, 573,857-573, 196,276-196),
                percent = c(State1_pre, 100-State1_pre, State2_pre, 100-State2_pre, State3_pre, 100-State3_pre),
                lable = c(paste0(round(State1_pre, digits = 2),'%'),paste0(round(100-State1_pre, digits = 2),'%'),
                          paste0(round(State2_pre, digits = 2),'%'),paste0(round(100-State2_pre, digits = 2),'%'),
                          paste0(round(State3_pre, digits = 2),'%'),paste0(round(100-State3_pre, digits = 2),'%')
                          ))

pvalue <- chisq.test(c(df$Freq,ncol=3))$p.value #卡方检验
library(plyr)
library(ggplot2)
"#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF""#FF7F7FFF""#8BCB3FFF"
ggplot(df,aes(Type,percent,fill=Mut))+
  geom_bar(stat="identity",position = position_stack())+
  # scale_color_igv()+
  scale_fill_manual(values = c("#DB423E","#42B540FF"),label=c("Mutation","Unmutated"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="State",y="Percent Weidght",
       fill="")+
  geom_text(aes(label=lable),vjust=1.5,size=6,color="black")+
  annotate(geom = "text",
           cex=6,
           x=2.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("Chi-Squared Test, P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = '堆积柱状图.pdf',width = 8,height = 6,path = './')

#获取共同突变基因------
rm(list = ls())
setwd("/home/datahup/syj/ESCC/SNP-TCGA-6.1/")
lusc.snv = read.delim(gzfile('/home/datahup/syj/NSCC/step8_1TCGA_snv/TCGA-ESCA.mutect2_snv.tsv.gz'))
load(file = "../../ESCC/拟时序-6/癌细胞_Statemarker.Rdata")
table(unique(Statemarker$gene) %in% unique(lusc.snv$gene))#用这个
743
load(file = "/home/datahup/syj/NSCC/step8_1TCGA_snv/snv.Rdata")
table(unique(Statemarker$gene) %in% unique(neu.state$Hugo_Symbol))
441 

Mut_genes <- unique(Statemarker$gene)[unique(Statemarker$gene) %in% unique(lusc.snv$gene)]
save(Mut_genes,file = "./突变基因.Rdata")










































