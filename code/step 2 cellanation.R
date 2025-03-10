rm(list = ls())
options(stringsAsFactors = F)
gc()

#设置工作路径
setwd("/home/datahup/syj/NSCC/step1_1NoteCell/")
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)

load(file = "/home/syj/NSCC/save/step 1/Step1_Finaldata.Rdata")

###细胞重命名（标识为NSCC_cell）
Idents(object = NSCC_list)
#levels(NSCC_list)
table(Idents(NSCC_list))
NSCC_list[["old.ident"]] <- Idents(object = NSCC_list)
Idents(object = NSCC_list) <- "NSCC_cell"

#再次查看细胞名字
Idents(object = NSCC_list)
levels(NSCC_list)
table(Idents(NSCC_list))
#保存重命名细胞
NSCC_list[["orig.ident"]] <- Idents(object = NSCC_list)


scRNA_harmony <- NSCC_list
scRNA_harmony
# QC
# 计算线粒体基因比例 线粒体是独立遗传的，不是染色体上基因控制的 默认<10%
Idents(scRNA_harmony) <- "NSCC_cell"
scRNA_harmony[["percent.mt"]] <- PercentageFeatureSet(scRNA_harmony, pattern = "^MT-")
table(scRNA_harmony[["percent.mt"]] < 10)
scRNA_harmony <- subset(scRNA_harmony, subset = nFeature_RNA > 1000 & 
                          nFeature_RNA < 5000 & percent.mt < 10 &
                           nCount_RNA > 500 & nCount_RNA < 35000)
# NSCLC.Integrate
#scRNA_harmony = subset(scRNA_harmony, subset = percent.mt < 10)
scRNA_harmony

#计算红血细胞基因比例 红细胞没有细胞核，没有转录组 默认<3%
rownames(scRNA_harmony)[grep("^HB[^(p)]", rownames(scRNA_harmony))]
scRNA_harmony <- PercentageFeatureSet(scRNA_harmony, "^HB[^(p)]", col.name = "percent_hb")
table(scRNA_harmony[["percent_hb"]] < 1)
scRNA_harmony = subset(scRNA_harmony, subset = percent_hb < 1)
scRNA_harmony

# 标准工作流
scRNA_harmony <- NormalizeData(scRNA_harmony)
scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA_harmony)
scRNA_harmony <- ScaleData(scRNA_harmony, features = all.genes)
scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(object = scRNA_harmony))

ElbowPlot(scRNA_harmony, ndims = 50)
ggsave(filename = "ElbowPlot.pdf",width = 12,height = 10,
       path = "./Fig/step 2/")

#保存Rdata
save(scRNA_harmony,
     file= "./save/step 2/scRNA_harmony_PCA.Rdata")

rm(list=ls())
load(file = "./save/step 2/scRNA_harmony_PCA.Rdata")

# Harmony整合
library(harmony)
Idents(scRNA_harmony) <- "patient"
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "patient")

#降维聚类 knn
scRNA_harmony = RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:25)###根据肘部图选择25###
scRNA_harmony = FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:25) ###选择25###
scRNA_harmony = FindClusters(object = scRNA_harmony, resolution = c(seq(0,1,by = 0.1))) #根据不同分辨率对细胞群聚类
library(clustree)
clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.") 
ggsave(filename = "resolution(0.3).pdf",width = 20,height = 14,
       path = "/home/syj/NSCC/Fig/step 2/")

save(scRNA_harmony,file = '/home/syj/NSCC/save/step 2/Step2_output_1.Rdata')

# 作图
rm(list = ls())
gc()
load(file = './save/step 2/Step2_output_1.Rdata')
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.4"
head(Idents(scRNA_harmony), 5)#查看前5个细胞的分类ID
table(Idents(scRNA_harmony))
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T, raster=T, label.box = T) 
plot1
ggsave(filename = "scRNA_harmony_0.4.pdf", plot = plot1, width = 8,height = 6,
       path = "/home/syj/NSCC/Fig/step 2/")

plot2 = DimPlot(scRNA_harmony, reduction = "umap", label=F, group.by = 'patient', raster=T) 
plot2
ggsave(filename = "scRNA_harmony_patient_0.4.pdf", plot = plot2, width = 14,height = 6,
       path = "/home/syj/NSCC/Fig/step 2/")

library(ggsci)
Idents(object = scRNA_harmony) <- "group"
head(Idents(scRNA_harmony), 5)#查看前5个细胞的分类ID
table(Idents(scRNA_harmony))
pal = pal_ucscgb(alpha = 0.8)(26)
plot3 = DimPlot(scRNA_harmony, reduction = "umap", label=F, raster=T, cols = pal) 
plot3
ggsave(filename = "scRNA_harmony_group_0.4.pdf", plot = plot3, width = 8,height = 6,
       path = "/home/syj/NSCC/Fig/step 2/")

# Find marker
rm(list = ls())
gc()
load(file = '/home/syj/NSCC/save/step 2/Step2_output_1.Rdata')
library(Seurat)
library(dplyr)
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.3"
table(Idents(scRNA_harmony))
markers = FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

table(top10$cluster)
write.csv(markers,file = './save/step 2/Step2_First_markers_1.csv')
write.csv(top10,file = './save/step 2/Step2_First_markers_top10_1.csv')

# SingleR
sce = scRNA_harmony
library(Seurat)
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.3
clusters=sce@meta.data$seurat_clusters

#加载注释包
load("/home/syj/NSCC/Data/Annotation packages/singleRref.Rdata")

#Blue.ref <- celldex::BlueprintEncodeData()
pred.Blue.ref <- SingleR(test = sce_for_SingleR, ref = Blue.ref, labels = Blue.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

#DICE.ref <- celldex::DatabaseImmuneCellExpressionData()
pred.DICE.ref <- SingleR(test = sce_for_SingleR, ref = DICE.ref, labels = DICE.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

#HPCA.ref <- celldex::HumanPrimaryCellAtlasData()
pred.HPCA.ref <- SingleR(test = sce_for_SingleR, ref = HPCA.ref, labels = HPCA.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

#Mona.ref <- celldex::MonacoImmuneData()
pred.Mona.ref <- SingleR(test = sce_for_SingleR, ref = Mona.ref, labels = Mona.ref$label.fine ,
                         method = "cluster", clusters = clusters, 
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")

#Nover.ref <- celldex::NovershternHematopoieticData()
pred.Nover.ref <- SingleR(test = sce_for_SingleR, ref = Nover.ref, labels = Nover.ref$label.fine ,
                          method = "cluster", clusters = clusters, 
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellType=data.frame(ClusterID=levels(sce@meta.data$seurat_clusters),
                    Blue=pred.Blue.ref$labels,
                    DICE=pred.DICE.ref$labels,
                    HPCA=pred.HPCA.ref$labels,
                    Mona=pred.Mona.ref$labels,
                    Nover=pred.Nover.ref$labels )
head(cellType)
sce@meta.data$singleR_Blue=cellType[match(clusters,cellType$ClusterID),'Blue']
sce@meta.data$singleR_DICE=cellType[match(clusters,cellType$ClusterID),'DICE']
sce@meta.data$singleR_HPCA=cellType[match(clusters,cellType$ClusterID),'HPCA']
sce@meta.data$singleR_Nover=cellType[match(clusters,cellType$ClusterID),'Nover']
sce@meta.data$singleR_Mona=cellType[match(clusters,cellType$ClusterID),'Mona']

scRNA_harmony = sce

#保存scRNA_harmony和机器注释结果
save(scRNA_harmony, file = "./save/step 2/Step2_output_1.Rdata")
write.csv(cellType,file = './save/step 2/Step2_First_celltype_1.csv')

# check marker在每个簇中的表达情况
library(Seurat)
library(ggplot2)
library(ggplot2)
load(file = '/home/syj/NSCC/save/step 2/Step2_output_1.Rdata')
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.4"
table(Idents(scRNA_harmony))
top10 <- read.csv('./code/singleR/Step2_First_markers_top10_0.4.csv',row.names = 1)

Gene <- subset(top10, top10$cluster==26) # 依次从0开始检查每个簇的top10marker基因
p_all_markers=DotPlot(scRNA_harmony,
                      features = Gene$gene,
                      cols = c("lightgrey", "purple"),
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
Gene$gene

# Gene <- c('FAP','THY1') # 检查感兴趣的基因
# p_all_markers1=DotPlot(scRNA_harmony,
#                        features = Gene,
#                        cols = c("lightgrey", "purple"),
#                        scale = T,assay='RNA' )+
#   theme(axis.text.x=element_text(angle=45,hjust = 1))
# p_all_markers1

#热图(先看一下结果）
library(pheatmap)
library(ggplot2)
load(file = '/home/syj/NSCC/save/step 2/Step2_output.Rdata')
table(Idents(scRNA_harmony))
top10 <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_celltype_0.4.csv')#,row.names = 1
# #1
# DoHeatmap(Lymph, features = top10$gene) + NoLegend()
# ggsave(filename = "Top10-MarkerGene.pdf", plot = plot, width = 12,height = 8,
#        path = "/home/syj/NSCC/Fig/step 2/")
#2
library(ComplexHeatmap)
library(ggunchull)
library(scRNAtoolVis)
library(ggplot2)
#show_col(pal_d3('category20')(20))
library(RColorBrewer)
pal <- brewer.pal(20, "Set3")
#pal = pal_d3('category20')(20)
library(viridis)
pal <- viridis(30, option = "D", alpha = 0.9)
pdf(file = './Fig/step 2/scRNA_harmony_pheatmap_0.4.pdf',width = 25,height = 25)
AverageHeatmap(object = scRNA_harmony,
               markerGene = top10$gene,
               column_names_rot = 0, 
               htCol = c("blue", "yellow", "red")
               )
dev.off()

#调整top_10行的排列顺序，重新画上方热图
table(top10$cluster)
custom_order  <- factor(top10$cluster, levels = c("0", "2", "1", "3","4","9","5","6","7","8"))  # 按照需要定义排序顺序
top10_1 <- top10[order(custom_order), ]
top10_1 <- top10_1[!(top10_1$gene %in% c("LTB", "IGFBP7","CALD1","HLA−DRA","HLA−DQB1","HLA−DQA1","SPARCL1")), ]

top10 <- top10_1

#人工注释
rm(list = ls())
options(stringsAsFactors = F)
gc()
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)
library(dplyr)
library(clustree)

setwd("/home/syj/NSCC/")

#load(file = "./save/step 2/Step2_output_1.Rdata")
load(file = "./code/singleR/Step2_output_0.4.Rdata")

#marker基因数据(csv文件中有汉字是无法读入的)
#markers <- read.csv('./code/singleR/Step2_First_markers_0.4.csv')
#top10 <- read.csv('./code/singleR/Step2_First_markers_top10.csv')
#自动注释数据
#cellType <- read.csv('./code/singleR/Step2_First_celltype_0.4.csv')#进行人工注释

###cellmarker###
# 结合Cellmarker数据库和SingleR结果，半监督注释细胞
celltype <- read.csv(file = './code/singleR/Step2_First_celltype_0.4.csv',header = T)
celltype$X
new.cluster.ids = celltype$X
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.4"#选用0.1的细胞分簇
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
table(Idents(scRNA_harmony))
scRNA_harmony@meta.data$FirstAnnotation = Idents(scRNA_harmony)

#plot
library(ggsci)
pal = pal_ucscgb(alpha = 0.8)(26)
pal = pal[-c(7,8)]
plot = DimPlot(scRNA_harmony, reduction = "umap", label=T, label.box = T, cols = pal, raster = T) + NoLegend()
plot
ggsave(filename = "scRNA_harmony_First_0.4.pdf", plot = plot, width = 12,height = 8,
       path = "/home/syj/NSCC/code/singleR/")

#保存
save(scRNA_harmony, file = "/home/syj/NSCC/code/singleR/Step2_output.Rdata")


#可视化--------
load(file = "/home/datahup/syj/NSCC/step1_1NoteCell/Step2_output.Rdata")
cell_cluster <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_celltype_0.4.csv')#,row.names = 1
top10 <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_markers_top10_0.4.csv',row.names = 1)

#热图(先看一下结果）
library(pheatmap)
library(ggplot2)
#需要数据单细胞数据、gene名
table(Idents(scRNA_harmony))
gene <- c("POSTN","COL1A1","COL3A1",
          "NKG7","GNLY","CCL5",
          "IGHG1","IGHG4","IGKC",
          "HLA-DPB1","HLA-DPA1","HLA-DQB1",
          "KRT6B","KRT16","KRT17",
          "ACTA2","TAGLN","CALD1",
          "IL7R","ICOS",
          "TPSAB1","TPSB2","CPA3",
          "MS4A1","CD79A","CD19",
          "EPCAM","KRT14","KRT5","KRT15",#"KLF5",#上皮
          "PLVAP","SELE","RAMP2",
          "TYROBP","CD14","CD68")
# #1
# DoHeatmap(Lymph, features = top10$gene) + NoLegend()
# ggsave(filename = "Top10-MarkerGene.pdf", plot = plot, width = 12,height = 8,
#        path = "/home/syj/NSCC/Fig/step 2/")
#2
library(ComplexHeatmap)
library(ggunchull)
library(scRNAtoolVis)
library(ggplot2)
#show_col(pal_d3('category20')(20))
library(RColorBrewer)
pal <- brewer.pal(20, "Set3")
#pal = pal_d3('category20')(20)
library(viridis)
pal <- viridis(30, option = "D", alpha = 0.9)
pdf(file = './code/singleR/scRNA_harmony_pheatmap_0.4.pdf',width = 25,height = 25)
AverageHeatmap(object = scRNA_harmony,
               markerGene = gene,
               column_names_rot = 0, 
               htCol = c("blue", "yellow", "red")
)
dev.off()

'#330066','#336699','#66CC66','#FFCC33'
# marker的气泡图
library(Seurat)
library(ggplot2)
library(ggplot2)
#需要数据：单细胞数据集、gene名
Idents(object = scRNA_harmony) <- scRNA_harmony@meta.data[["FirstAnnotation"]]
table(Idents(scRNA_harmony))

gene
pdf(file = './code/singleR/scRNA_harmony_DotPlot_0.4.pdf',width = 15,height = 15)
DotPlot(scRNA_harmony,
                      features = gene,
                      cols = c("lightgrey", "purple"),
                      scale = T,
                      assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
dev.off()

p <- DotPlot(pbmc, features = unique(top3pbmc.markers$gene) ,
             assay='RNA' ) + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
p

#美化图片----
load(file = "/home/datahup/syj/NSCC/step1_1NoteCell/Step2_output.Rdata")
cell_cluster <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_celltype_0.4.csv')#,row.names = 1
top10 <- read.csv('/home/datahup/syj/NSCC/step1_1NoteCell/Step2_First_markers_0.4.csv',row.names = 1)

gene_cell_pairs <- data.frame(
  gene = c("COL1A1", "COL3A1", "POSTN",
           "NKG7","CCL5", "GNLY",  
           "IGHG1","IGHG4","IGKC",
           "HLA-DPB1", "HLA-DPA1", "HLA-DQB1", 
           "KRT17","KRT6B", "KRT16",  
           "CALD1","ACTA2", "TAGLN",  
           "IL7R", "ICOS", 
           "TPSAB1", "TPSB2", "CPA3", 
           "MS4A1", "CD79A", "CD19", 
           "KRT5","EPCAM","KRT14","KRT15",
           "PLVAP", "RAMP2","SELE",  
           "TYROBP", "CD14", "CD68"),
  cell_type = c(rep("Fibroblast", 3),
                rep("NK cell", 3),
                rep("Plasma cell", 3),
                rep("Dendritic cell", 3),
                rep("Basal cell", 3),
                rep("Smooth muscle cell",3),
                rep("T cell", 2),
                rep("Mast cell", 3),
                rep("B cell", 3),
                rep("Epithelial cell", 4),
                rep("Endothelial cell", 3),
                rep("Macrophage", 3)))
#定义顺序
gene_cell_pairs$cell_type <- factor(gene_cell_pairs$cell_type, 
                                    levels = c("Fibroblast", "NK cell", "Plasma cell",
                                               "Dendritic cell", "Basal cell", "Smooth muscle cell",
                                               "T cell", "Mast cell", "B cell",
                                               "Epithelial cell", "Endothelial cell", 
                                               "Macrophage"))
#可视化
p1 <- DotPlot(scRNA_harmony,
              features = split(gene_cell_pairs$gene, gene_cell_pairs$cell_type),
              #cols = c("#ffffff", "#448444")
) +
  RotatedAxis() + # 来自Seurat
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
  )+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
p1
ggsave("细胞注释气泡图.png", plot = p1, width = 22, height = 11)
ggsave("细胞注释气泡图.pdf", plot = p1, width = 22, height = 11)



















