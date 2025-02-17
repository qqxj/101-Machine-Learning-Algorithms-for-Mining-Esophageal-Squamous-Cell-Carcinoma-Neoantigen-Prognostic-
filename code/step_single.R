

#2024.5.19(0.4)
# Find marker
rm(list = ls())
load(file = '/home/datahup/syj/NSCC/step1_1NoteCell/Step2_output.Rdata')
library(Seurat)
library(ggplot2)
library(Matrix)
library(stringr)
library(dplyr)
Idents(object = scRNA_harmony) <- "RNA_snn_res.0.4"
table(Idents(scRNA_harmony))
markers = FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

table(top10$cluster)
write.csv(markers,file = '/home/syj/NSCC/code/singleR/Step2_First_markers_0.4.csv')
write.csv(top10,file = '/home/syj/NSCC/code/singleR/Step2_First_markers_top10_0.4.csv')

# SingleR
sce = scRNA_harmony
library(Seurat)
library(celldex)
library(SingleR)
sce_for_SingleR <- GetAssayData(sce, slot="data")
sce@meta.data$seurat_clusters = sce@meta.data$RNA_snn_res.0.4
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
save(scRNA_harmony, file = "/home/syj/NSCC/code/singleR/Step2_output_0.4.Rdata")
write.csv(cellType,file = '/home/syj/NSCC/code/singleR/Step2_First_celltype_0.4.csv')


