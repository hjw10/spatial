suppressMessages({
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
})


sc_rds = '/public/workspace/stu21230110/scRNA/06_scRNA_harmony_singler202400810.Rdata'
spatial_rds = '/public/workspace/stu21230110/SPA_result/01CCA_Seurat/combined.seurat.rds'
outdir = '/public/workspace/stu21230110/SPA_result/06MIA/result'
normalization_method = 'SCT'
sc_min = 0.3
sc_logfc = 0.25
region = '/public/workspace/stu21230110/SPA_result/celltype.csv'
sp.pct = 0.3
sp.logfc = 0.25


if (!file.exists(outdir)) {dir.create(outdir,recursive = TRUE)}

load(sc_rds)

cortex_sc <- scRNA_harmony_singler

if (normalization_method == 'SCT'){

cortex_sc <- SCTransform(cortex_sc, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:20)}else{

cortex_sc <- NormalizeData(cortex_sc, normalization.method = "LogNormalize", scale.factor = 10000)

cortex_sc <- FindVariableFeatures(cortex_sc, selection.method = "vst", nfeatures = 2000)

cortex_sc <- ScaleData(cortex_sc, features = rownames(cortex_sc))

cortex_sc <- RunPCA(cortex_sc, features = VariableFeatures(object = cortex_sc))

cortex_sc <- RunUMAP(cortex_sc, dims = 1:20)

}

DefaultAssay(cortex_sc) = 'RNA'

sc.markers <- FindAllMarkers(cortex_sc, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)

sc.markers$d = sc.markers$pct.1 - sc.markers$pct.2

sc.main.marker = subset(sc.markers , avg_log2FC > 1 & p_val_adj < 0.01 & d > 0.3)

sc.main.marker = sc.main.marker %>% arrange(cluster,desc(avg_log2FC))

sc.main.marker = as.data.frame(sc.main.marker)

sc.main.marker$cluster = paste('sc',sc.main.marker$cluster,sep = '_')
#####空间转录组处理

cortex_sp = readRDS(spatial_rds)

cortex_sp <- SCTransform(cortex_sp, verbose = FALSE,assay = "Spatial")

cortex_sp <- RunPCA(cortex_sp, assay = "SCT", verbose = FALSE)

cortex_sp <- FindNeighbors(cortex_sp, reduction = "pca", dims = 1:20)

cortex_sp <- FindClusters(cortex_sp, verbose = FALSE,resolution = 0.5)

cortex_sp <- RunUMAP(cortex_sp, reduction = "pca", dims = 1:20)


region_marker=FindAllMarkers(cortex_sp,logfc.threshold = 0.25,only.pos = T,min.pct = 0.25)

region_marker$d=region_marker$pct.1 - region_marker$pct.2

region_main_marker = subset(region_marker,avg_log2FC > 0.6 & p_val_adj < 0.01 & d > 0.3)

region_main_marker = region_main_marker %>% arrange(cluster,desc(avg_log2FC))

region_main_marker = as.data.frame(region_main_marker)

region_main_marker$cluster = paste('spatial',region_main_marker$cluster,sep = '_')

#####MIA

region_specific = region_main_marker[,c("cluster","gene")]

colnames(region_specific)[1]="region"

celltype_specific = sc.main.marker[,c("cluster","gene")]

colnames(celltype_specific)[1]="celltype"

source("MIA.R")

Result = zhao_MIA(region_specific,celltype_specific,'adenomyosis',outdir)

saveRDS(region_specific,'region_specific.rds')

saveRDS(celltype_specific,'celltype_specific.rds')

saveRDS(cortex_sc,"cortex_sc.rds")

saveRDS(cortex_sp,"cortex_sp.rds")
