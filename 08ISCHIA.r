library(tidyverse)
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)

predictions = c('/public/workspace/stu21230110/SPA_result/03cell2loc/C1/C1.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/C2/C2.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/C3/C3.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/P1/P1.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/P2/P2.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/P3/P3.spatial.deconvolution.csv')

sp_data = c('/public/workspace/stu21230110/SPA_result/07Mistry/C1/C1.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Mistry/C2/C2.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Mistry/C3/C3.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Mistry/P1/P1.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Mistry/P2/P2.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Mistry/P3/P3.seurat.rds')
samples = c('C1','C2','C3','P1','P2','P3')

sample = samples[6]

outdir = paste0('~/SPA_result/08ISCHIA/',sample, '/')

setwd(outdir)

sp_predict = read.csv(predictions[6],header = TRUE)

rownames(sp_predict) <- sp_predict[,1]

sp_predict <- sp_predict[,-1]

slide <- readRDS(sp_data[6])

# 创建一个新的 assay 来存储 predict 信息
predict_assay <- CreateAssayObject(t(as.matrix(sp_predict)))

# 将此 assay 添加到之前创建的 Seurat 对象中
slide[['predictions']] <- predict_assay

Assays(slide)

slide@meta.data=cbind(slide@meta.data,sp_predict)

saveRDS(slide,paste0(sample,'.spatial.rds'))




#ISCHIA
samples = c('C1','C2','C3','P1','P2','P3')

sample = samples[1]

outdir = paste0('~/SPA_result/08ISCHIA/',sample, '/')

sp_data = readRDS(paste0(outdir,sample,'.spatial.rds'))

setwd(outdir)



deconv.mat <- as.matrix(sp_data@meta.data[,14:24])

head(deconv.mat)

pdf("k.cluster.pdf")
Composition.cluster.k(deconv.mat, 20)
dev.off()

# Composition clustering of the deconvoluted spatial spots
sp_data <- Composition.cluster(sp_data, deconv.mat, 5)
#not good...Idents(sp_data) <- "CompositionCluster_CC"

table(sp_data$CompositionCluster_CC)
#> 
#> CC1 CC2 CC3 CC4 CC5 CC6 CC7 
#> 356 193 343 547 125 222 392 
pdf("CC_DimPlot.pdf")
SpatialDimPlot(sp_data, pt.size.factor = 1.8, group.by = c("CompositionCluster_CC")) + scale_fill_manual(values = c("cyan", "orange", "purple","green","yellow","blue"))
dev.off()
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.

CC = c("CC1","CC2","CC3","CC4","CC5","CC6","CC7")

pdf(paste0(CC[1],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[1], deconv.mat)
dev.off()

pdf(paste0(CC[2],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[2], deconv.mat)
dev.off()

pdf(paste0(CC[3],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[3], deconv.mat)
dev.off()

pdf(paste0(CC[4],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[4], deconv.mat)
dev.off()

pdf(paste0(CC[5],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[5], deconv.mat)
dev.off()

pdf(paste0(CC[6],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[6], deconv.mat)
dev.off()

pdf(paste0(CC[7],".enriched.celltype.pdf"))
Composition_cluster_enrichedCelltypes(sp_data, CC[7], deconv.mat)
dev.off()

sp_data.umap <- Composition_cluster_umap(sp_data, deconv.mat)
#> Plotting scatterpies for 2185 pixels with 20 cell-types...this could take a while if the dataset is large.
pdf("umap.cluster.pdf")
sp_data.umap$umap.cluster.gg
dev.off()

pdf("umap.deconv.pdf")
sp_data.umap$umap.deconv.gg
dev.off()

CC_cooccur = CC[3]

celltype.cooccur <- spatial.celltype.cooccurence(spatial.object = sp_data,
                                                      deconv.prob.mat = deconv.mat,
                                                      COI = CC_cooccur, prob.th= 0.4, Condition=unique(sp_data$orig.ident))
pdf(paste0(CC_cooccur,".cooccur.pdf"))
plot.celltype.cooccurence(celltype.cooccur)
dev.off()

saveRDS(sp_data,paste0(sample,".ischia.spatial.rds"))





sample = samples[1]

sp_data = readRDS(paste0(outdir,sample,'.ischia.spatial.rds'))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
all.LR.network <- cbind(lr_network[,c("from","to")], LR_Pairs=paste(lr_network$from, lr_network$to, sep = "_"))
all.LR.network.exp <- all.LR.network[which(all.LR.network$from %in% rownames(sp_data) & all.LR.network$to %in% rownames(sp_data)),]

# To reduce the computation time for this example, we randomly sample from the whole dataset of LR interactions
all.LR.network.exp <- sample_n(all.LR.network.exp,500)
all.LR.genes <- unique(c(all.LR.network.exp$from, all.LR.network.exp$to))
all.LR.genes.comm <- intersect(all.LR.genes, rownames(sp_data))
LR.pairs <- all.LR.network.exp$LR_Pairs
LR.pairs.AllCombos <- combn(all.LR.genes.comm, 2, paste0, collapse = "_")


CC6.Enriched.LRs <- Enriched.LRs(sp_data, c("CC6"), unique(sp_data$orig.ident), all.LR.genes.comm, LR.pairs, 1, 0.2)

CC7.Enriched.LRs <- Enriched.LRs(sp_data, c("CC7"), unique(sp_data$orig.ident), all.LR.genes.comm, LR.pairs, 1, 0.2)

CC6vsCC7.Enriched.LRs.Specific <- Diff.cooc.LRs(CC6.Enriched.LRs, CC7.Enriched.LRs, 0.05, 0.1)

pdf("CC6.ChordPlot.pdf")
ChordPlot.Enriched.LRs(CC6.Enriched.LRs$COI.enrcihed.LRs[1:20,])
dev.off()

pdf("CC7.ChordPlot.pdf")
ChordPlot.Enriched.LRs(CC7.Enriched.LRs$COI.enrcihed.LRs[1:20,])
dev.off()


pdf("CC6_CC7.SankeyPlot.pdf")
SankeyPlot.Diff.LRs(CC6.Enriched.LRs$COI.enrcihed.LRs, CC7.Enriched.LRs$COI.enrcihed.LRs)
dev.off()
