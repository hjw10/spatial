##构建单样本数据
input_files = c("/public/workspace/stu21230110/SPA/C1",
                "/public/workspace/stu21230110/SPA/C2",
                "/public/workspace/stu21230110/SPA/C3",
                "/public/workspace/stu21230110/SPA/P1",
                "/public/workspace/stu21230110/SPA/P2",
                "/public/workspace/stu21230110/SPA/P3")
genename_file = c("/public/workspace/stu21230110/SPA/C1/filtered_feature_bc_matrix/features.tsv.gz",
                  "/public/workspace/stu21230110/SPA/C2/filtered_feature_bc_matrix/features.tsv.gz",
                  "/public/workspace/stu21230110/SPA/C3/filtered_feature_bc_matrix/features.tsv.gz",
                  "/public/workspace/stu21230110/SPA/P1/filtered_feature_bc_matrix/features.tsv.gz",
                  "/public/workspace/stu21230110/SPA/P2/filtered_feature_bc_matrix/features.tsv.gz",
                  "/public/workspace/stu21230110/SPA/P3/filtered_feature_bc_matrix/features.tsv.gz")
                  
sample_names = c("C1","C2","C3","P1","P2","P3")

image_files = file.path(input_files,"spatial")

sample = sample_names[6]

outdir = paste0("/public/workspace/stu21230110/SPA_result/07Misty/",sample)

input_format = "10x"

assay = "Spatial"
is_harmony = TRUE
use_colors = TRUE
## qc
#min_features = 200
low_feature = 200
high_feature = NULL
low_mt_percent = 0
low_hb_percent = 0
normalize_method = "SCT"
## pca
selection_method = "vst"
pc_dim = 50
x_low_cutoff = 0.125
x_high_cutoff = 5
dispersion_cutoff = 1
nfeatures = 3000
active.assay = assay
pc_number = NULL
## cluster
cluster.resolution = 0.5
annoy_metric = "euclidean"
## diff
logfc.diff = 0.5
padjust.diff = 0.05
logfc.marker = log(2)
padjust.marker = 0.05
test_method = "wilcox"
min_pct = 0.1
crop = FALSE

high_mt_percent = 0.99
high_hb_percent = 0.99
low_res = 70
mid_res = 150
high_res = 300
mt.genes0 = NULL
hb.genes0 = NULL

high_mt_percent = 0.99
high_hb_percent = 0.99


defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')

low_res = 70
mid_res = 150
high_res = 300

library(foreach)
library(grid)
library(dplyr)
library(RColorBrewer)
library(Seurat)
library(foreach)
library(reshape2)
library(ggplot2)
library(psych)
library(pheatmap)
library(reticulate)
#library(htmlwidgets)
library(plotly)
library(reticulate)
library(scikit-learn)

source("/public/workspace/stu21230110/SPA_result/utils.R")

setwd(outdir)
qc.dir = file.path(outdir, "qc")
pca.dir = file.path(outdir, "pca")
cluster.dir = file.path(outdir, "clustering")
tsne.dir = file.path(outdir, "tsne")
umap.dir = file.path(outdir, "umap")
diff.dir = file.path(outdir, "diffexp")
diff.src.dir = file.path(outdir, "diffexp","src")
diff.marker.dir = file.path(outdir, "diffexp","markers")
dir.create(qc.dir, recursive = TRUE)
dir.create(pca.dir, recursive = TRUE)
dir.create(cluster.dir, recursive = TRUE)
dir.create(tsne.dir, recursive = TRUE)
dir.create(umap.dir, recursive = TRUE)
dir.create(diff.dir, recursive = TRUE)
dir.create(diff.src.dir, recursive = TRUE)
dir.create(diff.marker.dir, recursive = TRUE)

## QC / meta data
## MT gene, HB genes; seruat_obj@meta.data
## gene name
GENENAME = Read_GeneName(input.file = input_file[6], input.format = input_format, genename.file = genename_file[6])
has_genename = TRUE
if(is.null(GENENAME)){
	has_genename = FALSE
}

SCrna.raw = Read_SC(
	input.file = input_files[6], 
	input.format = input_format, 
	image.file = image_files[6], 
	assay = assay, 
	GENENAME = GENENAME,
	project = sample, 
	min.cells = 200
)
SCrna.raw@meta.data$sample = sample

feature.names = paste(c("nCount_","nFeature_"),assay,sep="")
## MT
mt.genes0 = NULL

mt.genes = find_target_genes(SCrna.raw, mt.genes0, GENENAME, use.name=has_genename, mt=T)
if(length(mt.genes) != 0){
	feature.names = c(feature.names, "percent_MT")
	SCrna.raw[['percent_MT']] = PercentageFeatureSet(SCrna.raw, features = mt.genes)
	if(assay == "Spatial"){
		f1 = SpatialFeaturePlot(SCrna.raw, features = "percent_MT", crop=FALSE, pt.size.factor=1.2) + 
			theme(legend.position = "right") + coord_cartesian()
		dual.plot(f1, file.path(qc.dir, paste(sample,'.percent_MT.Spatial',sep="")), res=mid_res)
	}
}
## HB 
hb.genes0 = NULL

hb.genes = find_target_genes(SCrna.raw, hb.genes0, GENENAME, use.name=has_genename, mt=F)
if(length(hb.genes) != 0){
	feature.names = c(feature.names, "percent_HB")
	SCrna.raw[['percent_HB']] = PercentageFeatureSet(SCrna.raw, features = hb.genes)
}
print(head(SCrna.raw@meta.data))

SCrna = SCrna.raw
if(assay == "Spatial"){
	SCrna.raw@meta.data$filter = factor(colnames(SCrna.raw) %in% colnames(SCrna) + 0, levels=c(0,1))
	f = SpatialDimPlot(SCrna.raw, group.by="filter", label = FALSE, label.size = 3, crop=FALSE, pt.size.factor=1.2) + 
		coord_cartesian() + scale_fill_manual(values = c("0"="gray","1"="red"))
	dual.plot(f, file.path(outdir, 'qc', paste(sample, ".FilterCells.Spatial",sep="")),w=7,h=7,res=200)
}
rm(SCrna.raw)


SCrna = SCTransform(SCrna, assay = assay, verbose = FALSE, variable.features.n = nfeatures)
sel.features = VariableFeatures(object = SCrna)
active.assay = "SCT"


# pca
SCrna = RunPCA(object = SCrna, assay = active.assay, do.print = F, npcs = pc_dim)
## Determine the ‘dimensionality’ 
pc.stdev = data.frame(PC = 1:length(SCrna@reductions$pca@stdev), stdev = SCrna@reductions$pca@stdev)
pc.fit = nls(stdev ~ a*PC^b, data = pc.stdev, start = list(a=10, b= -0.5),trace = T)
if(!is.null(pc_number)){
	pc.num = min(20, pc_number) ## less than 20
	pc.num = max(5, pc.num)    ## greater than 5
}else{
	pc.num = determine_PCnum(fitted(pc.fit))
}
cat(paste("pc.number: ",pc.num,"\n"))

# cluster
cat("cluster ...\n")
SCrna = FindNeighbors(SCrna, reduction = "pca", dims = 1:pc.num, annoy.metric = annoy_metric)
SCrna = FindClusters(SCrna, resolution = cluster.resolution)
cluster.ids = levels(SCrna@meta.data$seurat_clusters)

# tSNE / UMAP
cat("tsne / umap ...\n")
SCrna = RunTSNE(object = SCrna, reduction = "pca", dims.use = 1:pc.num, do.fast = TRUE)
SCrna = RunUMAP(object = SCrna, reduction = "pca", dims = 1:pc.num, metric = "correlation")

saveRDS(SCrna, file=file.path(outdir,paste(sample,".seurat.rds",sep="")))

.libPaths(c('~/tools/seurat4/'))
library(tidyverse)
library(Seurat)
library(mistyR)
library(SeuratDisk)
source("misty_utilities.R")
library(hdf5r)

predictions = c('/public/workspace/stu21230110/SPA_result/03cell2loc/C1/C1.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/C2/C2.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/C3/C3.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/P1/P1.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/P2/P2.spatial.deconvolution.csv',
                '/public/workspace/stu21230110/SPA_result/03cell2loc/P3/P3.spatial.deconvolution.csv')

sp_data = c('/public/workspace/stu21230110/SPA_result/07Misty/C1/C1.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Misty/C2/C2.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Misty/C3/C3.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Misty/P1/P1.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Misty/P2/P2.seurat.rds',
            '/public/workspace/stu21230110/SPA_result/07Misty/P3/P3.seurat.rds')
samples = c('C1','C2','C3','P1','P2','P3')

sample = samples[6]

outdir = paste0('~/SPA_result/07Misty/', sample, '/Misty_results/')

#future::plan(future::multisession)

run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = outdir) {  ###输出目录大家自己制定
   
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 2,
                      "para" = 5)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}

sp_predict = read.csv(predictions[6],header = TRUE)

rownames(sp_predict) <- sp_predict[,1]

sp_predict <- sp_predict[,-1]

slide <- readRDS(sp_data[6])

slide_id = sample
# 创建一个新的 assay 来存储 predict 信息
predict_assay <- CreateAssayObject(t(as.matrix(sp_predict)))

# 将此 assay 添加到之前创建的 Seurat 对象中
slide[['predictions']] <- predict_assay

# 验证对象现在包含多个 assays
Assays(slide)

slide@meta.data=cbind(slide@meta.data,sp_predict)

DefaultAssay(slide) <- 'predictions'

useful_features <- rownames(slide)   ####也可以自我设定感兴趣的细胞类型

useful_features <- useful_features[! useful_features %in% "prolif"]

mout <- run_colocalization(slide = slide,
                     useful_features = useful_features,
                     out_label = slide_id,
                     assay = "predictions",
                     misty_out_alias = outdir)


misty_res_slide <- collect_results(mout)
  
  plot_folder <- paste0(mout, "/plots")
  
  system(paste0("mkdir ", plot_folder))
  
  pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
  
  mistyR::plot_improvement_stats(misty_res_slide)
  mistyR::plot_view_contributions(misty_res_slide)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "juxta_2", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "juxta_2", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "para_5", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "para_5", cutoff = 0.5)
  
  dev.off()