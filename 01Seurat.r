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
library(cowplot)
#library(htmlwidgets)
library(plotly)


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
image_files = file.path(input_files,"spatial")
input_format = "10x"
outdir = "./01Seurat/"
sample_names = c("C1","C2","C3","P1","P2","P3")
min.cells = 3
assay = "Spatial"
sample="combined"
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

source("~/SPA_result/utils.R")
library(RColorBrewer)
sample_cols = unique(c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")))
defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', sample_cols)
names(defined_cols) = 1:length(defined_cols) - 1


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


## gene name
has_genename = TRUE

SC.list = list()
feature.names = paste(c("nCount_","nFeature_"),assay,sep="")
HVGs = c()
image_names = sample_names ##paste("image",1:length(sample_names),sep="")
for(i in 1:length(sample_names)){
	GENENAME = Read_GeneName(input.file = input_files[i], input.format = input_format, genename.file = genename_file[i])
	sctmp = Read_SC(
		input.file = input_files[i],
		input.format = input_format,
		image.file = image_files[i],
		image.name = image_names[i],
		assay = assay,
		GENENAME = GENENAME,
		project = sample_names[i],
		min.cells = 200
	)
	sctmp@meta.data$sample = sample_names[i]
	## MT
	feature.names.tmp = paste(c("nCount_","nFeature_"),assay,sep="")
	mt.genes =mt.genes = grep("^MT-", rownames(sctmp), ignore.case = TRUE, value = T)
	if(length(mt.genes) != 0){
		feature.names.tmp = c(feature.names.tmp, "percent_MT")
		feature.names = unique(c(feature.names, "percent_MT"))
		sctmp[['percent_MT']] = PercentageFeatureSet(sctmp, features = mt.genes)
	}
	## HB 
	hb.genes =  rownames(sctmp)[grep("^Hb[^(p)]", rownames(sctmp),ignore.case = T)]
	if(length(hb.genes) != 0){
		feature.names.tmp = c(feature.names.tmp, "percent_HB")
		feature.names = unique(c(feature.names, "percent_HB"))
		sctmp[['percent_HB']] = PercentageFeatureSet(sctmp, features = hb.genes)
	}
	## filtration
		sctmp = FilterSeuratObject(sctmp,
		assay = assay,
		low.feature = low_feature,
		high.feature = high_feature,
		mt.genes = mt.genes,
		high.mt = high_mt_percent,
		hb.genes = hb.genes,
		high.hb = high_hb_percent
	)
	cat("normalize ... \n")
	if(normalize_method != "SCT"){
		scale_factor = median(sctmp@meta.data[,feature.names.tmp[1]])
		print(paste("scale_factor = ", scale_factor))
		sctmp = NormalizeData(sctmp, assay=assay, normalization.method=normalize_method, scale.factor=scale_factor)
		sctmp = FindVariableFeatures(object = sctmp,
			selection.method = selection_method,
			assay = assay,
			mean.function = ExpMean,
			dispersion.function = LogVMR,
			mean.cutoff = c(x_low_cutoff, x_high_cutoff),
			dispersion.cutoff = c(dispersion_cutoff, Inf),
			nfeatures = nfeatures
		)
	}else{
		sctmp = SCTransform(sctmp, assay = assay, verbose = FALSE, variable.features.n = 3000)
	}
	HVGs = c(HVGs, VariableFeatures(object = sctmp))
	SC.list[[i]] = sctmp
}
if(normalize_method == "SCT"){
	active.assay = "SCT"
}

SC.anchors <- FindIntegrationAnchors(object.list = SC.list, dims = 1:20)
SCrna <- IntegrateData(anchorset = SC.anchors, dims = 1:20)
active.assay = "integrated"
active.assay = "SCT"
DefaultAssay(SCrna) = active.assay
#跑harmony用SCT
#跑CCA用integrated 
#SCrna$integrated = NULL

VariableFeatures(SCrna) = unique(HVGs)
SCrna = ScaleData(object = SCrna, assay = active.assay)
## revised 11/03/2020
#SCrna = SCTransform(SCTransform, verbose = FALSE, variable.features.n = nfeatures)

# pca
cat("pca ...\n")
SCrna = RunPCA(object = SCrna, assay = active.assay, do.print = F, npcs = 100)


## Determine the ‘dimensionality’ 
pc.stdev = data.frame(PC = 1:length(SCrna@reductions$pca@stdev), stdev = SCrna@reductions$pca@stdev)
pc.fit = nls(stdev ~ a*PC^b, data = pc.stdev, start = list(a=10, b= -0.5),trace = T)
if(!is.null(pc_number)){
	pc.num = min(20, pc_number) ## less than 20
	pc.num = max(5, pc.num)    ## greater than 5
}else{
	pc.num = determine_PCnum(fitted(pc.fit))
}

library(harmony)
SCrna = RunHarmony(SCrna, group.by.vars = "orig.ident")
reduction = "harmony"
reduction = "pca"

# cluster
cat("cluster ...\n")
SCrna = FindNeighbors(SCrna,  reduction = "harmony", dims = 1:20, annoy.metric = annoy_metric)
SCrna = FindClusters(SCrna, resolution = 0.5)
cluster.ids = levels(SCrna@meta.data$seurat_clusters)

setwd("~/SPA_result/01seurat/")
pdf("CCA.pdf")
SpatialDimPlot(SCrna,label=T,label.size=3,combine=FALSE, crop=FALSE, pt.size.factor = 1.2)
dev.off()

pdf("Harmony.pdf")
SpatialDimPlot(SCrna,label=T,label.size=3,combine=FALSE, crop=FALSE, pt.size.factor = 1.2)
dev.off()

# tSNE / UMAP
cat("tsne / umap ...\n")
SCrna = RunTSNE(object = SCrna, reduction = "harmony", dims.use = 1:pc.num, do.fast = TRUE)
SCrna = RunUMAP(object = SCrna, reduction = "harmony", dims = 1:pc.num, umap.method = "umap-learn",metric = "correlation")

cell_type = read.csv("~/SPA_result/all_celltype.csv")
rownames(cell_type) = cell_type$X
cell_type = cell_type[,-1]
integrated_compositions <- cell_type[intersect(rownames(cell_type),rownames(SCrna@meta.data)),]

predict_assay <- CreateAssayObject(t(as.matrix(integrated_compositions)))
SCrna[["predictions"]] <- predict_assay
DefaultAssay(SCrna) <- "predictions"
SCrna@meta.data=cbind(SCrna@meta.data,integrated_compositions)


saveRDS(SCrna,"~/SPA_result/myRDS/filter_combined_spatial.rds")

##########################################################
###       figures and tables
##########################################################

write.xls(SCrna@meta.data[,c("orig.ident",feature.names)], file.path(qc.dir, paste(sample,".meta_data.xls",sep="")), first.colname="barcode")

f = VlnPlot(SCrna, features = feature.names, ncol = length(feature.names), group.by="orig.ident", col=sample_cols)
dual.plot(f, file.path(qc.dir, paste(sample,'.QC_metrics.violin',sep="")), w=4*length(feature.names), h=6, res=mid_res)

ncol.img = 3
if(length(sample_names)==1){
	ncol.img = 1
}else if (length(sample_names)<=4){
	 ncol.img = 2
}
nrow.img = ceiling(length(sample_names)/ncol.img)


if(assay == "Spatial"){
	f1 = SpatialFeaturePlot(SCrna, features = "nCount_Spatial",combine=FALSE, crop=crop, pt.size.factor = 1.2)
	for(x in 1:length(f1)){
		f1[[x]] = f1[[x]] + coord_cartesian()
	}
	p1 = CombinePlots(plots = f1, ncol = ncol.img)
	dual.plot(p1, file.path(qc.dir, paste(sample,'.nCount.Spatial',sep="")), w=5*ncol.img, h=6*nrow.img, res=mid_res)

	f2 = SpatialFeaturePlot(SCrna, features = "nFeature_Spatial",combine=FALSE, crop=crop, pt.size.factor = 1.2)
	for(x in 1:length(f2)){
		f2[[x]] = f2[[x]] + coord_cartesian()
	}
	p2 = CombinePlots(plots = f2, ncol = ncol.img)
	dual.plot(p2, file.path(qc.dir, paste(sample,'.nFeature.Spatial',sep="")), w=5*ncol.img, h=6*nrow.img, res=mid_res)

	if("percent_MT" %in% colnames(SCrna@meta.data)){
		f3 = SpatialFeaturePlot(SCrna, features = "percent_MT",combine=FALSE, crop=crop, pt.size.factor = 1.2)
		for(x in 1:length(f3)){
			f3[[x]] = f3[[x]] + coord_cartesian()
		}
		p3 = CombinePlots(plots = f3, ncol = ncol.img)
		dual.plot(p3, file.path(qc.dir, paste(sample,'.percent_MT.Spatial',sep="")), w=5*ncol.img, h=6*nrow.img, res=mid_res)
	}

	if("percent_HB" %in% colnames(SCrna@meta.data)){
		f4 = SpatialFeaturePlot(SCrna, features = "percent_HB", combine=FALSE, crop=crop, pt.size.factor = 1.2)
		for(x in 1:length(f4)){
			f4[[x]] = f4[[x]] + coord_cartesian()
		}
		p4 = CombinePlots(plots = f4, ncol = ncol.img)
		dual.plot(f4, file.path(qc.dir, paste(sample,'.percent_HB.Spatial',sep="")), w=5*ncol.img, h=6*nrow.img, res=mid_res)
	}
}

## FeatureScatter
scatter_plots = list()
for (i in 2:length(feature.names)){
	scatter_plots[[i-1]] = FeatureScatter(object = SCrna, feature1 = feature.names[1], feature2 = feature.names[i], group.by="orig.ident")
}
nscatter = length(scatter_plots)

f = CombinePlots(plots = scatter_plots, ncol = nscatter, legend = "none")
dual.plot(f, file.path(qc.dir,paste(sample,".QC_metrics.scatter",sep="")), w=nscatter*5, h=5, res=mid_res)

#counts = as(GetAssayData(SCrna,slot="counts"),"matrix")
#counts = data.frame(Gene=map.genename(rownames(counts),GENENAME,3,1), counts, check.names=F)
#write.table(counts,file=file.path(outdir,paste(sample,".genebar.tsv",sep="")),quote=F,sep="\t",row.names=F)

exprs = round(expm1(as(GetAssayData(SCrna,slot="data"),"matrix")),6)
exprs = data.frame(Gene=map.genename(rownames(exprs),GENENAME,3,1), GeneName=map.genename(rownames(exprs),GENENAME,3,2), exprs, check.names=F)
write.table(exprs,file=file.path(outdir,paste(sample,".gene_expression.xls",sep="")),quote=F,sep="\t",row.names=F)

avg_exp = AverageExpression(SCrna, assay = SCrna@active.assay, verbose=T)[[1]]
colnames(avg_exp) = paste("cluster", colnames(avg_exp), sep="")
avg_exp = data.frame(Gene=map.genename(rownames(avg_exp),GENENAME,3,1), GeneName=map.genename(rownames(avg_exp),GENENAME,3,2), avg_exp, check.names=F)
write.table(avg_exp,file=file.path(outdir,paste(sample,".gene_expression.clusters.xls",sep="")),quote=F,sep="\t",row.names=F)


if(FALSE){
	plot1 = VariableFeaturePlot(SCrna)
	plot2 = LabelPoints(plot = plot1, points = HVGs[1:10], repel = TRUE)
	#f = CombinePlots(plots = list(plot1, plot2))
	dual.plot(plot2, file.path(qc.dir,paste(sample,".VariableFeatures.volcano",sep="")), w=7, h=5, res=mid_res)

	counts.dat = round(as(SCrna[[active.assay]]@counts,"matrix"),5)
	write.xls.genename(round(as(SCrna[[active.assay]]@counts,"matrix"),5),
		file.path(qc.dir, paste(sample,".UMI_count.matrix.xls",sep="")),
		GENENAME)
	write.xls.genename(round(as(SCrna[[active.assay]]@data,"matrix"),5), 
		file.path(qc.dir, paste(sample,".normalize_expression.matrix.xls",sep="")), 
		GENENAME)

	pca_feat = data.frame(Gene=HVGs, Genename=map.genename(HVGs,GENENAME,3,1))
	write.table(pca_feat,file=file.path(pca.dir,paste(sample,".PCA.genes.txt",sep="")),quote=F,sep="\t",row.names=F)
	## pca out files
	f = VizDimLoadings(SCrna, dims = 1:2, reduction = "pca")
	dual.plot(f, file.path(pca.dir,paste(sample,".PCA.topGene",sep="")), w=8, h=6, res=mid_res)
}

f = DimHeatmap(SCrna, dims = 1:15, cells = 500, balanced = TRUE, fast=FALSE, raster=FALSE)
dual.plot(f, file.path(pca.dir,paste(sample,".PCA.heatmap",sep="")), w=12, h=20, res=mid_res)

p = ggplot(pc.stdev, aes(PC,stdev)) + geom_point(size=3) + geom_line(aes(PC,fitted(pc.fit)),col='red')
p = p+xlab("PC")+ylab("Standard Deviation of PC")
dual.plot(p, file.path(pca.dir,paste(sample,".PCA.sdev_fitted",sep="")), w=8, h=5, res=mid_res)

f = ElbowPlot(SCrna)
dual.plot(f, file.path(pca.dir,paste(sample,".PCA.ElbowPlot",sep="")), w=8, h=5, res=mid_res)

f = JackStrawPlot(SCrna, dim = 1:15)
dual.plot(f, file.path(pca.dir,paste(sample,".PCA.JackStrawPlot",sep="")), w=8, h=5, res=mid_res)

#########   non-linear dimensional reduction, tSNE / UMAP
tsne_project = data.frame(Barcode = rownames(SCrna@reductions$tsne@cell.embeddings), SCrna@reductions$tsne@cell.embeddings)
colnames(tsne_project) = c("Barcode", "TSNE-1", "TSNE-2")
write.csv(tsne_project, file.path(tsne.dir,paste(sample,".seurat.tSNE.csv",sep="")), quote=F, row.names=F)
use_colors=TRUE
if(use_colors){
	tsneplot = DimPlot(SCrna, reduction = "tsne",label = TRUE, cols=defined_cols[1:10])
}else{
	tsneplot = DimPlot(SCrna, reduction = "tsne",label = TRUE)
}
dual.plot(tsneplot, file.path(tsne.dir,paste(sample,".cluster.tSNE",sep="")), w=8, h=7, res=mid_res)
tsneplot = DimPlot(SCrna, reduction = "tsne", group.by = "orig.ident",cols=sample_cols)
dual.plot(tsneplot, file.path(tsne.dir,paste(sample,".sample.tSNE",sep="")), w=8, h=7, res=mid_res)

umap_project = data.frame(Barcode = rownames(SCrna@reductions$umap@cell.embeddings), SCrna@reductions$umap@cell.embeddings)
colnames(umap_project) = c("Barcode", "UMAP-1", "UMAP-2")
write.csv(umap_project, file.path(umap.dir,paste(sample,".seurat.UMAP.csv",sep="")), quote=F, row.names=F)

if(use_colors){
	umapplot = DimPlot(SCrna, reduction = "umap",label = TRUE, cols=defined_cols[1:10])
}else{
	umapplot = DimPlot(SCrna, reduction = "umap",label = TRUE)
}
dual.plot(umapplot, file.path(umap.dir,paste(sample,".cluster.UMAP",sep="")), w=8, h=7, res=mid_res)
umapplot = DimPlot(SCrna, reduction = "umap",group.by = "orig.ident",cols=sample_cols)
dual.plot(umapplot, file.path(umap.dir,paste(sample,".sample.UMAP",sep="")), w=8, h=7, res=mid_res)

#### cluster
cluster.res = data.frame(Barcode = rownames(SCrna@meta.data), Cluster = SCrna@meta.data[,"seurat_clusters"], Sample=SCrna@meta.data[,"orig.ident"])
#cluster.res = cluster.res[order(cluster.res$Cluster),]
write.csv(cluster.res, file.path(cluster.dir, paste(sample,".seurat.cluster.csv",sep="")), quote=F, row.names=F)
write.table(cluster.res, file.path(cluster.dir, paste(sample,".seurat.cluster.tsv",sep="")), quote=F, row.names=F, sep="\t")
for(each in sample_names){
	cluster.tmp = subset(cluster.res, Sample==each)
	cluster.tmp$Barcode = sapply(as.character(cluster.tmp$Barcode), function(x) strsplit(x,"_")[[1]][1])
	write.csv(cluster.tmp[,1:2], file.path(cluster.dir, paste(sample,".",each,".seurat.cluster.csv",sep="")), quote=F, row.names=F)
	write.table(cluster.tmp[,1:2], file.path(cluster.dir, paste(sample,".",each,".seurat.cluster.tsv",sep="")), quote=F, row.names=F, sep="\t")
}


cluster.matrix = data.frame(Cluster=SCrna@meta.data[,"seurat_clusters"], sample=SCrna@meta.data[,"orig.ident"])
clust_stat = data.frame(table(cluster.matrix))
cluster.bar = ggplot(data=clust_stat,aes(sample,weight=Freq,fill=Cluster))+
    geom_bar(position="fill") +
    ggtitle("Seurat") + coord_flip() +
    theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

if (assay == "Spatial"){
    cluster.bar = cluster.bar + ylab("Relative proportion of spots")
}else{
    cluster.bar = cluster.bar + ylab("Relative proportion of cells")
}
if(use_colors){
	cluster.bar = cluster.bar + scale_fill_manual(values = defined_cols[1:10])
}
dual.plot(cluster.bar, file.path(cluster.dir, paste(sample, ".cluster.barplot", sep="")),w=7,h=max(4,length(sample_names)/2),res=high_res)

cluster.bar2 = ggplot(data=clust_stat,aes(Cluster,weight=Freq,fill=sample))+
    geom_bar(position="fill") +
    ggtitle("Seurat") + coord_flip() +
    theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

if (assay == "Spatial"){
    cluster.bar2 = cluster.bar2 + ylab("Relative proportion of spots")
}else{
    cluster.bar2 = cluster.bar2 + ylab("Relative proportion of cells")
}
if(use_colors){
	cluster.bar2 = cluster.bar2 + scale_fill_brewer(values = defined_cols[1:7])
}
dual.plot(cluster.bar2, file.path(cluster.dir, paste(sample, ".cluster.barplot2", sep="")),w=7,h=max(4,length(sample_names)/2),res=high_res)

if(assay == "Spatial"){
	if(use_colors){
		spatial_plots = list()
		for(i in 1:length(sample_names)){
			spatial_plots[[i]] = plot_SpatialCluster(SCrna, sample.name = sample_names[i], 
				image.name = sample_names[i], legend.text = "Cluster", colours=defined_cols[1:10])
		}
		dual.plot_list(spatial_plots, file.path(cluster.dir, paste(sample,".cluster.Spatial",sep="")), w.each = 6, h.each = 6, res = mid_res)
	}else{
		f = SpatialDimPlot(SCrna, label = TRUE, label.size = 3, combine=FALSE, crop=FALSE, pt.size.factor = 1.2)
		for(x in 1:length(f)){
			f[[x]] = f[[x]] + coord_cartesian()
		}
		p = CombinePlots(plots = f, ncol = ncol.img)
		dual.plot(p, file.path(cluster.dir, paste(sample,".cluster.Spatial",sep="")), w=6*ncol.img, h=6*nrow.img, res=mid_res)
	}
}

########  biomarkers / diff analysis
diff.dir = file.path(outdir, 'diffexp')
dir.create(diff.dir, recursive = TRUE)
diff.exp = FindAllMarkers(
	object = SCrna, 
	only.pos = F, 
    assay = 'Spatial',
	min.pct = 0.25,
	test.use = test_method,
	logfc.threshold = 0,
	return.thresh = 1
)
if(has_genename){
	diff.exp$gene_id = map.genename(diff.exp$gene, GENENAME, 3, 1)
	diff.exp$gene_name = map.genename(diff.exp$gene_id, GENENAME, 1, 2)
}

diff.exp = diff.exp[,c(7,8,9,6,2,1,5,3,4)]
diff.sig.exp = subset(diff.exp, avg_log2FC>logfc.diff & p_val_adj<padjust.diff)
#diff.sig.exp = subset(diff.exp, abs(avg_logFC)>logfc.diff & p_val_adj<padjust.diff)
marker.exp = subset(diff.exp, avg_log2FC>logfc.marker & p_val_adj<padjust.marker & pct.1>0.5 & pct.2<0.5)
marker.avexp.data = AverageExpression(SCrna,features=unique(marker.exp$gene),assay=active.assay,slot="data",verbose=F)[[1]]
#marker.avexp.scale = AverageExpression(SCrna,features=unique(marker.exp$gene),assay=active.assay,slot="scale.data",verbose=F)[[1]]
marker.exp$avg_exp = foreach(gene = as.character(marker.exp$gene), clt = as.character(marker.exp$cluster), .combine="c") %do% {
	marker.avexp.data[gene,clt]
}

top4.markers <- diff.exp %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
top10.markers <- diff.exp %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

diff.file = file.path(diff.dir,paste(sample,".cluster.diffexp.xls",sep=""))
diff.sig.file = file.path(diff.dir,paste(sample,".cluster.diffexp.significant.xls",sep=""))
marker.file = file.path(diff.dir,paste(sample,".cluster.markers.xls",sep=""))
marker.top10.file = file.path(diff.dir,paste(sample,".cluster.top10.markers.xls",sep=""))
write.table(diff.exp[,-1], file=diff.file, quote=F, sep="\t", row.names=F)
write.table(diff.sig.exp[,-1], file=diff.sig.file, quote=F, sep="\t", row.names=F)
write.table(marker.exp[,-1], file=marker.file, quote=F, sep="\t", row.names=F)

if(nrow(marker.exp)>0){
	select_markers = unique(marker.exp$gene)[1:100]
	p = DotPlot(SCrna, features = select_markers) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	dual.plot(p, file.path(diff.dir,paste(sample,".cluster.markers.dotplot",sep="")), 
		w=length(select_markers)/5+4, h=length(cluster.ids)/2)
	write.table(p$data, file = file.path(diff.dir,paste(sample,".markers.ggplot_data.txt",sep="")),sep="\t",row.names=F,quote=F)
}



#f = DoHeatmap(SCrna, features = top10.markers$gene) #+ NoLegend()
#dual.plot(f, file.path(diff.dir,paste(sample,".cluster.top10.genename.heatmap",sep="")), w=10, h=8, res=high_res)
if(use_colors){
	f = DoHeatmap(SCrna, features = unique(top10.markers$gene), raster=FALSE, group.colors=defined_cols)
}else{
	f = DoHeatmap(SCrna, features = unique(top10.markers$gene), raster=FALSE)
}
dual.plot(f, file.path(diff.dir,paste(sample,".cluster.top10.heatmap",sep="")), w=10, h=6, res=high_res)

AverageExp = AverageExpression(SCrna, features=unique(top10.markers$gene))
if(assay == "RNA"){
	corelation = corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
}else{
	corelation = corr.test(AverageExp$Spatial,AverageExp$Spatial,method="spearman")
}
f = pheatmap(corelation$r)
dual.plot(f, file.path(diff.dir,paste(sample,".cluster.correlation",sep="")), res=high_res)

jam.pdfs = NULL
for(i in cluster.ids){
	## violin plot
	genes.plot = top4.markers$gene[top4.markers$cluster==i]
	print(genes.plot)
	v = VlnPlot(SCrna, genes.plot, ncol =2, pt.size = 0) + xlab("Cluster") + ylab("log(UMI)")
	dual.plot(v, file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.violin",sep="")), w=10, h=7, res=mid_res)
	## feature plot
	f1 = FeaturePlot(object = SCrna, genes.plot, cols = c("grey", "blue"), reduction = "tsne")
	dual.plot(f1, file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.tsne",sep="")), res=mid_res)
	f2 = FeaturePlot(object = SCrna, genes.plot, cols = c("grey", "blue"), reduction = "umap")
	dual.plot(f2, file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.umap",sep="")), res=mid_res)
	jam.pdfs = c(jam.pdfs, file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.tsne.pdf",sep="")),
		file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.umap.pdf",sep="")),
		file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.violin.pdf",sep="")))
	## spatial
	if(assay == "Spatial"){
		f = SpatialFeaturePlot(SCrna, features = c(genes.plot), pt.size.factor = 1, combine=FALSE, crop=crop)
		for(x in 1:length(f)){
			f[[x]] = f[[x]] + coord_cartesian()
		}
		p <- CombinePlots(plots = f, ncol = length(sample_names))
		dual.plot(p, file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.spatial",sep="")), w=length(sample_names)*3, h=4*4, res=mid_res)
		jam.pdfs = c(jam.pdfs, file.path(diff.marker.dir,paste(sample,".cluster",i,".marker.spatial.pdf",sep="")))
	}
}

violin.pngs = file.path("markers",paste(sample,".cluster",cluster.ids,".marker.violin.png",sep=""))
names(violin.pngs) = paste("cluster",cluster.ids,"_violin",sep="")
write.table(violin.pngs, file=file.path(diff.src.dir,"violin.list"), quote=F, sep="\t", col.names=F)

umap.pngs = file.path("markers",paste(sample,".cluster",cluster.ids,".marker.umap.png",sep=""))
names(umap.pngs) = paste("cluster",cluster.ids,"_umap",sep="")
write.table(umap.pngs, file=file.path(diff.src.dir,"umap.list"), quote=F, sep="\t", col.names=F)

tsne.pngs = file.path("markers",paste(sample,".cluster",cluster.ids,".marker.tsne.png",sep=""))
names(tsne.pngs) = paste("cluster",cluster.ids,"_tsne",sep="")
write.table(tsne.pngs, file=file.path(diff.src.dir,"tsne.list"), quote=F, sep="\t", col.names=F)
if(assay == "Spatial"){
	spatial.pngs = file.path("markers",paste(sample,".cluster",cluster.ids,".marker.spatial.png",sep=""))
	names(spatial.pngs) = paste("cluster",cluster.ids,"_spatial",sep="")
	write.table(spatial.pngs, file=file.path(diff.src.dir,"spatial.list"), quote=F, sep="\t", col.names=F)
}

if(TRUE){
	jamPDF(jam.pdfs,
		out.file = file.path(diff.dir,paste(sample,".cluster.markers",sep="")),
		layout = '2x2',
		delete.original = FALSE,
		ignore.stderr = FALSE)
}


