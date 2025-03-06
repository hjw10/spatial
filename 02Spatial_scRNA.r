#! usr/bin/R
### zhaoyunfei
### 20200111

suppressMessages({
library(Seurat)
library(argparse)
library(dplyr)
library(ggplot2)
library(SPOTlight)
})

parser = ArgumentParser()
parser$add_argument("--sc_rds", help="the sc data",required = T)
parser$add_argument("--spatial_rds", help="the sp data",required = T)
parser$add_argument("--sample", help="the sample name",required = T)
parser$add_argument("--outdir", help="the outdir",default = './')
parser$add_argument("--celltype", help="the annotation for celltype")
parser$add_argument("--normalization_method",default = 'SCT')
parser$add_argument("--reduction",help='Dimensional reduction to perform when finding anchors',choices = c('pcaproject','cca'),default = 'pcaproject')
parser$add_argument("--img",required = T)
parser$add_argument("--cell_interest", help="the interest of cell type,eg : 'T,B:Fobro,epi'")
args <- parser$parse_args()
str(args)


normalization_method = "SCT"
reduction = "cca"
img = c("/public/workspace/stu21230110/SPA/C1/spatial/tissue_hires_image.png",
        "/public/workspace/stu21230110/SPA/C2/spatial/tissue_hires_image.png",
        "/public/workspace/stu21230110/SPA/C3/spatial/tissue_hires_image.png",
        "/public/workspace/stu21230110/SPA/P1/spatial/tissue_hires_image.png",
        "/public/workspace/stu21230110/SPA/P2/spatial/tissue_hires_image.png",
        "/public/workspace/stu21230110/SPA/P3/spatial/tissue_hires_image.png")

setwd("/public/workspace/stu21230110/SPA_result")
outdir = "02Spatial_scRNA"
if (!file.exists(outdir)) {dir.create(outdir,recursive = TRUE)}
load("~/scRNA/06_scRNA_harmony_singler202400810.Rdata")
head(scRNA_harmony_singler$seurat_clusters)

##如果不是SCT而是RNA需要强行矫正
scRNA <- scRNA_harmony_singler
scRNA <- SCTransform(scRNA, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:20)

scRNA <- Seurat::FindNeighbors(scRNA, dims = 1:20, verbose = FALSE)
scRNA <- Seurat::FindClusters(scRNA, verbose = FALSE)
#scRNA <- Seurat::FindClusters(scRNA, verbose = FALSE, resolution = 0.5)
table(scRNA$seurat_clusters)
saveRDS(scRNA,"02Spatial_scRNA/scRNA.rds")

if ( ! is.null(celltype) ) {
anno = read.csv(celltype,header = T)
colnames(anno) = c('Barcode','Cluster')
if (length(anno$Cluster) != length(colnames(scRNA_harmony_singler))){print("warning ~~~ ,The length of celltype file is not match the sc data,subset will be acrry out ~~~")}
index = na.omit(match(colnames(scRNA_harmony_singler),anno$Barcode))
scRNA_harmony_singler = scRNA_harmony_singler[,index]
Seurat::Idents(object = scRNA_harmony_singler) <- anno$Cluster
scRNA_harmony_singler$seurat_clusters = anno$Cluster
if (! is.null(anno$Sample)){scRNA_harmony_singler$orig.ident = anno$Sample}
}else {Seurat::Idents(object = scRNA_harmony_singler) <- scRNA_harmony_singler$seurat_clusters
}


cortex_sp = readRDS("/public/workspace/stu21230110/SPA_result/01CCA_Seurat/combined.seurat.rds")
anchors <- FindTransferAnchors(reference = scRNA_harmony_singler, query = cortex_sp, normalization.method = normalization_method,reduction = reduction)
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_harmony_singler$seurat_clusters, prediction.assay = TRUE, 
    weight.reduction = cortex_sp[["pca"]],dims = 1:20)

cortex_sp[["predictions"]] <- predictions.assay
DefaultAssay(cortex_sp) <- "predictions"

cells = rownames(cortex_sp@assays$predictions@data)[which(rownames(cortex_sp@assays$predictions@data) != 'max')]
f = SpatialFeaturePlot(cortex_sp, features = cells,alpha = c(0.1,1),combine=TRUE,pt.size.factor = 1.3,stroke = 0,crop = FALSE,ncol = 3,min.cutoff = 0, max.cutoff = 1)
for(x in 1:length(cells)){
                f[[x]] = f[[x]] + coord_cartesian() + theme(legend.position = 'right') + labs(title = cells[x],fill = 'value') + theme(plot.title = element_text(hjust = 0.5,size = 12)) + scale_fill_gradient(low = 'white',high = '#ed1941') + theme(plot.title = element_text(size = 30))
        }


plot_list = list()

for (i in 1:length(cells)){plot_list[[i]] = f[[i]]}


decon_mtrx = t(cortex_sp@assays$predictions@data)

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "max")]

decon_df <- decon_mtrx %>%
  data.frame(check.names = F) %>%
  tibble::rownames_to_column("barcodes")

#decon_df$barcodes = rownames(tmp)

cortex_sp@meta.data <- cortex_sp@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

###plot dot
slice <- names(cortex_sp@images)[1]
metadata_ds <- data.frame(cortex_sp@meta.data)
colnames(metadata_ds) <- colnames(cortex_sp@meta.data)
cell_types_interest <- cell_types_all

metadata_ds <- metadata_ds %>% tibble::rownames_to_column("barcodeID") %>%
            dplyr::mutate(rsum = base::rowSums(.[, cell_types_interest,
                drop = FALSE])) %>% dplyr::filter(rsum != 0) %>%
            dplyr::select("barcodeID") %>% dplyr::left_join(metadata_ds %>%
            tibble::rownames_to_column("barcodeID"), by = "barcodeID") %>%
            tibble::column_to_rownames("barcodeID")


spatial_coord <- data.frame(cortex_sp@images[[slice]]@coordinates) %>%
        tibble::rownames_to_column("barcodeID") %>% dplyr::mutate(imagerow_scaled = imagerow *
        cortex_sp@images[[slice]]@scale.factors$lowres, imagecol_scaled = imagecol *
        cortex_sp@images[[slice]]@scale.factors$lowres) %>% dplyr::inner_join(metadata_ds %>%
        tibble::rownames_to_column("barcodeID"), by = "barcodeID")

i=1
img <- png::readPNG(img[i])
img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, 
        "npc"), height = grid::unit(1, "npc"))
for (cell in cell_types_all){
Max = max(spatial_coord[,cell])
spatial_coord[,cell] = spatial_coord[,cell]/Max
scatterpie_plt <- suppressMessages(ggplot2::ggplot() + ggplot2::annotation_custom(grob = img_grob, 
        xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) + 
        ggplot2::geom_point(data = spatial_coord, ggplot2::aes(x = imagecol_scaled, 
            y = imagerow_scaled,size = get(cell),alpha = get(cell)), color = '#FF4500') +
        ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img), 
        0) + ggplot2::xlim(0, ncol(img)) + cowplot::theme_half_open(11, 
        rel_small = 1) + ggplot2::theme_void() + ggplot2::coord_fixed(ratio = 1, 
        xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +ggplot2::scale_size_continuous(range=c(0,2))+ggplot2::scale_alpha_continuous(range=c(0,1))+labs(size = cell) + guides(alpha = "none"))
pdf(paste(outdir,paste(cell,'dot.pdf',sep = '.'),sep = '/'),width = 8,height = 7)
print(scatterpie_plt)
dev.off()
}


scatterpie_pie <- suppressMessages(ggplot2::ggplot() + ggplot2::annotation_custom(grob = img_grob,
        xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
        scatterpie::geom_scatterpie(data = spatial_coord, ggplot2::aes(x = imagecol_scaled,
            y = imagerow_scaled), cols = cell_types_all, color = NA,
            alpha = 1, pie_scale = 0.3) +
        ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img),
        0) + ggplot2::xlim(0, ncol(img)) + cowplot::theme_half_open(11,
        rel_small = 1) + ggplot2::theme_void() + ggplot2::coord_fixed(ratio = 1,
        xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))

pdf(paste(outdir,paste(sample,'pie.pdf',sep = '.'),sep = '/'),width = 8,height = 7)
print(scatterpie_pie)
dev.off()




p = CombinePlots(plots = plot_list, ncol = 4 ,legend = 'bottom')


pdf(paste(outdir,sprintf('%s.sp.predictions.combined.pdf',sample),sep ='/'),width = 16,height = 4*(ceiling((length(unique(Idents(scRNA_harmony_singler)))/3))))

print(p)

dev.off()


pdf(paste(outdir,sprintf('%s.sp.predictions.pdf',sample),sep ='/'),width = 16,height = 4*(ceiling((length(unique(Idents(scRNA_harmony_singler)))/3))))
print(f)
dev.off()

predictions.data = t(as.data.frame(predictions.assay@data))

write.csv(predictions.data,file = paste(outdir,sprintf('%s.sp.predictions.csv',sample),sep = '/'),quote = F)


if (!is.null(cell_interest)){
cells = strsplit(cell_interest,':')

for (len in 1:length(cell[[1]])){
cell = strsplit(cell[[1]][i],',')[[1]]
pdf(paste(outdir,paste(sample,paste(cell[[1]], collapse = '.'),'deconvolution.pdf',sep = '.'),sep = '/'))
scatterpie_pie <- suppressMessages(ggplot2::ggplot() + ggplot2::annotation_custom(grob = img_grob,
        xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
        scatterpie::geom_scatterpie(data = spatial_coord, ggplot2::aes(x = imagecol_scaled,
            y = imagerow_scaled), cols = cell, color = NA,
            alpha = 1, pie_scale = 0.3) +
        ggplot2::scale_y_reverse() + ggplot2::ylim(nrow(img),
        0) + ggplot2::xlim(0, ncol(img)) + cowplot::theme_half_open(11,
        rel_small = 1) + ggplot2::theme_void() + ggplot2::coord_fixed(ratio = 1,
        xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))

print(scatterpie_pie)
dev.off()

}
}



saveRDS(cortex_sp,file=paste(outdir,sprintf('%s.sp.predictions.rds',sample),sep = '/'))










