library(future)
plan("multiprocess",workers=20)
###future.globals.maxSize= X,x的单位是字节，下面这句代码是8个G
options(future.globals.maxSize= 80000*1024^2) 
library(data.table)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
#install.packages("devtools")
#library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
#BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(DoubletFinder)

set.seed(123456)
setwd("./03cell2loc")




#去除双细胞
data_directory=c("~/design/data/scRNA/C1/filtered_feature_bc_matrix/",
                 "~/design/data/scRNA/C2/filtered_feature_bc_matrix/",
                 "~/design/data/scRNA/C3/filtered_feature_bc_matrix/",
                 "~/design/data/scRNA/T1/filtered_feature_bc_matrix/",
                 "~/design/data/scRNA/T2/filtered_feature_bc_matrix/",
                 "~/design/data/scRNA/T3/filtered_feature_bc_matrix/")

project_name<-c("C1","C2","C3","P1","P2","P3")
samples <- project_name

colon.data <- Read10X(data.dir = data_directory[1])
currentSample <- CreateSeuratObject(counts = colon.data, project = project_name[1], min.cells = 3, min.features = 200)

seu_list <- currentSample
for (i in 2:length(samples)){
  colon.data <- Read10X(data.dir = data_directory[i])
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name[i], min.cells = 3, min.features = 200) 
  seu_list=merge(seu_list,currentSample)
  
}
scRNA=seu_list

load("~/scRNA/06_")
celltype <- cbind(rownames(scRNA_harmony_singler@meta.data),scRNA_harmony_singler@meta.data$singleRnew)
colnames(celltype) <- c("Barcode","celltype")
write.csv(celltype,"celltype.csv")


library(SeuratDisk)    
#seurat2h5seurat中间过渡	
SaveH5Seurat(scRNA,filename="scRNA.h5seurat", overwrite = TRUE)
#数据转为最终h5ad格式
Convert("scRNA.h5seurat", dest = "h5ad", overwrite = TRUE)

