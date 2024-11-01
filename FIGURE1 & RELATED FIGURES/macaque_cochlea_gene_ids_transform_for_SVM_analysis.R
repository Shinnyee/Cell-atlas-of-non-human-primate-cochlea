rm(list = ls())
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
# ============ Load data ============ 

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(zellkonverter)
#BiocManager::install("zellkonverter")
library(loomR)
library(Seurat)
library(SeuratDisk)
library(data.table)
adata_loom <- connect(filename = "macaque_coch.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('macaque_coch_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('macaque_coch_var.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene


macaque= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                          
                          min.cells = 0, 
                          min.features = 0)

macaque@assays[["RNA"]]@meta.features <- meta_feature
table(macaque$age, macaque$cell_type)
macaque_counts <- macaque@assays$RNA@counts
macaque_counts <- data.frame(gene=rownames(macaque_counts), macaque_counts, check.names = F)
#gene_id transform
# 转小鼠基因名为人基因名
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
genesV2 <- read.csv("orthologTable_human_macaque_mouse.csv", row.names = 1)
macaque_counts$Gene <- genesV2[match(macaque_counts$gene, genesV2[,5]),4]
macaque_counts <- subset(macaque_counts, Gene!='NA')
macaque_counts <- dplyr::select(macaque_counts, Gene, everything())
macaque_counts <- macaque_counts[, !(colnames(macaque_counts) %in% 'gene')]
dim(macaque_counts)
write.csv(macaque_counts,file="rawcounts_macaque_to_human_transformed.csv")
meta_macaque=macaque@meta.data

rownames(macaque_counts)<-macaque_counts[,1]
macaque_counts<-macaque_counts[,-1]
head(rownames(macaque_counts))

seurat_macaque <- CreateSeuratObject(counts = macaque_counts)

seurat_macaque
seurat_macaque <- AddMetaData(seurat_macaque, metadata = meta_macaque)
x_scvi = adata_loom$col.attrs$X_scVI[,]
x_umap = adata_loom$col.attrs$X_umap[,]
x_pca=adata_loom$col.attrs$X_pca[,]
x_scvi = t(x_scvi)
x_umap = t(x_umap)
x_pca = t(x_pca)
rownames(x_scvi) = barcode
rownames(x_umap) = barcode
rownames(x_pca) = barcode
colnames(x_scvi) = c("scVI_1","scVI_2","scVI_3","scVI_4","scVI_5",
                     "scVI_6","scVI_7","scVI_8","scVI_9","scVI_10")
colnames(x_umap) = c('UMAP_1','UMAP_2')


colnames(x_pca) = c('PCA_1','PCA_2','PCA_3','PCA_4','PCA_5','PCA_6','PCA_7','PCA_8',
                    'PCA_9','PCA_10','PCA_11','PCA_12','PCA_13','PCA_14','PCA_15','PCA_16',
                    'PCA_17','PCA_18','PCA_19','PCA_20','PCA_21','PCA_22','PCA_23','PCA_24',
                    'PCA_25','PCA_26','PCA_27','PCA_28','PCA_29','PCA_30','PCA_31','PCA_32',
                    'PCA_33','PCA_34','PCA_35','PCA_36','PCA_37','PCA_38','PCA_39','PCA_40',
                    'PCA_41','PCA_42','PCA_43','PCA_44','PCA_45','PCA_46','PCA_47','PCA_48',
                    'PCA_49','PCA_50')


sce1 = as.SingleCellExperiment(seurat_macaque)
sce1@int_colData@listData[["reducedDims"]]@listData[["scVI"]] = x_scvi
sce1@int_colData@listData[["reducedDims"]]@listData[["umap"]] = x_umap
sce1@int_colData@listData[["reducedDims"]]@listData[["PCA"]] = x_pca

seurat_macaque <- as.Seurat(sce1)
DimPlot(seurat_macaque,reduction = "umap",group.by = "cell_type")
DimPlot(seurat_macaque,reduction = "umap",group.by = "age")
# h5ad format transform
library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
sceasy::convertFormat(seurat_macaque, from="seurat", to="anndata",
                      outFile='Mac_coch_python.h5ad')



