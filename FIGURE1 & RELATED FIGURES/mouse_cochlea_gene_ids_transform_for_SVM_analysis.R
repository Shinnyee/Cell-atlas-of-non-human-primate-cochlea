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
adata_loom <- connect(filename = "mouse_coch.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('mouse_coch_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('mouse_coch_var.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene

#subsequential scripts are same as in the code from "macaque_cochlea_gene_ids_transform_for_SVM_analysis"









mouse= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                            
                            min.cells = 0, 
                            min.features = 0)

mouse@assays[["RNA"]]@meta.features <- meta_feature
table(mouse$batch, mouse$cell_type)
mouse_counts <- mouse@assays$RNA@counts
mouse_counts <- data.frame(gene=rownames(mouse_counts), mouse_counts, check.names = F)
#gene_id transform
# 转小鼠基因名为人基因名
genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
mouse_counts$Gene <- genesV2[match(mouse_counts$gene, genesV2[,1]),2]
mouse_counts <- subset(mouse_counts, Gene!='NA')
mouse_counts <- dplyr::select(mouse_counts, Gene, everything())
mouse_counts <- mouse_counts[, !(colnames(mouse_counts) %in% 'gene')]
dim(mouse_counts)
write.csv(mouse_counts,file="rawcounts_mouse_to_human_transformed.csv")
meta_mouse=mouse@meta.data

rownames(mouse_counts)<-mouse_counts[,1]
mouse_counts<-mouse_counts[,-1]
head(rownames(mouse_counts))




seurat_mouse <- CreateSeuratObject(counts = mouse_counts)

seurat_mouse












