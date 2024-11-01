rm(list = ls())
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
# In this script, we are going to anchor MACAQUE COE snRNA-seq data to mouse COE snRNA-seq data and transfer the mouse cell type labels to each of the MACAQUE cells
# Adapted from Seurat tutorial: https://satijalab.org/seurat/articles/integration_mapping.html

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
adata_loom <- connect(filename = "macaque_all.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('macaque_all_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('macaque_all_var.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene
macaque= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                
                                min.cells = 0, 
                                min.features = 0)

macaque@assays[["RNA"]]@meta.features <- meta_feature
table(macaque$age, macaque$cell_type)


DefaultAssay(macaque) <- "RNA"
macaque$species <- "macaque"
Idents(macaque) <- "age"
table(macaque$age)
macaque <- SCTransform(macaque, verbose = T, 
                           vars.to.regress = c("nCount_RNA"), 
                           conserve.memory = T)
macaque <- RunPCA(macaque,npcs = 50)
DimPlot(macaque,reduction = "pca")
ElbowPlot(macaque , ndims = 50)
macaque <- FindNeighbors(macaque, reduction = "pca", dims = 1:30)
macaque <- FindClusters(macaque)
macaque <- RunUMAP(macaque, reduction = "pca", dims = 1:20)
# Visualization
DimPlot(macaque, reduction = "umap",label = TRUE,group.by = 'age')
#subtract cochlear epithelium clusters
Idents(macaque) <- "cell_type"
table(Idents(macaque))
macaque_coe=subset(macaque,idents=c("HC","CC_ISC_OSC","IDC",
                                    "IBC_IPh_HeC","DC_PC","TBC"))
macaque_coe <- RunUMAP(macaque_coe, reduction = "pca", dims = 1:10)
DimPlot(macaque_coe, reduction = "umap",label = TRUE,
        group.by = 'cell_type')
library(future)
options(future.globals.maxSize = 10000 * 1024^2)
# JOINT CCA-EMBEDING 
DefaultAssay(macaque_coe) <- "RNA"
# reintegration of human cochlea epithelium by using cca-based integration
combined.list <- SplitObject(macaque_coe, split.by = "age")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], 
                                    vars.to.regress = c("nCount_RNA"),
                                    verbose = TRUE, method="glmGamPoi")
}
#combined.list<-lapply(X=combined.list,FUN = SCTransform, method="glmGamPoi")
features <- SelectIntegrationFeatures(object.list =combined.list, 
                                      nfeatures = 2000)
write.csv(features,file = "anchors_for_cca_integration_macaque_coe_v3.csv")


combined.list <- PrepSCTIntegration(object.list =combined.list, 
                                    anchor.features = features)
#combined.list <- lapply(X = combined.list, FUN = RunPCA, features = features)
cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, 
                                          normalization.method = "SCT",
                                          anchor.features = features,
                                          dims = 1:40)

cochlea.combined.sct <- IntegrateData(anchorset = cochlea.anchors,
                            normalization.method = "SCT",dims = 1:40)

cochlea.combined.sct <- RunPCA(cochlea.combined.sct,npcs=200)
ElbowPlot(cochlea.combined.sct , ndims = 200)

DefaultAssay(cochlea.combined.sct) <- "integrated"
cochlea.combined.sct <- FindNeighbors(cochlea.combined.sct, 
                                      reduction = "pca", 
                                      dims = 1:50)
cochlea.combined.sct <- FindClusters(cochlea.combined.sct)
cochlea.combined.sct <- RunUMAP(cochlea.combined.sct, reduction = "pca", 
                                dims = 1:15)
p1 <- DimPlot(cochlea.combined.sct, reduction = "umap", 
              group.by = "age")
p2 <- DimPlot(cochlea.combined.sct, reduction = "umap", 
              group.by = "cell_type",
              label = TRUE)
p1+p2

macaque_coe=cochlea.combined.sct
saveRDS(macaque_coe,file = "maaque_coe_cca_integration.rds")
#######################################=====================================
##########============================================================
adata_loom <- connect(filename = "mouse_coe_all.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('mouse_coe_all_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('mouse_coe_all_var.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene
mouse= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                            
                            min.cells = 0, 
                            min.features = 0)

mouse@assays[["RNA"]]@meta.features <- meta_feature
table( mouse$cell_type)


DefaultAssay(mouse) <- "RNA"
mouse$species <- "mouse"


mouse <- SCTransform(mouse, verbose = T, 
                       vars.to.regress = c("nCount_RNA"), 
                       conserve.memory = T)
mouse <- RunPCA(mouse,npcs = 50)
DimPlot(mouse,reduction = "pca")
ElbowPlot(mouse , ndims = 50)
mouse <- FindNeighbors(mouse, reduction = "pca", dims = 1:30)
mouse <- FindClusters(mouse)
mouse <- RunUMAP(mouse, reduction = "pca", dims = 1:10)
# Visualization
DimPlot(mouse, reduction = "umap",label = TRUE,group.by = 'cell_type')
mouse_coe=mouse

# JOINT CCA-EMBEDING 
DefaultAssay(mouse_coe) <- "RNA"
mouse_coe$batch
# reintegration of human cochlea epithelium by using cca-based integration
combined.list <- SplitObject(mouse_coe, split.by = "batch")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], 
                                    vars.to.regress = c("nCount_RNA"),
                                    verbose = TRUE, method="glmGamPoi")
}
#combined.list<-lapply(X=combined.list,FUN = SCTransform, method="glmGamPoi")
features <- SelectIntegrationFeatures(object.list =combined.list, 
                                      nfeatures = 3000)
write.csv(features,file = "anchors_for_cca_integration_mouse_coe_v3.csv")


combined.list <- PrepSCTIntegration(object.list =combined.list, 
                                    anchor.features = features)
#combined.list <- lapply(X = combined.list, FUN = RunPCA, features = features)
cochlea.anchors <- FindIntegrationAnchors(object.list = combined.list, 
                                          normalization.method = "SCT",
                                          anchor.features = features,
                                          dims = 1:30)

cochlea.combined.sct <- IntegrateData(anchorset = cochlea.anchors,
                                      normalization.method = "SCT",
                                      dims = 1:20)

cochlea.combined.sct <- RunPCA(cochlea.combined.sct,npcs=100)
ElbowPlot(cochlea.combined.sct , ndims = 100)

DefaultAssay(cochlea.combined.sct) <- "integrated"
cochlea.combined.sct <- FindNeighbors(cochlea.combined.sct, 
                                      reduction = "pca", 
                                      dims = 1:50)
cochlea.combined.sct <- FindClusters(cochlea.combined.sct)
cochlea.combined.sct <- RunUMAP(cochlea.combined.sct, reduction = "pca", 
                                dims = 1:8)
p1 <- DimPlot(cochlea.combined.sct, reduction = "umap", 
              group.by = "batch")
p2 <- DimPlot(cochlea.combined.sct, reduction = "umap", 
              group.by = "cell_type",
              label = TRUE)
p1+p2
dev.off()
mouse_coe=cochlea.combined.sct
saveRDS(mouse_coe,file = "mouse_coe_cca_integration.rds")
#mouse_coe=subset(mouse_coe,idents="TBC",invert=T)
#macaque_coe=subset(macaque_coe,idents="TBC",invert=T)
################################################################################
################################################################################
################################################################################
################################################################################
#KEEP CELL TYPE CONSISTENCE
Idents(mouse_coe) <- "cell_type"
table(Idents(mouse_coe))
Idents(macaque_coe) <- "cell_type"
table(Idents(macaque_coe))
seurat_macaque_coe=macaque_coe
seurat_mouse_coe=mouse_coe
p1 <- DimPlot(seurat_macaque_coe, reduction = "umap", 
              group.by = "cell_type")
p2 <- DimPlot(seurat_mouse_coe, reduction = "umap", 
              group.by = "cell_type",
              label = TRUE)
p1+p2
p1 <- DimPlot(seurat_macaque_coe, reduction = "umap", 
              group.by = "age")
p2 <- DimPlot(seurat_mouse_coe, reduction = "umap", 
              group.by = "batch",
              label = TRUE)
p1+p2
# Names of macaque genes 
#matrix_macaque=seurat_macaque_coe@assays$RNA@counts


#meta_macaque=seurat_macaque_coe@meta.data

# create a seurat object from macaque dataset that is compatible with mouse dataset
#seurat_macaque_coe <- CreateSeuratObject(counts = matrix_macaque)
seurat_macaque_coe
seurat_mouse_coe
#preprocessing of macaque data
#seurat_macaque_coe <- AddMetaData(seurat_macaque_coe, metadata = meta_macaque)
seurat_macaque_coe <- SCTransform(object = seurat_macaque_coe,
                                  vars.to.regress = c("nCount_RNA"),
                                  conserve.memory = T)
#seurat_macaque_coe <- RunPCA(seurat_macaque_coe)
#seurat_mouse_coe <- PrepSCTFindMarkers(seurat_mouse_coe)
seurat_mouse_coe <- SCTransform(object = seurat_mouse_coe,vars.to.regress = c("nCount_RNA"),conserve.memory = T)
#seurat_mouse_coe <- RunPCA(seurat_mouse_coe)
# ============ projection of reference data onto query object ============ 
# We use all default parameters here for identifying anchors
transfer.anchors <- FindTransferAnchors(reference = seurat_mouse_coe,
                                        query = seurat_macaque_coe, 
                                        normalization.method = "SCT",
                                        reference.assay = "SCT", 
                                        query.assay = "SCT", 
                                        reduction = "cca",dims = 1:30)
#macaque_coe <- RunPCA(macaque_coe,npcs = 30)
#mouse_coe <- RunPCA(mouse_coe,npcs = 30)
# TransferData returns a matrix with predicted IDs and prediction scores, which are written into meta.data
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                        refdata = seurat_mouse_coe$cell_type, 
                        weight.reduction =seurat_macaque_coe[["pca"]] ,
                                     dims =1:20)#seurat_macaque_coe[["pca"]]
seurat_macaque_coe <- AddMetaData(seurat_macaque_coe, 
                                  metadata = celltype.predictions)


# Cells with low anchoring scores (<=0.25) are either low-quality cells or cells with ambigous identity.
# They are marked as low-confidence cells and removed later.
cutoff=0.25
seurat_macaque_coe$hiConfi=ifelse(seurat_macaque_coe$prediction.score.max > cutoff,'TRUE','FALSE')

# summary report of anchoring results
table(seurat_macaque_coe$predicted.id == seurat_macaque_coe$cell_type)
seurat_macaque_coe@meta.data%>%mutate(agree.cell_type=
                                        (predicted.id==cell_type))%>%
  group_by(cell_type)%>%summarise(n=sum(agree.cell_type)/n())


# ============ Intergration of two datasets for visualization ============ 

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_mouse_coe)
refdata <- GetAssayData(seurat_mouse_coe, assay = "SCT", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, 
                           refdata = refdata, 
                           weight.reduction = seurat_macaque_coe[["pca"]],
                           dims = 1:10)
# this line adds the imputed data matrix to the seurat_MACAQUE object
seurat_macaque_coe[["SCT"]] <- imputation
seurat_macaque_coe$species = "Macaque"
seurat_mouse_coe$species= "Mouse"
coembed <- merge(x = seurat_mouse_coe, y = seurat_macaque_coe)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)

coembed <- RunUMAP(coembed, dims = 1:30)
# umap cluster colors
Idents(coembed) <- "cell_type"
table(Idents(coembed))
levels(Idents(coembed))

my_cols <- c('TBC'='#DDA0DD','DC_PC'='#D2691E','IBC_IPh_HeC'='#53ad7e',
             'CC_ISC_OSC'='#872657',
             'HC'='#2891a8','IDC'='#9932cc'
             
)
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
DimPlot(coembed,
        cols = my_cols2, label=TRUE , repel=TRUE,reduction = "umap",pt.size = 1)
my_cols3 <- c( "royalblue1","maroon4")

p1 <- DimPlot(coembed,reduction = "umap",label = TRUE,repel = TRUE,
              group.by = "species",cols = my_cols3,pt.size = 1.5)
p2 <- DimPlot(coembed,reduction = "umap",label = TRUE,repel = TRUE,
              group.by = "cell_type",cols = my_cols2,pt.size = 1.5)
p1+p2

cells <- rownames(seurat_macaque_coe@meta.data[seurat_macaque_coe@meta.data$prediction.score.max>0.25,])

DimPlot(seurat_macaque_coe, cells= cells, group.by="predicted.id", 
        label=T,reduction = "umap",pt.size = 1)

p1 <- DimPlot(seurat_macaque_coe,reduction = "umap",label = TRUE,repel = TRUE,cols = my_cols2,
              group.by = "cell_type",pt.size = 1.2)
p2 <- DimPlot(seurat_macaque_coe, cells= cells, group.by="predicted.id", 
              label=T,reduction = "umap",pt.size = 1.2)
p3 <- DimPlot(seurat_macaque_coe, cells= cells, group.by="age", 
              label=T,reduction = "umap",pt.size = 1.2)
library(patchwork)
patchwork::wrap_plots(p1+p2+p3)
p1+p2+p3
dev.off()

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets

macaque_coe
mouse_coe
VlnPlot(mouse_coe,features = "nFeature_RNA")
coembed = AddMetaData(coembed,Embeddings(coembed[["umap"]]),
                      colnames(Embeddings(coembed[["umap"]])))

coembed$species=ifelse(!is.na(coembed$predicted.id), 'Macaque', 'Mouse')

# summary plots
d1 <- DimPlot(coembed, group.by = "species")
d2 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Mouse')%>%pull(V1),group.by = "cell_type", label = TRUE, repel = TRUE)+labs(title='Mouse subclasses')

d3 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Macaque')%>%pull(V1),group.by = "cell_type", label = TRUE, repel = TRUE)+labs(title='Macaque subclasses')
d4 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Macaque')%>%pull(V1),group.by = "predicted.id", label = TRUE, repel = TRUE)+labs(title='Macaque predicted.id')
library(viridis)
d5 = FeaturePlot(coembed, 'prediction.score.max',reduction = "umap",cols =viridis(100))
d6 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'))


grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)
DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%pull(V1),group.by = "cell_type", label = TRUE, repel = TRUE)+labs(title='subclasses')

saveRDS(coembed,file='Seurat_coe_LABEL_TRNSFER_v2.rds')

table(seurat_macaque_coe$cell_type)
table(seurat_macaque_coe$cell_type,seurat_macaque_coe$predicted.id)

dev.off()

results <- table(seurat_macaque_coe@meta.data$cell_type, seurat_macaque_coe@meta.data$predicted.id)

# IF NO macaque CELL MAPS TO A CERTAIN MOUSE IDENTITY, THAT MOUSE IDENTITY WILL BE MISSING IN THE FINAL HEATMAP!!!
library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(loomR)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)

results.norm <- 100*(results/rowSums(results))
mus.ids <- colnames(results.norm)
results.norm <- results.norm[,mus.ids]
#saveRDS(results.norm, file=paste0("macaque_query_mouse_reference_normalized_ClusterID_vargenes_",nVarGenes,"_dims_",nDims_integration_anchor,".rds"))
#label_transfer_matrix_mus_reference <- readRDS(file="macaque_query_mouse_reference_normalized_ClusterID_vargenes_3000_dims_30.rds")

#pdf(file=paste0("macaque_query_mouse_reference_normalized_ClusterID_vargenes_",nVarGenes,"_dims",nDims_integration_anchor,".pdf", sep=""), width=6, height=10)
library(corrplot)
hist(seurat_macaque_coe@meta.data$prediction.score.max, breaks=20)
corrplot(results.norm, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()
mus_names <- c("TBC","IDC","IBC_IPh_HeC","HC","DC_PC","CC_ISC_OSC")
mac_names_anatomical <- c("TBC","IDC","IBC_IPh_HeC","HC","DC_PC","CC_ISC_OSC")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered <- results.norm[order(match(rownames(results.norm),
                                                   mac_names_anatomical_rev)), 
                                       order(match(colnames(results.norm), mus_names))]
#pdf(file=paste0("macaque_query_mouse_reference_normalized_ClusterID_vargenes_",nVarGenes,"_dims",nDims_integration_anchor,"_PCA_re-ordered_anatomical-rev.pdf", sep=""), width=10, height=8)
corrplot(results.norm_reordered, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()
#filtering based on >25% matches
results.norm_reordered_cropped <- results.norm_reordered[rowMaxs(results.norm_reordered) > 0.25, colMaxs(results.norm_reordered) > 0.25]
#pdf(file=paste0("macaque_query_mouse_reference_normalized_ClusterID_vargenes_",nVarGenes,"_dims",nDims_integration_anchor,"_reordered_dendrogram_names_compact-25-rev.pdf", sep=""), width=10, height=8)
corrplot(results.norm_reordered_cropped, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()
################################################################################################

##symphany

Idents(macaque_coe) <- "cell_type"
table(Idents(macaque_coe))
P1=DimPlot(macaque_coe, reduction = "umap",label = TRUE,repel = TRUE)
Idents(mouse_coe) <- "cell_type"
table(Idents(mouse_coe))
P2=DimPlot(mouse_coe, reduction = "umap",label = TRUE,repel = TRUE)
P1+P2
library(symphony)
library(Matrix)
library(ggplot2)
library(Seurat)
library(BGmix)

mouse_coe$age=mouse_coe$batch
ref_exp = mouse_coe@assays[["RNA"]]@counts
ref_metadata = mouse_coe@meta.data
ref_UMAP = mouse_coe@reductions[["umap"]]@cell.embeddings
colnames(ref_UMAP) = c("UMAP1","UMAP2")
ref_metadata = cbind(ref_UMAP,ref_metadata)
colnames(ref_metadata)
ref_metadata = ref_metadata[,c("UMAP1","UMAP2","age","cell_type")]
ref_metadata$Sample = NA
ref_metadata$Sample = "reference"

# Build reference
set.seed(0)
reference = symphony::buildReference(
  ref_exp,
  ref_metadata,
  #vars = c('batch'),         # variables to integrate over
  K = 100,                    # number of Harmony clusters
  verbose = TRUE,             # verbose output
  do_umap = TRUE,             # can set to FALSE if want to run umap separately later
  do_normalize = TRUE,        # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',    # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'batch', # metadata column specifying groups for variable gene selection 
  topn = 2000,                # number of variable genes to choose per group
  d = 20,                     # number of PCs
  save_uwot_path = './testing_uwot_model_2'  
)

testing_uwot_model_2 = uwot::load_uwot("testing_uwot_model_2", verbose = FALSE)
head(testing_uwot_model_2[["embedding"]])
head(mouse_coe@reductions[["umap"]]@cell.embeddings)
colnames(mouse_coe@reductions[["umap"]]@cell.embeddings) = c("UMAP1","UMAP2")
testing_uwot_model_2[["embedding"]] = mouse_coe@reductions[["umap"]]@cell.embeddings

uwot::save_uwot(model = testing_uwot_model_2,file = "testing_uwot_model_2",verbose = FALSE)

head(testing_uwot_model_2[["embedding"]])
head(mouse_coe@reductions[["umap"]]@cell.embeddings)

###########
reference$normalization_method = 'log(CP10k+1)' # optionally save normalization method in custom slot
# Save reference (modify with your desired output path)
#saveRDS(reference, './testing_reference2.rds')

head(reference[["umap"]][["embedding"]])
head(ref_UMAP)
reference[["umap"]][["embedding"]] = ref_UMAP


umap_labels = ref_metadata
#plotBasic(reference,ybar=umap_labels,ss =umap_labels)
ggplot(data=umap_labels,aes(UMAP1,UMAP2,colour=cell_type))+
  geom_point(size=1,alpha=0.7)
#ggsave(file = "reference_umap.jpeg", width= 6, height = 4)
#ggsave(file = "reference_umap.pdf", width= 6, height = 2)



query_exp = macaque_coe@assays[["RNA"]]@counts
query_meta = macaque_coe@meta.data

table(query_meta$cell_type)
query_exp[1:5,1:5]
query_meta = query_meta[,c("age","cell_type")]


query_meta$Sample = NA
query_meta$Sample = "query"
colnames(ref_metadata)
colnames(query_meta) 
colnames(query_meta) = c("age","cell_type","Sample")
# Map query
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_meta,            # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = TRUE,   # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP

query = knnPredict(query, reference, reference$meta_data$cell_type, k = 5)

#################
table(query$meta_data$cell_type_pred_knn)
reference$meta_data$cell_type_pred_knn = NA
reference$meta_data$cell_type_pred_knn_prob = NA
umap=query$meta_data
umap=cbind(query$umap,umap)
query$meta_data = umap
colnames(reference[["meta_data"]])[1:2] = c("UMAP1","UMAP2")
meta_data_combined = rbind(query$meta_data,reference$meta_data)

#################

p1=ggplot(data=query$meta_data,aes(UMAP1,UMAP2,colour=cell_type))+
  geom_point(size=2,alpha=0.7)
#ggsave(file = "query_umap.jpeg", width= 6, height = 4)
#ggsave(file = "query_umap.pdf", width= 6, height = 4)

p2=ggplot(data=query$meta_data,aes(UMAP1,UMAP2,colour=query$meta_data$cell_type_pred_knn))+
  geom_point(size=2,alpha=0.7)
#ggsave(file = "query_umap_pred.jpeg", width= 6, height = 4)
#ggsave(file = "query_umap_pred.pdf", width= 6, height = 4)
p1+p2

########################################

table(query$meta_data$cell_type_pred_knn)
table(query$meta_data$cell_type)
plot_heatmap = query$meta_data
write.csv(plot_heatmap,file = "plot_heatmap_integrated_symphony_v3.csv")
results <- table(query$meta_data$cell_type,query$meta_data$cell_type_pred_knn)


results.norm <- 100*(results/rowSums(results))
mus.ids <- colnames(results.norm)
results.norm <- results.norm[,mus.ids]
library(corrplot)
hist(coembed@meta.data$prediction.score.max, breaks=20)
corrplot(results.norm, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()


mus_names <- c("TBC","IDC","IBC_IPh_HeC","HC","DC_PC","CC_ISC_OSC")
mac_names_anatomical <- c("TBC","IDC","IBC_IPh_HeC","HC","DC_PC","CC_ISC_OSC")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered <- results.norm[order(match(rownames(results.norm),
                                                   mac_names_anatomical_rev)), 
                                       order(match(colnames(results.norm), mus_names))]
#pdf(file=paste0("macaque_query_mouse_reference_normalized_ClusterID_vargenes_",nVarGenes,"_dims",nDims_integration_anchor,"_PCA_re-ordered_anatomical-rev.pdf", sep=""), width=10, height=8)
corrplot(results.norm_reordered, order="original",tl.pos="lt", method="color", tl.col="black", is.corr=F)
dev.off()


##################################################
##############################################
## Calculate intersection between macaque  vs mouse DE genes for each cell type according to Jaccard similarity index 
Idents(macaque_coe) <- "cell_type"
Idents(mouse_coe) <- "cell_type"
macaque_coe<- PrepSCTFindMarkers(macaque_coe)
mouse_coe<- PrepSCTFindMarkers(mouse_coe)
deg_macaque <- FindAllMarkers(macaque_coe, assay = "SCT", slot = "data",
                              test.use = "roc")

deg_mouse <- FindAllMarkers(mouse_coe, assay = "SCT", slot = "data",
                            test.use = "roc")
deg_macaque_TOP100 <- deg_macaque %>% group_by(cluster) %>% top_n(n=100,wt=myAUC)
deg_mouse_TOP100 <- deg_mouse %>% group_by(cluster) %>% top_n(n=100,wt=myAUC)

##we also used scanpy-based tt-test for calculating the DEGs within each cell type from macaque species or mouse species, respectively.

deg_macaque_2 <- read.csv(file = "DEG_MACAQUE_COE_JACARD_INDEX.csv")
deg_mouse_2 <- read.csv(file = "DEG_MOUSE_COE_JACARD_INDEX.csv")





library(tidyverse)
library(readxl)
library(data.table)
library (ggplot2)
#install.packages("ggpubr")
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

cluster1=names(table(deg_macaque$cluster))
cluster2=names(table(deg_mouse$cluster))

## jaccard
df_ja=c()

for (i in cluster1) {
  ja=c()
  for (j in cluster2) {
    a=deg_macaque[deg_macaque$cluster==i,]
    a=a$gene
    b=deg_mouse[deg_mouse$cluster==j,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    ja=c(ja,jaccard)
  }
  df_ja=rbind(df_ja,ja)
}

rownames(df_ja)=paste0(cluster1)
colnames(df_ja)=paste0(cluster2)
library(pheatmap)
pheatmap(df_ja,cluster_cols = T,show_rownames = T, show_colnames =T,
         color = colorRampPalette(c(rep("blue",1), "white", rep("red",1)))(100),
         breaks = seq(0,0.25,length.out = 100),border_color = NA)

library(RColorBrewer)
## Set the common color scale
breaksList = seq(0, 0.2, by = 0.00001)


heatmap_up_noleg = pheatmap(mat = df_ja, 
                            cluster_rows = T,
                            cluster_cols = T,
                            display_numbers = F,
                            number_format = "%.1f",
                            scale = "none",
                            na_col = "#DDDDDD", 
                            border_color = NA, 
                            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
                            breaks = breaksList, legend = FALSE)


heatmap_up_noleg = pheatmap(mat = df_ja, 
                            cluster_rows = F,
                            cluster_cols = F,
                            display_numbers = F,
                            number_format = "%.1f",
                            scale = "none",
                            na_col = "#DDDDDD", 
                            border_color = NA, 
                            color = colorRampPalette(rev(brewer.pal(n = 7,
                                                                    name = "RdYlBu")))(length(breaksList)), 
                            breaks = breaksList, legend = FALSE)

#cumulative gene family expression by subclasses
library(Seurat)
library(matrixStats)
library(ggplot2)
Seurat_coe_LABEL_TRNSFER_v2 <- readRDS("F:/PROJECTS/PROJECT_MONKEY SC-RNA SEQ/WORKPLACE/R/Seurat_coe_LABEL_TRNSFER_v2.rds")
Seurat_coe_LABEL_TRNSFER=Seurat_coe_LABEL_TRNSFER_v2
DimPlot(Seurat_coe_LABEL_TRNSFER,group.by = "cell_type")

Seurat_coe_LABEL_TRNSFER$subclasses.species <- paste(Seurat_coe_LABEL_TRNSFER$cell_type,
                                                     Seurat_coe_LABEL_TRNSFER$species,sep = ".")
Idents(Seurat_coe_LABEL_TRNSFER) <- "subclasses.species"
DimPlot(Seurat_coe_LABEL_TRNSFER, reduction = "umap",label = TRUE,split.by = c("species"))
Idents(Seurat_coe_LABEL_TRNSFER) <- "cell_type"
DefaultAssay(Seurat_coe_LABEL_TRNSFER) <- "SCT"
object <- Seurat_coe_LABEL_TRNSFER
table(Idents(Seurat_coe_LABEL_TRNSFER))
clusters <-levels(Idents(Seurat_coe_LABEL_TRNSFER))
genes.use = rownames(object)
object.raw.data <- as.matrix(GetAssayData(object, slot="data",assay = "SCT"))
pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
rownames(pct.matrix) <- genes.use
colnames(pct.matrix) <- clusters
thresh.min=0
for (i in clusters){
  cells.cluster <- WhichCells(object=object, idents=i)
  data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
  pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  pct.matrix[,i] <- pct.cluster
}
pct.df <- as.data.frame(pct.matrix)
Human_TFs <- read.csv("Human_TFs.csv", header=T, row.names = 1)
Human_TFs$gene_family <- "transcription factors"
Human_TFs <- Human_TFs[,c(1,3,2)]
names(Human_TFs) <- c("ncbi_gene_symbol","gene_family","TF_family")

Human_GPCRs <- read.delim("GPCRs_gene_symbols.txt", header=T, sep = "\t")
head(Human_GPCRs)
Human_GPCRs <- Human_GPCRs[,c(1:2)]
Human_GPCRs$gene_family <- "GPCR"
Human_GPCRs$TF_family <- NA

Human_ICs <- read.delim("ion-channels_gene_symbols.txt", header=T, sep = "\t")
head(Human_ICs)
Human_ICs <- Human_ICs[,c(1:2)]
Human_ICs$gene_family <- "ion channels"
Human_ICs$TF_family <- NA

Human_CAMs <- read.delim("CAMs_gene_symbols.txt", header=T, sep = "\t")
head(Human_CAMs)
Human_CAMs <- Human_CAMs[,c(1:2)]
Human_CAMs$gene_family <- "cell adhesion molecules"
Human_CAMs$TF_family <- NA

Human_PGs <- read.delim("proteoglycans_gene_symbols.txt", header=T, sep = "\t")
head(Human_PGs)
Human_PGs <- Human_PGs[,c(1:2)]
Human_PGs$gene_family <- "proteoglycans"
Human_PGs$TF_family <- NA

Human_ribo <- read.delim("ribosome_gene_symbols.txt", header=T, sep = "\t")
head(Human_ribo)
Human_ribo <- Human_ribo[,c(1:2)]
Human_ribo$gene_family <- "ribosome"
Human_ribo$TF_family <- NA

#transcription factors
TF_intersection <- intersect(rownames(pct.df),Human_TFs$ncbi_gene_symbol)
pct.df_TFs <- pct.df[rownames(pct.df) %in% TF_intersection,]
Human_TFs_cropped <- Human_TFs[Human_TFs$ncbi_gene_symbol %in% TF_intersection,]

Human_TFs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_TFs)),ncol=2)
Human_TFs_cluster_counts <- as.data.frame(Human_TFs_cluster_counts)
names(Human_TFs_cluster_counts) <- c("gene","cluster_count")
for (i in 1:length(rownames(pct.df_TFs))){
  gene_row <- pct.df_TFs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_TFs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_TFs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_TFs_cluster_counts_matched <- Human_TFs_cropped[order(match(Human_TFs_cropped$ncbi_gene_symbol, 
                                                                  Human_TFs_cluster_counts$gene)),]
Human_TFs_cluster_counts$gene_family <- Human_TFs_cluster_counts_matched$gene_family
Human_TFs_cluster_counts$TF_family <- Human_TFs_cluster_counts_matched$TF_family
TF_counts <- Human_TFs_cluster_counts


#GPCRs
GPCR_intersection <- intersect(rownames(pct.df),Human_GPCRs$ncbi_gene_symbol)
pct.df_GPCRs <- pct.df[rownames(pct.df) %in% GPCR_intersection,]
Human_GPCRs_cropped <- Human_GPCRs[Human_GPCRs$ncbi_gene_symbol %in% GPCR_intersection,]

Human_GPCRs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_GPCRs)),ncol=2)
Human_GPCRs_cluster_counts <- as.data.frame(Human_GPCRs_cluster_counts)
names(Human_GPCRs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_GPCRs))){
  gene_row <- pct.df_GPCRs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_GPCRs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_GPCRs_cluster_counts$cluster_count[i] <- cluster_count
}
Human_GPCRs_cluster_counts_matched <- Human_GPCRs_cropped[order(match(Human_GPCRs_cropped$ncbi_gene_symbol, 
                                                                      Human_GPCRs_cluster_counts$gene)),]
Human_GPCRs_cluster_counts$gene_family <- Human_GPCRs_cluster_counts_matched$gene_family
GPCR_counts <- Human_GPCRs_cluster_counts

#ion channels
IC_intersection <- intersect(rownames(pct.df),Human_ICs$ncbi_gene_symbol)
pct.df_ICs <- pct.df[rownames(pct.df) %in% IC_intersection,]
Human_ICs_cropped <- Human_ICs[Human_ICs$ncbi_gene_symbol %in% IC_intersection,]

Human_ICs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_ICs)),ncol=2)
Human_ICs_cluster_counts <- as.data.frame(Human_ICs_cluster_counts)
names(Human_ICs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ICs))){
  gene_row <- pct.df_ICs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_ICs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_ICs_cluster_counts$cluster_count[i] <- cluster_count
}


Human_ICs_cluster_counts_matched <- Human_ICs_cropped[order(match(Human_ICs_cropped$ncbi_gene_symbol,
                                                                  Human_ICs_cluster_counts$gene)),]
Human_ICs_cluster_counts$gene_family <- Human_ICs_cluster_counts_matched$gene_family
IC_counts <- Human_ICs_cluster_counts

#cell adhesion molecules
CAM_intersection <- intersect(rownames(pct.df),Human_CAMs$ncbi_gene_symbol)
pct.df_CAMs <- pct.df[rownames(pct.df) %in% CAM_intersection,]
Human_CAMs_cropped <- Human_CAMs[Human_CAMs$ncbi_gene_symbol %in% CAM_intersection,]

Human_CAMs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_CAMs)),ncol=2)
Human_CAMs_cluster_counts <- as.data.frame(Human_CAMs_cluster_counts)
names(Human_CAMs_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_CAMs))){
  gene_row <- pct.df_CAMs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_CAMs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_CAMs_cluster_counts$cluster_count[i] <- cluster_count
}

Human_CAMs_cluster_counts_matched <- Human_CAMs_cropped[order(match(Human_CAMs_cropped$ncbi_gene_symbol, 
                                                                    Human_CAMs_cluster_counts$gene)),]
Human_CAMs_cluster_counts$gene_family <- Human_CAMs_cluster_counts_matched$gene_family
CAM_counts <- Human_CAMs_cluster_counts

#proteoglycans
PG_intersection <- intersect(rownames(pct.df),Human_PGs$ncbi_gene_symbol)
pct.df_PGs <- pct.df[rownames(pct.df) %in% PG_intersection,]
Human_PGs_cropped <- Human_PGs[Human_PGs$ncbi_gene_symbol %in% PG_intersection,]

Human_PGs_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_PGs)),ncol=2)
Human_PGs_cluster_counts <- as.data.frame(Human_PGs_cluster_counts)
names(Human_PGs_cluster_counts) <- c("gene","cluster_count")
for (i in 1:length(rownames(pct.df_PGs))){
  gene_row <- pct.df_PGs[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_PGs_cluster_counts$gene[i] <- rownames(gene_row)
  Human_PGs_cluster_counts$cluster_count[i] <- cluster_count
}


Human_PGs_cluster_counts_matched <- Human_PGs_cropped[order(match(Human_PGs_cropped$ncbi_gene_symbol, 
                                                                  Human_PGs_cluster_counts$gene)),]
Human_PGs_cluster_counts$gene_family <- Human_PGs_cluster_counts_matched$gene_family
PG_counts <- Human_PGs_cluster_counts

#ribosomes
ribo_intersection <- intersect(rownames(pct.df),Human_ribo$ncbi_gene_symbol)
pct.df_ribos <- pct.df[rownames(pct.df) %in% ribo_intersection,]
Human_ribos_cropped <- Human_ribo[Human_ribo$ncbi_gene_symbol %in% ribo_intersection,]
Human_ribos_cropped <- Human_ribos_cropped[unique(Human_ribos_cropped$ncbi_gene_symbol),]

Human_ribos_cluster_counts <- matrix(data=NA,nrow=length(rownames(pct.df_ribos)),ncol=2)
Human_ribos_cluster_counts <- as.data.frame(Human_ribos_cluster_counts)
names(Human_ribos_cluster_counts) <- c("gene","cluster_count")

for (i in 1:length(rownames(pct.df_ribos))){
  gene_row <- pct.df_ribos[i,]
  cluster_count <- length(gene_row[gene_row > 0.2])
  Human_ribos_cluster_counts$gene[i] <- rownames(gene_row)
  Human_ribos_cluster_counts$cluster_count[i] <- cluster_count
}


Human_ribos_cluster_counts_matched <- Human_ribos_cropped[order(match(Human_ribos_cropped$ncbi_gene_symbol, 
                                                                      Human_ribos_cluster_counts$gene)),]
Human_ribos_cluster_counts$gene_family <- Human_ribos_cluster_counts_matched$gene_family
Human_ribos_cluster_counts$ribo_family <- Human_ribos_cluster_counts_matched$ribo_family
ribo_counts <- Human_ribos_cluster_counts

#Homeodomain 
countsHD <- TF_counts[TF_counts$TF_family=="Homeodomain",]
countsHD <- countsHD[order(countsHD$cluster_count, decreasing = F),]
countsHD <- countsHD[countsHD$cluster_count>0,]
tHD <- table(countsHD$cluster_count)
xHD <- as.numeric(names(tHD))
yHD <- as.numeric(tHD)
plot(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.5, type="o", pch=19)
#bHLH 
countsbHLH <- TF_counts[TF_counts$TF_family=="bHLH",]
countsbHLH <- countsbHLH[order(countsbHLH$cluster_count, decreasing = F),]
countsbHLH <- countsbHLH[countsbHLH$cluster_count>0,]
tbHLH <- table(countsbHLH$cluster_count)
xbHLH <- as.numeric(names(tbHLH))
ybHLH <- as.numeric(tbHLH)
plot(xbHLH,(cumsum(ybHLH)/max(cumsum(ybHLH))), cex=0.5, type="o", pch=19)
#bZIP 
countsZIP <- TF_counts[TF_counts$TF_family=="bZIP",]
countsZIP <- countsZIP[order(countsZIP$cluster_count, decreasing = F),]
countsZIP <- countsZIP[countsZIP$cluster_count>0,]
tZIP <- table(countsZIP$cluster_count)
xZIP <- as.numeric(names(tZIP))
yZIP <- as.numeric(tZIP)
plot(xZIP,(cumsum(yZIP)/max(cumsum(yZIP))), cex=0.5, type="o", pch=19)

#Forkhead 
countsFH <- TF_counts[TF_counts$TF_family=="Forkhead",]
countsFH <- countsFH[order(countsFH$cluster_count, decreasing = F),]
countsFH <- countsFH[countsFH$cluster_count>0,]
tFH <- table(countsFH$cluster_count)
xFH <- as.numeric(names(tFH))
yFH <- as.numeric(tFH)
plot(xFH,(cumsum(yFH)/max(cumsum(yFH))), cex=0.5, type="o", pch=19)

#all TFs 
counts <- TF_counts[order(TF_counts$cluster_count, decreasing = F),]
countsTFs <- counts[counts$cluster_count>0,]
tTFs <- table(countsTFs$cluster_count)
xTFs <- as.numeric(names(tTFs))
yTFs <- as.numeric(tTFs)
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.5, type="o", pch=19)

#Nuclear receptor 
countsNR <- TF_counts[TF_counts$TF_family=="Nuclear receptor",]
countsNR <- countsNR[order(countsNR$cluster_count, decreasing = F),]
countsNR <- countsNR[countsNR$cluster_count>0,]
tNR <- table(countsNR$cluster_count)
xNR <- as.numeric(names(tNR))
yNR <- as.numeric(tNR)
plot(xNR,(cumsum(yNR)/max(cumsum(yNR))), cex=0.5, type="o", pch=19)

#C2H2 ZF 
countsC2F <- TF_counts[TF_counts$TF_family=="C2H2 ZF",]
countsC2F <- countsC2F[order(countsC2F$cluster_count, decreasing = F),]
countsC2F <- countsC2F[countsC2F$cluster_count>0,]
tC2F <- table(countsC2F$cluster_count)
xC2F <- as.numeric(names(tC2F))
yC2F <- as.numeric(tC2F)
plot(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.5, type="o", pch=19)

plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=0.5, type="o", pch=19, xlab = "Number of cochlear epithelium cell types expressing", ylab = "Cumulative fraction")
lines(xNR,(cumsum(yNR)/max(cumsum(yNR))), cex=0.5, type="o", pch=19, col="magenta")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=0.5, type="o", pch=19, col="yellow")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=0.5, type="o", pch=19, col="darkgreen")
lines(xbHLH,(cumsum(ybHLH)/max(cumsum(ybHLH))), cex=0.5, type="o", pch=19, col="blue")
lines(xFH,(cumsum(yFH)/max(cumsum(yFH))), cex=0.5, type="o", pch=19, col="red")
lines(xZIP,(cumsum(yZIP)/max(cumsum(yZIP))), cex=0.5, type="o", pch=19, col="orange")
#all ion channels 
counts_ICs <- IC_counts[order(IC_counts$cluster_count, decreasing = F),]
countsICs <- counts_ICs[counts_ICs$cluster_count>0,]
tICs <- table(countsICs$cluster_count)
xICs <- as.numeric(names(tICs))
yICs <- as.numeric(tICs)
plot(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=0.5, type="o", pch=19)

#all GPCRs 
counts_GPCRs <- GPCR_counts[order(GPCR_counts$cluster_count, decreasing = F),]
countsGPCRs <- counts_GPCRs[counts_GPCRs$cluster_count>0,]
tGPCRs <- table(countsGPCRs$cluster_count)
xGPCRs <- as.numeric(names(tGPCRs))
yGPCRs <- as.numeric(tGPCRs)
plot(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=0.5, type="o", pch=19)

#all CAMs 
counts_CAMs <- CAM_counts[order(CAM_counts$cluster_count, decreasing = F),]
countsCAMs <- counts_CAMs[counts_CAMs$cluster_count>0,]
tCAMs <- table(countsCAMs$cluster_count)
xCAMs <- as.numeric(names(tCAMs))
yCAMs <- as.numeric(tCAMs)
plot(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=0.5, type="o", pch=19)

#all proteoglycans
counts_PGs <- PG_counts[order(PG_counts$cluster_count, decreasing = F),]
countsPGs <- counts_PGs[counts_PGs$cluster_count>0,]
tPGs <- table(countsPGs$cluster_count)
xPGs <- as.numeric(names(tPGs))
yPGs <- as.numeric(tPGs)
plot(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=0.5, type="o", pch=19)
#all ribosomal genes 
counts_ribos <- ribo_counts[order(ribo_counts$cluster_count, decreasing = F),]
countsribos <- counts_ribos[counts_ribos$cluster_count>0,]
tRibos <- table(countsribos$cluster_count)
xRibos <- as.numeric(names(tRibos))
yRibos <- as.numeric(tRibos)
plot(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=1, type="o", pch=19)
plot(xTFs,(cumsum(yTFs)/max(cumsum(yTFs))), cex=1, type="o", pch=19, 
     xlab = "Number of neuron types expressing",
     ylab = "Cumulative fraction",xlim=c(0,6),ylim=c(0,1), col="black")
lines(xC2F,(cumsum(yC2F)/max(cumsum(yC2F))), cex=1, type="o", pch=19, col="magenta")
lines(xHD,(cumsum(yHD)/max(cumsum(yHD))), cex=1, type="o", pch=19, col="darkgreen")
lines(xGPCRs,(cumsum(yGPCRs)/max(cumsum(yGPCRs))), cex=1, type="o", pch=19, col="blue")
lines(xICs,(cumsum(yICs)/max(cumsum(yICs))), cex=1, type="o", pch=19, col="red")
lines(xCAMs,(cumsum(yCAMs)/max(cumsum(yCAMs))), cex=1, type="o", pch=19, col="orange")
lines(xPGs,(cumsum(yPGs)/max(cumsum(yPGs))), cex=1, type="o", pch=19, col="darkblue")
lines(xRibos,(cumsum(yRibos)/max(cumsum(yRibos))), cex=1, type="o", pch=19, col="darkred")

legend("bottomright", legend=c("all TFs (292)", "zinc finger-C2H2 (89)","homeodomain (14)","GPCRs (11)","ion channels (66)","cell adhesion molecules (50)","proteoglycans (17)","ribosomal genes (40)"),
       col=c("black", "magenta","darkgreen","blue","red","orange","darkblue","darkred"), lty=1, cex=0.4)
####################################################################################################
#METANEIGHBOR ANALYSIS
library(Seurat)
library(gridExtra)
library(ggplot2)
library(sctransform)
library(tidyverse)
library(dplyr)
library(Matrix)
library(matrixStats)
library(gplots)
library(ggplot2)
library(feather)
library(SingleCellExperiment)
options(stringsAsFactors = FALSE)
sce=Seurat_coe_LABEL_TRNSFER
################## Prep SingleCellExperiment object for analysis  #############################################
library(SingleCellExperiment)
library(Matrix)
sce$class_label = "CoE"
table(sce$class_label)
Idents(sce) <- "class_label"
sce$class.species <- paste(sce$class_label, sce$species, sep = "_")
table(sce$class.species)
Idents(sce) <- "species"
table(sce$species)
mac <- subset(sce, idents ="Macaque")
ms <- subset(sce, idents ="Mouse")

mac <- GetAssayData(object = mac[["RNA"]], slot = "counts")
ms <- GetAssayData(object = ms[["RNA"]], slot = "counts")
gs=read.csv("orthologTable_human_macaque_mouse.csv",header=T,stringsAsFactors=F)  ### read in ortholog mapping table


mac <- as.matrix(mac)
ms <- as.matrix(ms)

p2 <- subset(sce, idents ="Macaque")
p3<- subset(sce, idents ="Mouse")

p2$sample_id <- colnames(p2)
p3$sample_id <- colnames(p3)

p2$study_id="Macaque"
p3$study_id="Mouse"
p2<- as.matrix(p2@meta.data)
p3<- as.matrix(p3@meta.data)

m<-match(as.character(rownames(mac)),as.character(gs$macaque_symbol))  ### match gene IDs from expression data to those in the ortholog table
f.a=!is.na(m)
f.b=m[f.a]
mac=mac[f.a,]
rownames(mac)=gs[f.b,"human_symbol"]  ## convert to human IDs (16K genes)
m<-match(as.character(rownames(ms)),as.character(gs$human_symbol))
f.a=!is.na(m)
f.b=m[f.a]
ms=ms[f.a,]
rownames(ms)=as.character(gs[f.b,"human_symbol"]) ## 14350 genes
rownames(p3)=colnames(ms)
m<-match(as.character(rownames(mac)),as.character(rownames(ms)))
f.a=!is.na(m)
f.b=m[f.a]

dat=cbind(mac[f.a,],ms[f.b,])
x=as.vector(rownames(dat))
rownames(dat)=NULL
rownames(dat)=x
p_pri2=rbind(p2,p3)
rownames(p_pri2)=colnames(dat)
sce_all=SingleCellExperiment(assays=list(counts=dat),colData=p_pri2)  
sce_all

saveRDS(sce_all,file="mac_ms_sce_V3.rds")

sce_all$sample_id_append <- sce_all$sample_id
head(sce_all$sample_id_append)
sce$sample_id <-colnames(sce)
head(sce$sample_id)
m<-match(sce_all$sample_id_append,sce$sample_id)
sum(!is.na(m))
f.a=!is.na(m)
f.b=m[f.a]

sce_all$cell_type[f.a]=as.character(sce$cell_type[f.b])
#sce_all$final_integrated_cluster_color[f.a]=cochlea_epithelium$[f.b]
head(sce_all$study_id)
head(sce_all$cell_type[f.a])
sce_all=sce_all[,f.a]
head(sce_all$study_id)
saveRDS(sce_all,file="mac_ms_sce_v3.rds")
sce_all
########## Within- and cross-species classification with highly variable genes ###################

library(SingleCellExperiment)
library(Matrix)
source("1v1_analysis.R")
classes=unique(sce_all$class_label)
#classes <- classes[-1]
vgs=vector("list",length=length(classes))
names(vgs)=names(classes)
for(i in seq_along(classes)){
  f=sce_all$class_label==classes[i]
  vgs[[i]]=get_variable_genes(sce_all[,f.b])
}

mn_1v1_subclass=vector("list",length=length(classes))
for(i in seq_along(classes)){
  f=sce_all$class_label==classes[i]
  
  mn_1v1_subclass[[i]]=compute_best_hits(sce_all[vgs[[i]],f.b],
                                         sce_all$cell_type[f.b],one_vs_one=TRUE)
}
saveRDS(mn_1v1_subclass,file="coe_1v1_subclass_level_v3.Rdata")

b<- as.data.frame(mn_1v1_subclass)
write.csv(b,file = "mn_1v1_coe_subclasses_v3.csv")
options(stringsAsFactors = FALSE)

d<- read.csv("mn_1v1_coe_subclasses_v3.csv",row.names = 1)
library(pheatmap)
library(viridis)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE)
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = viridis(50))
#color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = cividis(50))

########## Within- and cross-species classification using HGNC and SynGO gene sets ###################
library(SingleCellExperiment)
library(Matrix)
options(stringsAsFactors = FALSE)
source("metaneighbor.R")
gs=readRDS("hgnc_syngo.rds")
################ within mouse
head(sce_all$study_id)
spe_f=sce_all$study_id=="Mouse"
sce_ms=sce_all[,spe_f]

rand3cv=sample(c(1:3),size=sum(spe_f),replace=T)
sce_ms$study_id=rand3cv

classes=unique(sce_ms$class_label)
res=vector("list",length=length(classes))
names(res)=classes
head(sce_ms$cell_type)
for(i in seq_along(classes)){
  sub_f=sce_ms$class_label==classes[i]
  x=model.matrix(~(as.character(sce_ms$cell_type[sub_f]))+0)
  colnames(x)=names(table(sce_ms$cell_type[sub_f]))
  res[[i]]=call_my_metaneighbor(sce_ms[,sub_f],gs,x)
}

saveRDS(res,file="within_mouse_MN_cluster_level_v3.rds")
res_mouse=res
rm(sce_ms)
################ within macaque
spe_f=sce_all$study_id=="Macaque"
sce_mac=sce_all[,spe_f]
rand3cv=sample(c(1:3),size=sum(spe_f),replace=T)
sce_mac$study_id=rand3cv
head(sce_mac$study_id)
classes=unique(sce_mac$class_label)
res=vector("list",length=length(classes))
names(res)=classes
head(sce_mac$cell_type,n=20)
for(i in seq_along(classes)){
  sub_f=sce_mac$class_label==classes[i]
  x=model.matrix(~(as.character(sce_mac$cell_type[sub_f]))+0)
  colnames(x)=names(table(sce_mac$cell_type[sub_f]))
  res[[i]]=call_my_metaneighbor(sce_mac[,sub_f],gs,x)
}
saveRDS(res,file="within_macaque_MN_cluster_level_v3.rds")
res_macaque=res
rm(sce_mac)
################ cross macaque and mouse
sce_all2=sce_all
sce_all2$study_original=sce_all2$study_id
classes=unique(sce_all2$class_label)
res=vector("list",length=length(classes))
names(res)=classes

for(i in seq_along(classes)){
  sub_f=sce_all2$class_label==classes[i]
  x=model.matrix(~(as.character(sce_all2$cell_type[sub_f]))+0)
  colnames(x)=names(table(sce_all2$cell_type[sub_f]))
  res[[i]]=call_my_metaneighbor(sce_all2[,sub_f],gs,x)
}


save(res,file="cross_macaque_mouse_cluster_level2_v3.rds")
res_macaque_mouse=res


write.csv(res_macaque,file = "1within_macaque_mn_cluster_level_v3.csv")
write.csv(res_mouse,file = "1within_mouse_mn_cluster_level_v3.csv")
write.csv(res_macaque_mouse,file = "1cross_macaque_mouse_mn_cluster_level_v3.csv")



b <- read.csv("1within_macaque_mn_cluster_level_v3.csv")
b
b<- aggregate(b$CoE.auroc,by=list(type=b$CoE.gene_set),mean)
c <- read.csv("1within_mouse_mn_cluster_level_v3.csv")
c
c<- aggregate(c$CoE.auroc,by=list(type=c$CoE.gene_set),mean)

f <- read.csv("1cross_macaque_mouse_mn_cluster_level_v3.csv")
f
f<- aggregate(f$CoE.auroc,by=list(type=f$CoE.gene_set),mean)

g <- merge(x=b,y=c,by="type")
g <- merge(x=g,y=f,by="type")
write.csv(g, file="MN_CLUSTER_meanROC_coe_v3.csv")
g<-read.csv("MN_CLUSTER_meanROC_coe_v3.csv")
# Library
library(ggplot2)
library(hrbrthemes)
p3 <- ggplot(g, aes(g$within_species_meanROC,g$macaque_meanROC)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()+
  xlim(0.4,1)+
  ylim(0.4,1)
p3
p4 <- ggplot(g, aes(g$within_species_meanROC,g$mouse_meanROC)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()+
  xlim(0.4,1)+
  ylim(0.4,1)
p4
p7 <- ggplot(g, aes(g$within_species_meanROC,g$macaque_mouse_meanROC)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=FALSE) +
  theme_ipsum()+xlim(0.4,1)+
  ylim(0.4,1)
p7
coef(lm(g$macaque_mouse_meanROC~g$within_species_meanROC))[2]
g<-read.csv("MN_CLUSTER_meanROC_coe_v3.csv")
h<-read.csv("Supplementary Table 9_2.csv")
i<- merge(x=g,y=h,by="gene_set_label")
j<- i[!duplicated(i$gene_set_label), ]
k<- merge(x=j,y=g, all = TRUE)
write.csv(k,file = "coe_mn_cluster_meanROC_v3.csv")

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
mn.df <- read.csv("coe_mn_cluster_meanROC_v3.csv")
mn.df$gene_set_type <- factor(mn.df$gene_set_type, 
                              levels = c("Other", "Signaling", 
                                         "Cell Adhesion", "Ion Channel"))
roc.tests <- c("macaque_meanROC", "mouse_meanROC", 
               "macaque_mouse_meanROC")
mn.l <- mn.df %>% 
  subset(cell_class == "CoE") %>%
  arrange(gene_set_type) %>% 
  gather(roc.tests, key = "comp", value = "ROC")
mn.l$comp <- factor(mn.l$comp, levels = roc.tests)
levels(mn.l$comp) <- c("Macaque", "Mouse", 
                       
                       "Macaque vs. Mouse")
mn.l$comp_type <- ifelse(mn.l$comp %in% c("Macaque", "Mouse"), 
                         "Within-species", "Cross-species")
geneset.pal <- c("#f70ff1", "#2626d3", "#bbce32", "#4ce03d","#998a82")
g.scatter <- ggplot(mn.l, aes(x = within_species_meanROC, y = ROC)) +
  facet_grid(. ~ comp) +
  geom_hline(yintercept = 0.5, color = "grey", size = 0.25) +
  geom_vline(xintercept = 0.5, color = "grey", size = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 0.25) +
  geom_point(alpha = 0.5, aes(color = gene_set_type)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_colour_manual(values=geneset.pal) +
  coord_fixed(xlim = c(0.3, 1), ylim = c(0.3, 1)) +
  scale_x_continuous(breaks = c(0.3, 0.5,0.7, 0.9)) +
  scale_y_continuous(breaks = c(0.3, 0.5,0.7, 0.9)) +
  xlab("Within-species mean AUROC") +
  ylab("AUROC") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot(g.scatter)
mn.l %>% 
  group_by(comp) %>%
  group_modify(~ broom::tidy(lm(ROC ~ within_species_meanROC, data = .x)))

#######################################################################################
##########Gene ontology analysis of conserved cluster markers
#For the gene ontology analysis of maker genes we first calculated conserved marker genes for 
#label-transfer clustersusing Seurat’s FindConservedMarkers 
#(grouping.var = "species", assay="SCT", only.pos= TRUE). We removed genes that were not 
#detected in at least 20% of a single cluster in both species from the marker genes list. We used 
#gprofiler2 (69) to calculate enriched GO-terms (gost: user_threshold = 0.01, correction_method = 
#"g_SCS"), as background we used the intersection of expressed one-to-one orthologs between 
#mouse and macaque
table(sce$cell_type,sce$species)
Idents(sce) <- "cell_type"
# for hc
conserved_markers_hc_go <- FindConservedMarkers(sce,ident.1 = "HC",
                                                ident.2=NULL,grouping.var = "species",
                                                assay="RNA", only.pos= TRUE)
NoisyGenes <- function(object, min.pct=0.2, clusters.use) {
  clusters <- clusters.use
  genes.use = rownames(object)
  object.raw.data <- as.matrix(GetAssayData(object, slot="counts"))
  pct.matrix = matrix(data=NA, nrow=length(genes.use), ncol=length(clusters))
  rownames(pct.matrix) <- genes.use
  colnames(pct.matrix) <- clusters
  thresh.min=0
  for (i in clusters){
    cells.cluster <- WhichCells(object=object, idents=i)
    data.cluster <- object.raw.data[,colnames(object.raw.data) %in% cells.cluster]
    pct.cluster <- round(apply(object.raw.data[genes.use, cells.cluster, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
    pct.matrix[,i] <- pct.cluster
  }
  pct.max <- rowMaxs(pct.matrix)
  names(pct.max) <- genes.use
  noisy.genes <- names(pct.max[pct.max < min.pct])
  return(noisy.genes)
}
noisy.liz <- NoisyGenes(sce, 0.1, levels(Idents(sce)))
conserved_markers_hc_go_filtered <- conserved_markers_hc_go[!rownames(conserved_markers_hc_go) %in% noisy.liz,]
# for dc_pc
conserved_markers_dc_pc_go <- FindConservedMarkers(sce,ident.1 = "DC_PC",
                                                ident.2=NULL,grouping.var = "species",
                                                assay="RNA", only.pos= TRUE)

conserved_markers_dc_pc_filtered <- conserved_markers_dc_pc_go[!rownames(conserved_markers_dc_pc_go) %in% noisy.liz,]

# for IPHIBC
conserved_markers_iphibc_go <- FindConservedMarkers(sce,ident.1 = "IBC_IPh_HeC",
                                                    ident.2=NULL,grouping.var = "species",
                                                    assay="RNA", only.pos= TRUE)

conserved_markers_iphibc_filtered <- conserved_markers_iphibc_go[!rownames(conserved_markers_iphibc_go) %in% noisy.liz,]

# for TBC
conserved_markers_tbc_go <- FindConservedMarkers(sce,ident.1 = "TBC",
                                                ident.2=NULL,grouping.var = "species",
                                                assay="RNA", only.pos= TRUE)

conserved_markers_tbc_filtered <- conserved_markers_tbc_go[!rownames(conserved_markers_tbc_go) %in% noisy.liz,]

# for idc
conserved_markers_idc_go <- FindConservedMarkers(sce,ident.1 = "IDC",
                                                 ident.2=NULL,grouping.var = "species",
                                                 assay="RNA", only.pos= TRUE)

conserved_markers_idc_filtered <- conserved_markers_idc_go[!rownames(conserved_markers_idc_go) %in% noisy.liz,]

# for CC_ISC_OSC
conserved_markers_cc_go <- FindConservedMarkers(sce,ident.1 = "CC_ISC_OSC",
                                                ident.2=NULL,grouping.var = "species",
                                                assay="RNA", only.pos= TRUE)

conserved_markers_cc_filtered <- conserved_markers_cc_go[!rownames(conserved_markers_cc_go) %in% noisy.liz,]

# 
a<- rownames(conserved_markers_hc_go)
b<- rownames(conserved_markers_dc_pc_go)
c<- rownames(conserved_markers_iphibc_go)
d<- rownames(conserved_markers_tbc_go)
e<- rownames(conserved_markers_idc_go)
f<- rownames(conserved_markers_cc_go)

conserved_marker_go_all <- data.frame(genes = unique(c(deg_macaque_2$names,
                                                       deg_mouse_2$names)))
all.genes <- data.frame(genes = unique(c(a, b,c,d,e,f)))

tf_gene_list <- read.csv("human_tf_gene_list.csv",row.names = 1)

all.genes$HGNC.symbol <- all.genes$genes
conserved_tf <- merge(x=all.genes,y=tf_gene_list, by="HGNC.symbol", all=FALSE)#for plot TF types
write.csv(conserved_tf,file = "conserved_TFs_across_macaque_and_mouse_v3.csv")
# we use these conserved tfs for bar plot in graphpad


write.csv(all.genes,file = "conserved_marker_genes_across_macaque_and_mouse_v3.csv")


