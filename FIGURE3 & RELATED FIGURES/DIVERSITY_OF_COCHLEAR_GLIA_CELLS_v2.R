rm(list=ls())
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(zoo)
library(ggplot2)
library(loomR)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)
library(eulerr)
library(viridis)
library(pvclust)
library(parallel)
library(dendextend)
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize= 62914560000) # 60GB (60000*1024^2)
library(SingleCellExperiment)
library(Matrix)
# REMOVE OLIGODENDROCYTES THAT IS DURING DISSECTION CONTAINING A PORITON OF CENTRAL
#NERVE AND ITS ASSOCIATED OLIGODENDROCYTES, WHICH EXPRESS PLP1, SPECIFIC MOBP+
# substract glia cells from mouse
seurat_mouse_new <- readRDS("F:/PROJECTS/PROJECT_MONKEY SC-RNA SEQ/WORKPLACE/R/seurat_mouse_new.rds")
Idents(seurat_mouse_new) <- "subclasses_label"
table(seurat_mouse_new$subclasses_label)
gc_mouse <- subset(seurat_mouse_new, idents="OL_SC")
gc_mouse
# re-cluster glia cells 
DefaultAssay(gc_mouse) <- "RNA"
gc_mouse <- SCTransform(gc_mouse, vars.to.regress = c("nCount_RNA","pMT"),
                        conserve.memory = T)
gc_mouse <- RunPCA(gc_mouse,npcs = 50)
DimPlot(gc_mouse, reduction = "pca")
gc_mouse <- RunUMAP(gc_mouse, dims = 1:10)
DimPlot(gc_mouse, reduction = "umap")
gc_mouse <- FindNeighbors(gc_mouse,dims = 1:10)
gc_mouse <- FindClusters(gc_mouse,resolution = 1.2)
DimPlot(gc_mouse, label = TRUE,split.by = "orig.ident")
DefaultAssay(gc_mouse) <- "SCT"
markers.to.plot<-c("MBP","PLP1","SOX10","MOBP",
                   "SCN7A",      #NMSC
                   "PMP22","PRX","NCMAP", # MSC
                   "IGFBP6","GLUL","KCNJ10","SLC1A3","SLC6A15","GATA2" # SGC
                   
                   
)
VlnPlot(gc_mouse, features = markers.to.plot)
### RENAME GLIA CELL TYPE ACCORDING TO THESE WELL-CHARACTERIZED MARKERS
#rename celltype
Idents(gc_mouse) <- "seurat_clusters"
table(Idents(gc_mouse))
new.cluster.ids <- c("SGC","MSC","SGC","NMSC","SGC")
names(new.cluster.ids) <- levels(gc_mouse)
gc_mouse<- RenameIdents(gc_mouse, new.cluster.ids)
DimPlot(gc_mouse, reduction = "umap",label = TRUE,repel = TRUE)
gc_mouse$celltype <- Idents(gc_mouse)
DefaultAssay(gc_mouse) <- "SCT"
markers.to.plot<-c("MBP","PLP1","PMP22","MPZ","GFAP"
                   
                   
)
VlnPlot(gc_mouse, features = markers.to.plot)
#rename celltype
Idents(gc_mouse) <- "celltype"
table(Idents(gc_mouse))
new.cluster.ids <- c("SGC_mouse","MSC_mouse","NMSC_mouse")
names(new.cluster.ids) <- levels(gc_mouse)
gc_mouse<- RenameIdents(gc_mouse, new.cluster.ids)
DimPlot(gc_mouse, reduction = "umap",label = TRUE,repel = TRUE)
gc_mouse$celltype_species <- Idents(gc_mouse)
saveRDS(gc_mouse,file = "glia_mouse.rds")


#transfer seurat object into annadata-readable format h5ad

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
sce=gc_mouse
DefaultAssay(sce) <- "RNA"
sceasy::convertFormat(sce, from="seurat", to="anndata",
                      outFile='Mm_gc_python.h5ad')

#######################################################################################
#########################################################################################
#substract glia cells from macaque cochlea
#transfer H5AD into SEURAT OBJECT

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
adata_loom <- connect(filename = "macaque_gc.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('macaque_gc_obs.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('macaque_gc_var.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene
gc_macaque= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                                           
                                           min.cells = 0, 
                                           min.features = 0)

gc_macaque@assays[["RNA"]]@meta.features <- meta_feature
table(gc_macaque$age, gc_macaque$cell_type)
# re-cluster glia cells by joint CCA-embedding
DefaultAssay(gc_macaque) <- "RNA"
gc_macaque$species <- "macaque"
# RENAME CLUSTERS WITH "celltype'
Idents(gc_macaque) <- "age"
table(gc_macaque$age)
new.cluster.ids <- c("Macaque_1Y", "Macaque_5Y","Macaque_11Y")
names(new.cluster.ids) <- levels(gc_macaque)
gc_macaque<- RenameIdents(gc_macaque, new.cluster.ids)
gc_macaque$age <- Idents(gc_macaque)
table(Idents(gc_macaque))
gc_macaque <- SCTransform(gc_macaque, verbose = T, vars.to.regress = c("nCount_RNA"), 
                   conserve.memory = T)
gc_macaque <- RunPCA(gc_macaque,npcs = 50)
DimPlot(gc_macaque,reduction = "pca")
ElbowPlot(gc_macaque , ndims = 50)
gc_macaque <- FindNeighbors(gc_macaque, reduction = "pca", dims = 1:30)
gc_macaque <- FindClusters(gc_macaque)
gc_macaque <- RunUMAP(gc_macaque, reduction = "pca", dims = 1:10)
# Visualization
DimPlot(gc_macaque, reduction = "umap",label = TRUE,group.by = 'age')

combined.list <- SplitObject(gc_macaque, split.by = "age")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], 
                                    vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}
library(scrattch.hicat)
table(gc_macaque$cluster3)
table(gc_macaque$age)
Var.genes.macaque_1y <- select_markers(combined.list$Macaque_11Y@assays$SCT@counts, 
                                            combined.list$Macaque_11Y$cluster3, 
                                    n.markers = 100)
Var.genes.macaque_1y.markers <- Var.genes.macaque_1y$markers

Var.genes.macaque_11y <- select_markers(combined.list$Macaque_11Y@assays$SCT@counts, 
                                       combined.list$Macaque_11Y$seurat_clusters, 
                                       n.markers = 100)
Var.genes.macaque_11y.markers <- Var.genes.macaque_11y$markers
Var.genes.macaque_5y <- select_markers(combined.list$Macaque_5Y@assays$SCT@counts, 
                                       combined.list$Macaque_5Y$cluster3, 
                                       n.markers = 100)
Var.genes.macaque_5y.markers <- Var.genes.macaque_5y$markers

total.Var.genes <- unique(c(Var.genes.macaque_1y$markers,
                            Var.genes.macaque_5y$markers,
                            Var.genes.macaque_11y$markers))

total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Macaque_1Y@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Macaque_11Y@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Macaque_5Y@assays$SCT@counts))]
write.csv(total.Var.genes,file = "features used for cca anchors glia cells macaque V2.csv")

features=total.Var.genes
###consider low common genes across 3 samples calculated by scrattch.hicat, we switch to 
#standar pipeline in seurat
features <- SelectIntegrationFeatures(object.list =combined.list, nfeatures = 2000)
combined.list <- PrepSCTIntegration(object.list =combined.list,
                                    anchor.features = features)
gc_macaque.anchors <- FindIntegrationAnchors(object.list = combined.list,
                                          normalization.method = "SCT", 
                                          anchor.features = features,dims = 1:40,
                                          reduction="cca")
gc_macaque.combined.sct <- IntegrateData(anchorset = gc_macaque.anchors, 
                                      normalization.method = "SCT",dims = 1:40)
gc_macaque.combined.sct <- RunPCA(gc_macaque.combined.sct, features = features, npcs=100)
ElbowPlot(gc_macaque.combined.sct , ndims = 100)
gc_macaque.combined.sct <- FindNeighbors(gc_macaque.combined.sct, 
                                      reduction = "pca", dims = 1:30,nn.eps = 0)
gc_macaque.combined.sct <- FindClusters(gc_macaque.combined.sct,resolution = 0.4)
saveRDS(gc_macaque.combined.sct,file="gc_macaque_anchors_integrated_V2.rds")
gc_macaque_anchors_integrated_V2 <- readRDS("F:/PROJECTS/PROJECT_MONKEY SC-RNA SEQ/WORKPLACE/R/gc_macaque_anchors_integrated_V2.rds")
gc_macaque.combined.sct=gc_macaque_anchors_integrated_V2
#gc_macaque.combined.sct=gc_macaque_anchors_integrated
gc_macaque.combined.sct <- FindNeighbors(gc_macaque.combined.sct, reduction = "pca",
                                         dims = 1:30)
gc_macaque.combined.sct <- FindClusters(gc_macaque.combined.sct,resolution = 0.3)
gc_macaque.combined.sct <- RunUMAP(gc_macaque.combined.sct, reduction = "pca", dims = 1:20)

gc_macaque.combined.sct$seurat_clusters.new <- as.integer(gc_macaque.combined.sct$seurat_clusters)
# Visualization
p1 <- DimPlot(gc_macaque.combined.sct, reduction = "umap", group.by = "age")
p2 <- DimPlot(gc_macaque.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(gc_macaque.combined.sct, reduction = "umap", label = TRUE,group.by = "seurat_clusters")
p1+p2+p3
table(gc_macaque.combined.sct$orig.ident,gc_macaque.combined.sct$seurat_clusters.new)

#########################################################################################
DefaultAssay(gc_macaque.combined.sct) <- "SCT"
markers.to.plot<-c("MBP","PLP1","SOX10","MOBP","GFAP",
                   "SCN7A" ,   #NMSC
                   "PMP22","PRX","NCMAP",    # MSC
                   "IGFBP6","GLUL","KCNJ10","SLC1A3","GATA2" # SGC "SLC6A15",
                   
                   
)

VlnPlot(gc_macaque.combined.sct, features = markers.to.plot)
VlnPlot(gc_macaque.combined.sct, features = "GFAP")

Idents(gc_macaque.combined.sct) <- "seurat_clusters"
table(Idents(gc_macaque.combined.sct))
new.cluster.ids <- c("GC_1","GC_2","GC_3")
names(new.cluster.ids) <- levels(gc_macaque.combined.sct)
gc_macaque.combined.sct<- RenameIdents(gc_macaque.combined.sct, new.cluster.ids)
DimPlot(gc_macaque.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
gc_macaque.combined.sct$celltype <- Idents(gc_macaque.combined.sct)
DimPlot(gc_macaque.combined.sct, reduction = "umap", label = TRUE,
        group.by = "celltype",pt.size = 2)



# Visualization
p1 <- DimPlot(gc_macaque.combined.sct, reduction = "umap", 
              group.by = "age",pt.size = 1.5)
p2 <- DimPlot(gc_macaque.combined.sct, reduction = "umap", label = TRUE,
              group.by = "celltype",pt.size = 1.5)
p2+p1
# CHANGE CLUSTER_LABEL COLORS
# umap cluster colors
table(gc_macaque.combined.sct$celltype)
levels(Idents(gc_macaque.combined.sct))


my_cols <- c(
  "GC_1"="#bcaf26","GC_2"="#c09ecd","GC_3"="#8faae8"
  
)


my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
p2=DimPlot(gc_macaque.combined.sct,
           cols = my_cols2, label=TRUE , repel=TRUE,reduction = "umap",pt.size = 1.5
           ,group.by = "celltype")
p2
p2+p1
p2
###
FeaturePlot(gc_macaque.combined.sct,features = c("PLP1","SOX10","SCN7A","NCMAP","IGFBP6"),
            label = TRUE,pt.size = 1)
FeaturePlot(gc_macaque.combined.sct,features = "PLP1",
            label = TRUE,pt.size = 1)
library(dittoSeq)
dittoDimPlot(gc_macaque.combined.sct,"MBP",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_macaque.combined.sct,"MPZ",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_macaque.combined.sct,"SCN7A",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_macaque.combined.sct,"NCMAP",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_macaque.combined.sct,"IGFBP6",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_macaque.combined.sct,"GFAP",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
p1=multi_dittoDimPlot(gc_macaque.combined.sct,c("MBP","SCN7A","NCMAP","IGFBP6"),size = 2,
                   min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
                   nrow = 1,ncol = 4
                   )
p1
DefaultAssay(gc_mouse) <- "SCT"
dittoDimPlot(gc_mouse,"MBP",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_mouse,"MPZ",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_mouse,"SCN7A",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_mouse,"NCMAP",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_mouse,"IGFBP6",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
dittoDimPlot(gc_mouse,"GFAP",size = 2,
             min.color = "#E0FFFF",max.color = "#C71585",order = "increasing")
p2=multi_dittoDimPlot(gc_mouse,c("MBP","SCN7A","NCMAP","IGFBP6"),size = 2,
                   min.color = "#E0FFFF",max.color = "#C71585",order = "increasing",
                   nrow = 1,ncol = 4
)
p2
grid.arrange(p1,p2)
saveRDS(gc_macaque.combined.sct,file = "glia_macaque.rds")
#rename celltype
Idents(gc_macaque.combined.sct) <- "celltype"
table(Idents(gc_macaque.combined.sct))
new.cluster.ids <- c("GC_1_macaque","GC_2_macaque","GC_3_macaque")
names(new.cluster.ids) <- levels(gc_macaque.combined.sct)
gc_macaque.combined.sct<- RenameIdents(gc_macaque.combined.sct, new.cluster.ids)
DimPlot(gc_macaque.combined.sct, reduction = "umap",label = TRUE,repel = TRUE)
gc_macaque.combined.sct$celltype_species <- Idents(gc_macaque.combined.sct)
saveRDS(gc_macaque.combined.sct,file = "glia_macaque.rds")


#transfer seurat object into annadata-readable format h5ad

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
sce=gc_macaque.combined.sct
DefaultAssay(sce) <- "RNA"
sceasy::convertFormat(sce, from="seurat", to="anndata",
                      outFile='Mac_gc_python.h5ad')
#########################################################################################
##############cca-based label transfer for projecting macaque glia to mouse glia 
table(gc_macaque.combined.sct$celltype)
table(gc_mouse$celltype)
library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)
# ============ Load data ============ 
seurat_macaque_gc=gc_macaque.combined.sct # query
seurat_mouse_gc=gc_mouse # reference
table(seurat_macaque_gc$celltype)

table(seurat_macaque_gc$celltype)
table(seurat_mouse_gc$celltype)
# Names of macaque genes 
matrix_macaque=seurat_macaque_gc@assays$RNA@counts


meta_macaque=seurat_macaque_gc@meta.data
# create a seurat object from macaque dataset that is compatible with mouse dataset
seurat_macaque_gc <- CreateSeuratObject(counts = matrix_macaque)
seurat_macaque_gc
#preprocessing of macaque data
seurat_macaque_gc <- AddMetaData(seurat_macaque_gc, metadata = meta_macaque)
seurat_macaque_gc <- SCTransform(object = seurat_macaque_gc,
                                  vars.to.regress = c("nCount_RNA"),
                                  conserve.memory = T)
seurat_macaque_gc <- RunPCA(seurat_macaque_gc)
seurat_mouse_gc <- SCTransform(object = seurat_mouse_gc,
                                vars.to.regress = c("nCount_RNA"),
                                conserve.memory = T)
seurat_mouse_gc <- RunPCA(seurat_mouse_gc)
# ============ projection of reference data onto query object ============ 
# We use all default parameters here for identifying anchors
transfer.anchors <- FindTransferAnchors(reference = seurat_mouse_gc,
                                        query = seurat_macaque_gc, 
                                        normalization.method = "SCT",
                                        reference.assay = "SCT", query.assay = "SCT", 
                                        reduction = "cca",dims = 1:30)


# TransferData returns a matrix with predicted IDs and prediction scores, which are written into meta.data
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = seurat_mouse_gc$celltype, 
                                     weight.reduction = seurat_macaque_gc[["pca"]],
                                     dims = 1:30)


seurat_macaque_gc <- AddMetaData(seurat_macaque_gc, metadata = celltype.predictions)

# Cells with low anchoring scores (<=0.25) are either low-quality cells or cells with ambigous identity.
# They are marked as low-confidence cells and removed later.
cutoff=0.25
seurat_macaque_gc$hiConfi=ifelse(seurat_macaque_gc$prediction.score.max > cutoff,'TRUE','FALSE')



# ============ Intergration of two datasets for visualization ============ 

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_mouse_gc)
refdata <- GetAssayData(seurat_mouse_gc, assay = "SCT", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, 
                           refdata = refdata, 
                           weight.reduction = seurat_macaque_gc[["pca"]],dims = 1:30)

# this line adds the imputed data matrix to the seurat_MACAQUE object
seurat_macaque_gc[["SCT"]] <- imputation
seurat_macaque_gc$species = "Macaque"
seurat_mouse_gc$species= "Mouse"
coembed <- merge(x = seurat_mouse_gc, y = seurat_macaque_gc)

# Finally, we run PCA and UMAP on this combined object, 
#to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)

coembed <- RunUMAP(coembed, dims = 1:10
)
p1 <- DimPlot(coembed,reduction = "umap",label = TRUE,repel = TRUE, 
              group.by = "species",pt.size = 1.5)
p2 <- DimPlot(coembed,reduction = "umap",label = TRUE,repel = TRUE,
              group.by = "celltype",
              pt.size = 1.5)
p1+p2

# Finally, we run PCA and UMAP on this combined object, 
#to visualize the co-embedding of both
# datasets


coembed = AddMetaData(coembed,Embeddings(coembed[["umap"]]),
                      colnames(Embeddings(coembed[["umap"]])))

coembed$species=ifelse(!is.na(coembed$predicted.id), 'Macaque', 'Mouse')
# summary plots
d1 <- DimPlot(coembed, group.by = "species")
d2 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')
              %>%filter(species=='Mouse')%>%pull(V1),group.by = "celltype",
              label = TRUE, repel = TRUE)+labs(title='Mouse celltype')

d3 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')
              %>%filter(species=='Macaque')%>%pull(V1),group.by = "celltype",
              label = TRUE, repel = TRUE)+labs(title='Macaque celltype')
d4 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')
              %>%filter(species=='Macaque')%>%pull(V1),
              group.by = "predicted.id", label = TRUE, repel = TRUE)+labs(title='Macaque predicted.id')
library(viridis)
d5 = FeaturePlot(coembed, 'prediction.score.max',reduction = "umap",cols =viridis(100))
d6 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'))


grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)
DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%
          pull(V1),group.by = "celltype",
        label = TRUE, repel = TRUE)+labs(title='celltype')

saveRDS(coembed,file='Seurat_macaque_mouse_gc_LABEL_TRNSFER_V2.rds')

table(seurat_macaque_gc$celltype)
table(seurat_macaque_gc$celltype,seurat_macaque_gc$predicted.id)

dev.off()


results <- table(seurat_macaque_gc@meta.data$celltype, 
                 seurat_macaque_gc@meta.data$predicted.id)

######################################################################################
############integration of macaque and mouse glia cell by using the Seurat cca pipeline
#with different numbers of high variable genes (hvgs)
seurat_macaque_gc=gc_macaque.combined.sct 
seurat_mouse_gc=gc_mouse 
table(seurat_macaque_gc$celltype)
table(seurat_mouse_gc$celltype)

DefaultAssay(seurat_macaque_gc) <- "RNA"
DefaultAssay(seurat_mouse_gc) <- "RNA"
seurat_macaque_gc$species <- "macaque"
seurat_mouse_gc$species <- "mouse"
table(seurat_macaque_gc$age,seurat_macaque_gc$species)
table(seurat_mouse_gc$age,seurat_mouse_gc$species)
seurat_mouse_gc$age=seurat_mouse_gc$orig.ident

# normalized with SCTransform v2

all.data <- merge(x = seurat_macaque_gc, y = seurat_mouse_gc)

combined.list <- SplitObject(all.data, split.by = "species")

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi",vst.flavor="v2")
}
rm(all.data)
ngenes=3000 # shift to 500, 1000, 3000
comb.features <- SelectIntegrationFeatures(object.list = combined.list, 
                                           nfeatures = ngenes)
comb.list.prep <- PrepSCTIntegration(object.list = combined.list, 
                                     anchor.features = comb.features, 
                                     verbose = FALSE)

comb.anchors <- FindIntegrationAnchors(object.list = comb.list.prep, 
                                       normalization.method = "SCT", 
                                       anchor.features = comb.features,
                                       reduction="cca",
                                       k.filter = 120)
comb.integrated <- IntegrateData(anchorset = comb.anchors, 
                                 normalization.method = "SCT", 
                                  dims=1:30)
comb.integrated <- RunPCA(comb.integrated, npcs=30)
ElbowPlot(comb.integrated , ndims = 100)
comb.integrated <- RunUMAP(comb.integrated, dims = 1:30)
comb.integrated <- FindNeighbors(comb.integrated, dims = 1:20)
comb.integrated <- FindClusters(comb.integrated)
comb.integrated$seurat_clusters.new <- as.integer(comb.integrated$seurat_clusters)

# renaming clusters starting from 1 instead of 0
DefaultAssay(comb.integrated) <- "integrated"
new.ids <- seq(1, length(levels(Idents(comb.integrated))), by=1)
names(new.ids) <- levels(Idents(comb.integrated))
comb.integrated <- RenameIdents(comb.integrated, new.ids)
species_annot <- comb.integrated@meta.data$species
meta_df <- comb.integrated@meta.data
names(species_annot) <- rownames(comb.integrated@meta.data)
meta_df$speciesTree <- species_annot

comb.integrated <- AddMetaData(comb.integrated, meta_df)

# Visualization

p1 <- DimPlot(comb.integrated, reduction = "umap", group.by = "species")
p2 <- DimPlot(comb.integrated, reduction = "umap", label = TRUE)
p3 <- DimPlot(comb.integrated, reduction = "umap", label = TRUE,group.by = "celltype")
p1+p2+p3
table(comb.integrated$seurat_clusters,comb.integrated$celltype)

comb.integrated_500 = comb.integrated
comb.integrated_1000 = comb.integrated
comb.integrated_2000 = comb.integrated
comb.integrated_3000 = comb.integrated

comb.integrated_500$seurat_clusters.new <- comb.integrated_2000$seurat_clusters.new
comb.integrated_1000$seurat_clusters.new <- comb.integrated_2000$seurat_clusters.new
comb.integrated_3000$seurat_clusters.new <- comb.integrated_2000$seurat_clusters.new

table(comb.integrated_500$seurat_clusters, comb.integrated_500$seurat_clusters.new)
table(comb.integrated_1000$seurat_clusters, comb.integrated_1000$seurat_clusters.new)
# species integration
p1 <- DimPlot(comb.integrated_500, reduction = "umap", group.by = "species",
              cols = c("#F08080","#4169E1"),pt.size = 1.5)
p1
p2 <- DimPlot(comb.integrated_1000, reduction = "umap", group.by = "species",
              cols = c("#F08080","#4169E1"),pt.size = 1.5)
p2
p3 <- DimPlot(comb.integrated_2000, reduction = "umap", group.by = "species",
              cols = c("#F08080","#4169E1"),pt.size = 1.5)
p3
p4 <- DimPlot(comb.integrated_3000, reduction = "umap", group.by = "species",
              cols = c("#F08080","#4169E1"),pt.size = 1.5)
p4


grid.arrange(p1,p2,p3,p4, ncol = 4, nrow = 1)

# integrated_clusetrs from 2000 hvgs
p5 <- DimPlot(comb.integrated_500, reduction = "umap", group.by = "seurat_clusters.new",
              pt.size = 1.5,label = TRUE,repel = TRUE)
p6 <- DimPlot(comb.integrated_1000, reduction = "umap", group.by = "seurat_clusters.new",
              pt.size = 1.5,label = TRUE,repel = TRUE)
p7 <- DimPlot(comb.integrated_2000, reduction = "umap", group.by = "seurat_clusters.new",
              pt.size = 1.5,label = TRUE,repel = TRUE)
p8 <- DimPlot(comb.integrated_3000, reduction = "umap", group.by = "seurat_clusters.new",
              pt.size = 1.5,label = TRUE,repel = TRUE)


grid.arrange(p5,p6,p7,p8, ncol = 4, nrow = 1)


p5 <- DimPlot(comb.integrated_500, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 1.5,label = TRUE,repel = TRUE)
p6 <- DimPlot(comb.integrated_1000, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 1.5,label = TRUE,repel = TRUE)
p7 <- DimPlot(comb.integrated_2000, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 1.5,label = TRUE,repel = TRUE)
p8 <- DimPlot(comb.integrated_3000, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 1.5,label = TRUE,repel = TRUE)


grid.arrange(p5,p6,p7,p8, ncol = 4, nrow = 1)


# celltype



my_cols <- c("GC_1"="#87CEFA","GC_2"="#6495ED","GC_3"="#6A5ACD",
             "NMSC"="#DC143C",
             "MSC"="#B22222","SGC"="#FFA07A"
             
)



my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
p9 <- DimPlot(comb.integrated_500, reduction = "umap", group.by = "celltype",
              pt.size = 1.5, cols = my_cols2,label = TRUE)
p10 <- DimPlot(comb.integrated_1000, reduction = "umap", group.by = "celltype",
              pt.size = 1.5, cols = my_cols2,label = TRUE)
p11 <- DimPlot(comb.integrated_2000, reduction = "umap", group.by = "celltype",
              pt.size = 1.5, cols = my_cols2,label = TRUE)
p12 <- DimPlot(comb.integrated_3000, reduction = "umap", group.by = "celltype",
              pt.size = 1.5, cols = my_cols2,label = TRUE)


grid.arrange(p9,p10,p11,p12, ncol = 4, nrow = 1)

####% CELL PERCENTAGE ACROSS INTEGRTAED CLUSTERS

results_500 <- table(comb.integrated_500@meta.data$celltype, 
                 comb.integrated_500@meta.data$seurat_clusters)

results_1000 <- table(comb.integrated_1000@meta.data$celltype, 
                     comb.integrated_1000@meta.data$seurat_clusters)

results_2000 <- table(comb.integrated_2000@meta.data$celltype, 
                      comb.integrated_2000@meta.data$seurat_clusters)

results_3000 <- table(comb.integrated_3000@meta.data$celltype, 
                      comb.integrated_3000@meta.data$seurat_clusters)


library(zoo)
library(loomR)
library(scater)
library(stringr)
library(corrplot)
library(matrixStats)

results.norm_500 <- 100*(results_500/rowSums(results_500))
mus.ids_500 <- colnames(results.norm_500)
results.norm_500 <- results.norm_500[,mus.ids_500]

corrplot(results.norm_500, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)


mus_names <- c("0","1","2","3","4","5")
mac_names_anatomical <- c("MSC","NMSC","SGC","GC_1","GC_2","GC_3")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered_500 <- results.norm_500[order(match(rownames(results.norm_500),
                                                   mac_names_anatomical_rev)),
                                               order(match(colnames(results.norm_500), 
                                                           mus_names))]
corrplot(results.norm_reordered_500, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)

#1000
results.norm_1000 <- 100*(results_1000/rowSums(results_1000))
mus.ids_1000 <- colnames(results.norm_1000)
results.norm_1000 <- results.norm_1000[,mus.ids_1000]

corrplot(results.norm_1000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)


mus_names <- c("0","1","2","3","4","5")
mac_names_anatomical <- c("MSC","NMSC","SGC","GC_1","GC_2","GC_3")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered_1000 <- results.norm_1000[order(match(rownames(results.norm_1000),
                                                           mac_names_anatomical_rev)),
                                               order(match(colnames(results.norm_1000), 
                                                           mus_names))]
corrplot(results.norm_reordered_1000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)

#2000
results.norm_2000 <- 100*(results_2000/rowSums(results_2000))
mus.ids_2000 <- colnames(results.norm_2000)
results.norm_2000 <- results.norm_2000[,mus.ids_2000]

corrplot(results.norm_2000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)


mus_names <- c( "0",  "1","2","3","4","5","6","7")
mac_names_anatomical <- c("MSC","NMSC","SGC","GC_1","GC_2","GC_3")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered_2000 <- results.norm_2000[order(match(rownames(results.norm_2000),
                                                             mac_names_anatomical_rev)),
                                                 order(match(colnames(results.norm_2000), 
                                                             mus_names))]
corrplot(results.norm_reordered_2000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)


#3000

results.norm_3000 <- 100*(results_3000/rowSums(results_3000))
mus.ids_3000 <- colnames(results.norm_3000)
results.norm_3000 <- results.norm_3000[,mus.ids_3000]

corrplot(results.norm_3000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)


mus_names <- c("0",   "1","2","3","4","5","6","7")
mac_names_anatomical <- c("MSC","NMSC","SGC","GC_1","GC_2","GC_3")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered_3000 <- results.norm_3000[order(match(rownames(results.norm_3000),
                                                             mac_names_anatomical_rev)),
                                                 order(match(colnames(results.norm_3000), 
                                                             mus_names))]
corrplot(results.norm_reordered_3000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)


# % CELL PERCENTAGE ACROSS INTEGRTAED CLUSTERS
p1 = corrplot(results.norm_reordered_500, order="original",tl.pos="lt",
               method="color", tl.col="black", is.corr=F)
p2 = corrplot(results.norm_reordered_1000, order="original",tl.pos="lt",
               method="color", tl.col="black", is.corr=F)
p3 = corrplot(results.norm_reordered_2000, order="original",tl.pos="lt",
               method="color", tl.col="black", is.corr=F)
p4 = corrplot(results.norm_reordered_3000, order="original",tl.pos="lt",
               method="color", tl.col="black", is.corr=F)

### average gene expression profiles were computed for the clusters identified 
#in the integrated space (2000 hvgs), and then used for hierarchical clustering.
library(pvclust)
library(parallel)
library(dendextend)
library(future)
library(Seurat)
library(dplyr)
library(cowplot)
library(zoo)
library(ggplot2)
library(scater)
library(stringr)
library(gplots)
library(matrixStats)
sce=comb.integrated_2000
table(sce$seurat_clusters,sce$seurat_clusters.new)
Idents(sce) <- "seurat_clusters"
DefaultAssay(sce)<-"SCT"
neurons_avg <- AverageExpression(sce)
neurons_avg_data <- neurons_avg[["SCT"]]
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
noisy.liz <- NoisyGenes(sce, 0.2, levels(Idents(sce)))
neurons_avg_data_filtered <- neurons_avg_data[!rownames(neurons_avg_data) %in% noisy.liz,]
dist_func <- function(x){x <- as.dist(1-cor(x, method="s"))}
pvclust.TFs_dist <- pvclust(neurons_avg_data_filtered, method.dist=dist_func, 
                            method.hclust="ward.D2", nboot=10000, parallel=T)
plot(pvclust.TFs_dist)
pvrect(pvclust.TFs_dist, alpha=0.90)
#plotting dendrogram colore by significance
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>% pvclust_show_signif_gradient(pvclust.TFs_dist, signif_col_fun = colorRampPalette(c("lightgrey","grey","black")), signif_type="au") %>%
  plot(main = "GC clusters branch confidence")
dend <- rotate(dend)
dendrogram_names <- labels(dend)
dend <- as.dendrogram(pvclust.TFs_dist)
dend %>%
  pvclust_show_signif(pvclust.TFs_dist, signif_value = c("black", "black"), show_type = "col") %>%
  plot(main = "Cluster dendrogram with AU/BP values (%)")

pvrect2(pvclust.TFs_dist, alpha=0.90)
plot(dend)

# building the Dendrogram with leaves colores by species' mixing
# get the data from integrated assay
library(Seurat)
library(corrplot)
library(MatrixGenerics)
library(ggplot2)
library(dplyr)
library(gplots)
library(Hmisc)
library(stringr)
library(speciesTree)
# change comb.integrated number 500 1000 2000 3000
table(comb.integrated_3000$seurat_clusters)
Idents(comb.integrated_3000) <- "seurat_clusters"
expression.matrix <- GetAssayData(comb.integrated_3000, 
                                  slot="data", 
                                  assay="integrated")
meta_clusters <- as.integer(Idents(comb.integrated_3000))
names(meta_clusters) <- names(Idents(comb.integrated_3000))
upperlevelinfo = NULL
species_annot <- comb.integrated_2000@meta.data$species
names(species_annot) <- rownames(comb.integrated_3000@meta.data)

# obtaining distance matrix and performing hierarchical clustering
d <- cluster.matrix.expression.distances(t(expression.matrix),
                                         groups=meta_clusters, dist="cor", 
                                           useVariablegenes=FALSE, 
                                         use.scaled.data=TRUE)


dendr <- hclust(as.dist(d), method='ward.D2')
dend <- as.dendrogram(dendr)

# coloring dendrogram leaves according to the proportion of cells from each species in each cluster
dendr <- TransferDend(dend, renameCluster = FALSE, cls.groups = meta_clusters)
cls.groups <- dendr$new.groups
dend <- dendr$dendrogram
leafcontent <- dendr$leafcontent
stability.measurements = NULL
dend <- AddTreeAttribute(dend, species_annot, leafcontent)
dend <- dendSetWidthBysize(dend, scale = 8)
colorpallete <- colorRampPalette(c("blue", "grey", "grey",  "grey", "red"))(101)

upperLevelnodes = NULL
fac <- as.factor(species_annot)
totalCells <- table(fac)
cc2col <- function(cc, rate=15, base=0.001){
  cc <- round((cc[2:3]/totalCells)/sum(cc[2:3]/totalCells), 2)
  cv <- cc
  cv <- dexp(cv, rate)
  cv <- cv/rate * (1-base)
  col <- adjustcolor(rgb(cv[1],cv[2], 1), offset = c(0.5, 0.5, 0.5, 0.1))
  return(col)
}

cbm <- function(d,fac) {
  if(is.leaf(d)) {
    #lc <- fac[leafContent[[as.numeric(attr(d,'label'))]]]
    cc <- attr(d, "cc")
    col <- cc2col(cc)
    attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
    return(d);
  } else {
    oa <- attributes(d);
    d <- lapply(d,cbm,fac=fac);
    attributes(d) <- oa;
    cc <- attr(d, "cc")
    col <- cc2col(cc)
    attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
    return(d);
  }
}
dend <- cbm(dend,species_annot)
tree <- dend

plot(tree)


mus_names <- c("0","3","1","6","2","7","4","5")
mac_names_anatomical <- c("MSC","NMSC","SGC","GC_3","GC_2","GC_1")
mac_names_anatomical_rev <- rev(mac_names_anatomical)
results.norm_reordered_3000 <- results.norm_3000[order(match(rownames(results.norm_3000),
                                               mac_names_anatomical_rev)),
                                 order(match(colnames(results.norm_3000), 
                                                           mus_names))]
corrplot(results.norm_reordered_3000, order="original",tl.pos="lt",
         method="color", tl.col="black", is.corr=F)
dev.off()
#######################################################################################
######========================META-NEIGHBOR ANALYIS=================================

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
sce=comb.integrated_500 # change 1000 2000 3000 
table(sce$celltype,sce$seurat_clusters)
################## Prep SingleCellExperiment object for analysis  #############################################
library(SingleCellExperiment)
library(Matrix)
sce$class_label = "GC"
table(sce$class_label)
Idents(sce) <- "class_label"
sce$class.species <- paste(sce$class_label, sce$species, sep = "_")
table(sce$class.species)
Idents(sce) <- "species"

mac <- subset(sce, idents ="macaque")
ms <- subset(sce, idents ="mouse")

mac <- GetAssayData(object = mac[["RNA"]], slot = "counts")
ms <- GetAssayData(object = ms[["RNA"]], slot = "counts")

gs=read.csv("orthologTable_human_macaque_mouse.csv",header=T,stringsAsFactors=F)  ### read in ortholog mapping table


mac <- as.matrix(mac)
ms <- as.matrix(ms)

p2 <- subset(sce, idents ="macaque")
p3<- subset(sce, idents ="mouse")

p2$sample_id <- colnames(p2)
p3$sample_id <- colnames(p3)

p2$study_id="macaque"
p3$study_id="mouse"
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

sce_all$sample_id_append <- sce_all$sample_id
head(sce_all$sample_id_append)
sce$sample_id <-colnames(sce)
head(sce$sample_id)
m<-match(sce_all$sample_id_append,sce$sample_id)
sum(!is.na(m))
f.a=!is.na(m)
f.b=m[f.a]

sce_all$celltype[f.a]=as.character(sce$celltype[f.b])
#sce_all$final_integrated_cluster_color[f.a]=cochlea_epithelium$[f.b]
head(sce_all$study_id)
head(sce_all$celltype[f.a])
sce_all=sce_all[,f.a]
head(sce_all$study_id)
saveRDS(sce_all,file="mac_ms_gc_hvgs_500_V2.rds")
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
                                         sce_all$seurat_clusters[f.b],
                                         one_vs_one=TRUE)
}

b<- as.data.frame(mn_1v1_subclass)
write.csv(b,file = "mn_1v1_gc_celltype_500_V2.csv")
options(stringsAsFactors = FALSE)

d<- read.csv("mn_1v1_gc_celltype_500_V2.csv",row.names = 1)
library(pheatmap)
library(viridis)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE)
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = viridis(50))
#color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap::pheatmap(d,cluster_rows = FALSE,cluster_cols = FALSE,border=FALSE,color = cividis(50))

###########################################################################################
DefaultAssay(comb.integrated_2000) <- "SCT"
FeaturePlot(comb.integrated_2000, features = "GFAP")
VlnPlot(comb.integrated_2000,features = "S100B",group.by = "celltype")

gc=comb.integrated_2000
gc <- SCTransform(gc,vars.to.regress = "nCount_RNA")
genes_to_plot <- c("GLUL","SLC1A2","SLC7A11","SLC1A4","SLC7A10","SLC6A1","SLC6A11",
                   "SLC6A9","SLC2A3","SLC44A3","SLC44A5",
                   "MERTK","MEGF10","ADGRB1","LGALS3","PROS1","LRP1","C1QL1","C1QL4",
                   "EZR","RDX","TNR","TNC","CNTNAP2","EPB41L3","CNTN1",
                   "DPP10","ATP1B1","ATP1A1","KCNQ5","KCNK2","KCNN2","KCNMA1","ATP1A2","ATP1B2",
                   "FXYD1","KCNIP1","KCNJ16","KCNJ10","KCNJ6","KCND2","KCND3",
                   "GRINA","GRIA2","GRIA4","GRIK2","GRIK1","GRID1","GRM3","GRM5","GRM7","GRM8",
                   "GABBR1","GABBR2","GABRB1","GABRG3",
                   "FGF1","BMP7","WNT4","WNT7B",
                   "VEGFA","TJP2","CLDN10","GAN","FGD4","GPC1","GPC3","GPC5","SPARC",
                   "SPARCL1")
Idents(gc) <- "celltype"
Idents(gc) <- factor(Idents(gc), 
                               levels = c("GC_1","GC_2","GC_3","MSC","NMSC","SGC"))
table(Idents(gc))
gc$celltype2 <- Idents(gc)
DotPlot(gc, features = genes_to_plot, dot.scale = 8,group.by = "celltype2",
        cols  =c("white", "#ad9300")) + RotatedAxis()

genes_to_plot <- c("GLUL",
                   "MERTK","MEGF10","PROS1","LRP1","CNTN1",
                   "DPP10","ATP1B1","ATP1A1","KCNQ5","KCNMA1","ATP1A2",
                  "KCNJ10",
                 "GRIA2","GRIK2","GABRG3",
                   "FGF1","BMP7",
                   "TJP2","GAN","FGD4","GPC5","SPARC",
                   "SPARCL1")

DotPlot(gc, features = genes_to_plot, dot.scale = 8,group.by = "celltype2",
        cols  =c("white", "#ad9300")) + RotatedAxis()


genes_to_plot <- c("LRP1","CLU","PROS1","AXL","MEGF10","MERTK","ITGAV","GULP1", # GENES RELATED TO CLEARANCE
                   "GLUL","SLC7A11","SLC44A3","SLC44A5","CRYAB", #NT TRANSPORT AND SYNTHESIS
                   "CNTN1","RDX","EPB41L3", #PARANODAL FUNCTION
                   "DPP10","ATP1B1","ATP1A1","ATP1A2","KCNQ5","KCNMA1","KCNIP4", "KCNJ10", #POTASSIUM BUFFERING
                   "GRIA2", "GRIK2", "GABRG3", #NEUROTRANSMITTER RECEPTOR
                   "FGF1","BMP7",#MORPHOGEN
                   "GPC5","SPARC","SPARCL1" # SYNAPSE FORMATION
)

DotPlot(gc, features = genes_to_plot, dot.scale = 10,group.by = "celltype2",
        cols  =c("white", "#2141B5")) + RotatedAxis()

saveRDS(comb.integrated_500,file = "gc_hvgs_500.rds")
saveRDS(comb.integrated_1000,file = "gc_hvgs_1000.rds")
saveRDS(comb.integrated_2000,file = "gc_hvgs_2000.rds")
saveRDS(comb.integrated_3000,file = "gc_hvgs_3000.rds")

# alternative integration methods
####################################
########### Harmony ################
####################################
library(harmony)
sce=comb.integrated_2000
sce <- ScaleData(sce)
sce <- RunHarmony(sce, group.by.vars = "species", assay.use = "SCT", project.dim = F)
# careful because these commands overwrite existing ones
sce <- FindNeighbors(sce, reduction="harmony", dims=1:10)
# sce <- FindClusters(sce, res=0.8)
sce <- RunUMAP(sce, reduction="harmony", dims=1:3)
# save results
harmony.embeddings <- Embeddings(sce, reduction="harmony")
harmony.umap <- Embeddings(sce, reduction="umap")
saveRDS(list(harmony.embeddings, harmony.umap), file="Final_2sp_harmony.rds")
# save seurat object
saveRDS(sce, file="Final_2sp_integration_harmony.rds")

# to re-attach Harmony results to seurat object:
# harmony.data <- readRDS("Final_2sp_harmony.rds")
# sce[["harmony"]] <- CreateDimReducObject(embeddings = harmony.data[[2]], key = "harmony_", assay = "RNA")

# plots 
# species_cols <- c("#EECC66", "#6699CC","#EE99AA","#FFA500")
#species_cols <- c("#EECC66", "#EE99AA", "#6699CC")
species_cols <- c("#EE99AA", "#6699CC")
DimPlot(sce, reduction="harmony", group.by = "species",
        cols = adjustcolor(species_cols, alpha=0.7), shuffle=T, pt.size=1.5)

DimPlot(sce, reduction="harmony", group.by = "species", 
        cols = adjustcolor(species_cols, alpha=0.7), shuffle = T, pt.size=1.5) + NoLegend()
my_cols <- c("SC_1"="#87CEFA","SC_2"="#6495ED","SC_3"="#6A5ACD",
             "NMSC"="#DC143C",
             "MSC"="#B22222","SGC"="#FFA07A"
             
)
Idents(sce) <- "celltype"
DimPlot(sce, reduction="harmony", group.by = "celltype",
         shuffle=T, pt.size=1.5)
DimPlot(sce, reduction="harmony", 
        cols = adjustcolor(my_cols, alpha=0.7), shuffle=T, pt.size=2)


#############################
sce <- merge(x=seurat_macaque_gc, y=seurat_mouse_gc)
table(sce$celltype)

sce <- sce %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
sce <- RunPCA(sce, features = VariableFeatures(sce), npcs = 30)



sce <- RunHarmony(sce, group.by.vars = "species", assay.use = "SCT", project.dim = F)
# careful because these commands overwrite existing ones
sce <- FindNeighbors(sce, reduction="harmony", dims=1:30)
# sce <- FindClusters(sce, res=0.8)
sce <- RunUMAP(sce, reduction="harmony", dims=1:30)
# save results
harmony.embeddings <- Embeddings(sce, reduction="harmony")
harmony.umap <- Embeddings(sce, reduction="umap")
saveRDS(list(harmony.embeddings, harmony.umap), file="Final_2sp_harmony.rds")
# save seurat object
saveRDS(sce, file="Final_2sp_integration_harmony.rds")

# to re-attach Harmony results to seurat object:
# harmony.data <- readRDS("Final_2sp_harmony.rds")
# sce[["harmony"]] <- CreateDimReducObject(embeddings = harmony.data[[2]], key = "harmony_", assay = "RNA")

# plots 
# species_cols <- c("#EECC66", "#6699CC","#EE99AA","#FFA500")
#species_cols <- c("#EECC66", "#EE99AA", "#6699CC")
species_cols <- c("#EE99AA", "#6699CC")
DimPlot(sce, reduction="harmony", group.by = "species",
        cols = adjustcolor(species_cols, alpha=0.7), shuffle=T, pt.size=1.5)

DimPlot(sce, reduction="harmony", group.by = "species", 
        cols = adjustcolor(species_cols, alpha=0.7), shuffle = T, pt.size=1.5) + NoLegend()
my_cols <- c("SC_1"="#87CEFA","SC_2"="#6495ED","SC_3"="#6A5ACD",
             "NMSC"="#DC143C",
             "MSC"="#B22222","SGC"="#FFA07A"
             
)
Idents(sce) <- "celltype"
DimPlot(sce, reduction="harmony", group.by = "celltype",
        shuffle=T, pt.size=1.5)
DimPlot(sce, reduction="harmony", 
        cols = my_cols, shuffle=T, pt.size=1.5)
###############################################################################
sce <- merge(x=seurat_macaque_gc, y=seurat_mouse_gc)
table(sce$celltype)

DefaultAssay(sce) <- "RNA"
sce <- sce %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
sce <- SCTransform(sce, verbose = T, conserve.memory = T,vars.to.regress = 
                         c("nCount_RNA"))
sce <- RunPCA(sce, features = VariableFeatures(sce), npcs = 50)

sce_all_harmony <- RunHarmony(sce,group.by.vars = 'species',reduction = "pca",
                              dims.use = 1:50,assay.use = "SCT")
sce[["harmony"]] <- sce_all_harmony[["harmony"]]
sce <- RunUMAP(sce,dims = 1:50,
                   reduction = "harmony",reduction.name = "umap_harmony")

p1 <- DimPlot(sce, reduction = "umap_harmony", group.by = "species",pt.size = 2) + 
  ggtitle("UMAP Harmony")

p1
p2 <- DimPlot(sce, reduction = "umap_harmony",pt.size = 2,
              group.by = "celltype",repel = TRUE, label = TRUE) + 
  ggtitle("UMAP Harmony")
p2

p1+p2

species_cols <- c("#EE99AA", "#6699CC")
p1=DimPlot(sce, reduction="umap_harmony", group.by = "species",
        cols = species_cols, shuffle=T, pt.size=2)
my_cols <- c("SC_1"="#87CEFA","SC_2"="#6495ED","SC_3"="#6A5ACD",
             "NMSC"="#DC143C",
             "MSC"="#B22222","SGC"="#FFA07A"
             
)
Idents(sce) <- "celltype"
p2=DimPlot(sce, reduction="umap_harmony", cols = my_cols,
        group.by = "celltype",
        shuffle=T, pt.size=2)

p1+p2

saveRDS(sce, file="Final_2sp_gc_integration_harmony.rds")
######################################################################################
####################################
########### SCVI ################
####################################
## Preparing next object for scvi integration
sce <- merge(x=seurat_macaque_gc, y=seurat_mouse_gc)
table(sce$celltype)
combined.list <- SplitObject(sce, split.by = "species")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi",vst.flavor="v2")
}
macaque.vf <- VariableFeatures(combined.list[[1]])[1:3000] 
mouse.vf <- VariableFeatures(combined.list[[2]])[1:3000]
comb.features <- intersect(macaque.vf, mouse.vf)

BiocManager::install("renv")
renv::init()
renv::install('reticulate')
library(reticulate)
renv::install('Seurat')
renv::install('ggplot2')
renv::install('png')
renv::install('BiocManager')
renv::install('devtools')
renv::install('cowplot')
renv::install("sctransform")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("multtest")
BiocManager::install("glmGamPoi")

devtools::install_github('cellgeni/sceasy')


# Installing Python requirements

reticulate::py_install('scanpy')
reticulate::py_install('python-igraph')
reticulate::py_install('louvain')

renv::snapshot()
install.packages("reticulate")
library(ggplot2)
library(Seurat)
library(SeuratData)
library(SingleCellExperiment)
library(cowplot)
library(sceasy)
library(scater)
library(glmGamPoi)
sce <- merge(x=seurat_macaque_gc, y=seurat_mouse_gc)
table(sce$celltype)
scvi.int <- sce
scvi.int <- scvi.int[comb.features]

adata <- convertFormat(scvi.int, from="seurat", to="anndata", 
                       main_layer="counts", drop_single_values=FALSE)
print(adata)
adata$obs$head()
class(adata$obs)
head(py_to_r(adata$obs))
class(adata)

library(reticulate)
library(sceasy)

sc <- import("scanpy", convert = FALSE)

py_install("pandas")
conda_install("r-reticulate", "scipy")
conda_install("r-reticulate", "scvi-tools")

scvi <- import("scvi", convert = FALSE)
# run setup_anndata
scvi$model$SCVI$setup_anndata(adata)



