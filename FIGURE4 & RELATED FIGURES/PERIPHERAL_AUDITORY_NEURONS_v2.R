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



################RENAME MOUSE SGN DIVIDED INTO TYPE I II AND ia ib ic
# substract spiral ganglion neurons from mouse
# integrate our data with P28 published data

seurat_mouse_mod <- readRDS("F:/PROJECTS/PROJECT_MONKEY SC-RNA SEQ/WORKPLACE/R/seurat_mouse_mod.rds")
sce=seurat_mouse_mod

Idents(sce) <- "subclass_label"
table(sce$subclass_label)
sgn_mouse <- subset(sce, idents=c("Type I","Type II"))
sgn_mouse

table(sgn_mouse$orig.ident, sgn_mouse$subclass_label)
# re-cluster glia cells by joint CCA-embedding
DefaultAssay(sgn_mouse) <- "RNA"
sgn_mouse$species <- "mouse"
combined.list <- SplitObject(sgn_mouse, split.by = "orig.ident")
for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], vars.to.regress = "nCount_RNA",
                                    verbose = TRUE, method="glmGamPoi")
}


library(scrattch.hicat)
table(sgn_mouse$cluster_label)
table(sgn_mouse$seurat_clusters)
Var.genes.mouse_p25_1 <- select_markers(combined.list$Mouse_P25_1@assays$SCT@counts, 
                                       combined.list$Mouse_P25_1$seurat_clusters, 
                                       n.markers = 100)
Var.genes.mouse_p25_1.markers <- Var.genes.mouse_p25_1$markers

Var.genes.mouse_p25_2 <- select_markers(combined.list$Mouse_P25_2@assays$SCT@counts, 
                                        combined.list$Mouse_P25_2$seurat_clusters, 
                                        n.markers = 100)
Var.genes.mouse_p25_2.markers <- Var.genes.mouse_p25_2$markers

Var.genes.mouse_p28_1 <- select_markers(combined.list$Mouse_P28_1@assays$SCT@counts, 
                                        combined.list$Mouse_P28_1$seurat_clusters, 
                                        n.markers = 100)
Var.genes.mouse_p28_1.markers <- Var.genes.mouse_p28_1$markers

Var.genes.mouse_p28_2 <- select_markers(combined.list$Mouse_P28_2@assays$SCT@counts, 
                                        combined.list$Mouse_P28_2$seurat_clusters, 
                                        n.markers = 100)
Var.genes.mouse_p28_2.markers <- Var.genes.mouse_p28_2$markers

total.Var.genes <- unique(c(Var.genes.mouse_p25_1$markers,Var.genes.mouse_p25_2$markers,
                            Var.genes.mouse_p28_1$markers,Var.genes.mouse_p28_2$markers))


total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Mouse_P25_1@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Mouse_P25_2@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Mouse_P28_1@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$Mouse_P28_2@assays$SCT@counts))]
write.csv(total.Var.genes,file = "features used for cca anchors sgn mouse.csv")

features=total.Var.genes

combined.list <- PrepSCTIntegration(object.list =combined.list,
                                    anchor.features = features)
sgn_mouse.anchors <- FindIntegrationAnchors(object.list = combined.list,
                                             normalization.method = "SCT", 
                                             anchor.features = features,dims = 1:30,
                                             reduction="cca")
sgn_mouse.combined.sct <- IntegrateData(anchorset = sgn_mouse.anchors, 
                                         normalization.method = "SCT",dims = 1:30)
sgn_mouse.combined.sct <- RunPCA(sgn_mouse.combined.sct, features = features, npcs=30)
ElbowPlot(sgn_mouse.combined.sct , ndims = 30)
sgn_mouse.combined.sct <- FindNeighbors(sgn_mouse.combined.sct, 
                                         reduction = "pca", dims = 1:30,nn.eps = 0)
sgn_mouse.combined.sct <- FindClusters(sgn_mouse.combined.sct)

saveRDS(sgn_mouse.combined.sct,file="sgn_mouse_anchors_integrated.rds")
sgn_mouse.combined.sct=sgn_mouse_anchors_integrated

sgn_mouse.combined.sct <- FindNeighbors(sgn_mouse.combined.sct, reduction = "pca",
                                         dims = 1:10)
sgn_mouse.combined.sct <- FindClusters(sgn_mouse.combined.sct,resolution = 0.8)
sgn_mouse.combined.sct <- RunUMAP(sgn_mouse.combined.sct, reduction = "pca", dims = 1:10)

sgn_mouse.combined.sct$seurat_clusters.new <- as.integer(sgn_mouse.combined.sct$seurat_clusters)
# Visualization
p1 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE)
p3 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE,group.by = "cluster_label")
p1+p2+p3
table(sgn_mouse.combined.sct$orig.ident,sgn_mouse.combined.sct$seurat_clusters.new)
table(sgn_mouse.combined.sct$cluster_label,sgn_mouse.combined.sct$seurat_clusters.new)
saveRDS(sgn_mouse.combined.sct,file = "sgn_mouse_combined.rds")
sgn_mouse.combined.sct=sgn_mouse_combined
#########################################################################################
DefaultAssay(sgn_mouse.combined.sct) <- "SCT"
markers.to.plot<-c("NEFL","SNAP25","TUBB3","NEFH",  #TYPE I SGN 
                   "CALB2","MDGA1","RXRG","B3GAT1",  #IA
                   "SEMA3E", "CALB1","NTNG1","LRRC52",    #IB
                   "RUNX1", "POU4F1","LYPD1","GRM8",            #IC
                   "GATA3","SMAD6","ATP2B4","TH","PRPH" #TYPE II
                  
                  
                   
                   
)
VlnPlot(mouse, features = markers.to.plot,group.by = "celltype")
VlnPlot(sgn_mouse.combined.sct, features = markers.to.plot,group.by = "cluster_label")
VlnPlot(sgn_mouse.combined.sct, features = "KDM5B",group.by = "cluster_label")
DotPlot(sgn_mouse.combined.sct, features = markers.to.plot, dot.scale = 8,
        group.by = "cluster_label",
        cols  =c("white", "#ad9300")) + RotatedAxis()
# Visualization
p1 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", 
              group.by = "orig.ident",pt.size = 2)
p2 <- DimPlot(sgn_mouse.combined.sct, reduction = "umap", label = TRUE,
              group.by = "cluster_label",pt.size = 2)
p2+p1

FeaturePlot(sgn_mouse.combined.sct,features = "TLE4")
saveRDS(sgn_mouse.combined.sct, file = "P25_SGN_MOUSE.rds")
P25_SGN_MOUSE <- readRDS("F:/PROJECTS/PROJECT_MONKEY SC-RNA SEQ/WORKPLACE/R/P25_SGN_MOUSE.rds")
#rename celltype
Idents(P25_SGN_MOUSE) <- "cluster_label"
table(Idents(P25_SGN_MOUSE))
p1 <- DimPlot(P25_SGN_MOUSE, reduction = "umap", 
              group.by = "orig.ident",pt.size = 2)
p2 <- DimPlot(P25_SGN_MOUSE, reduction = "umap", label = TRUE,
              group.by = "cluster_label",pt.size = 2)
p2+p1
new.cluster.ids <- c("TypeII_mouse","TypeIC_mouse",
                     "TypeIA_mouse","TypeIB_mouse")
names(new.cluster.ids) <- levels(P25_SGN_MOUSE)
P25_SGN_MOUSE<- RenameIdents(P25_SGN_MOUSE, new.cluster.ids)
DimPlot(P25_SGN_MOUSE, reduction = "umap",label = TRUE,repel = TRUE)
P25_SGN_MOUSE$celltype_species <- Idents(P25_SGN_MOUSE)

library(sceasy)
library(reticulate)
use_condaenv('EnvironmentName')
sce=P25_SGN_MOUSE
DefaultAssay(sce) <- "RNA"
sceasy::convertFormat(sce, from="seurat", to="anndata",
                      outFile='Mm_sgn_python.h5ad')
# we used scvi for integration.
#########################################################################################
#substract sgn from macaque cochlea
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
adata_loom <- connect(filename = "macaque_sgn2.loom",
                      mode = "r+",skip.validate = TRUE)
matrix=adata_loom[["matrix"]][,]
matrix=t(matrix)
dim(matrix)
gene = adata_loom$row.attrs$var_names[]
barcode = adata_loom$col.attrs$obs_names[]

meta_data = read.csv('macaque_sgn_obs2.csv',row.names = 1) # as form as dataframe format
meta_feature = read.csv('macaque_sgn_var2.csv',row.names = 1)
colnames(matrix)= barcode
row.names(matrix)= gene
sgn_macaque= CreateSeuratObject(counts = matrix,meta.data = meta_data,
                               
                               min.cells = 0, 
                               min.features = 0)

sgn_macaque@assays[["RNA"]]@meta.features <- meta_feature
table(sgn_macaque$age, sgn_macaque$cell_type)

# re-cluster sgn by joint CCA-embedding
DefaultAssay(sgn_macaque) <- "RNA"
sgn_macaque$species <- "macaque"
# RENAME CLUSTERS WITH "celltype'
Idents(sgn_macaque) <- "age"
table(sgn_macaque$age)
new.cluster.ids <- c("Macaque_1Y", "Macaque_11Y","Macaque_5Y")
names(new.cluster.ids) <- levels(sgn_macaque)
sgn_macaque<- RenameIdents(sgn_macaque, new.cluster.ids)
sgn_macaque$age <- Idents(sgn_macaque)
table(Idents(sgn_macaque))
sgn_macaque <- SCTransform(sgn_macaque, verbose = T, 
                           vars.to.regress = c("nCount_RNA"), 
                          conserve.memory = T)
sgn_macaque <- RunPCA(sgn_macaque,npcs = 50)
DimPlot(sgn_macaque,reduction = "pca")
ElbowPlot(sgn_macaque , ndims = 50)
sgn_macaque <- FindNeighbors(sgn_macaque, reduction = "pca", dims = 1:30)
sgn_macaque <- FindClusters(sgn_macaque)
sgn_macaque <- RunUMAP(sgn_macaque, reduction = "pca", dims = 1:10)
# Visualization
DimPlot(sgn_macaque, reduction = "umap",label = TRUE,group.by = 'age')



#########################################################################################

####################################################################################
# sgn from macaque can be divided into 4 subtypes that is remarked from mouse.
table(P25_SGN_MOUSE$celltype_species)
table(sgn_macaque$celltype_species)


# RENAME CLUSTERS WITH "celltype'
Idents(sgn_macaque) <- "celltype_species"
table(sgn_macaque$celltype_species)
new.cluster.ids <- c("TypeIB", "TypeIA","TypeIC","TypeII")
names(new.cluster.ids) <- levels(sgn_macaque)
sgn_macaque<- RenameIdents(sgn_macaque, new.cluster.ids)
sgn_macaque$celltype <- Idents(sgn_macaque)
table(Idents(sgn_macaque))

# RENAME CLUSTERS WITH "celltype'
Idents(sgn_macaque) <- "celltype"
table(sgn_macaque$celltype)
new.cluster.ids <- c("TypeI", "TypeI","TypeI","TypeII")
names(new.cluster.ids) <- levels(sgn_macaque)
sgn_macaque<- RenameIdents(sgn_macaque, new.cluster.ids)
sgn_macaque$subclass <- Idents(sgn_macaque)
table(Idents(sgn_macaque))


# RENAME CLUSTERS WITH "celltype'
Idents(P25_SGN_MOUSE) <- "celltype_species"
table(P25_SGN_MOUSE$celltype_species)
new.cluster.ids <- c("TypeII", "TypeIC","TypeIA","TypeIB")
names(new.cluster.ids) <- levels(P25_SGN_MOUSE)
P25_SGN_MOUSE<- RenameIdents(P25_SGN_MOUSE, new.cluster.ids)
P25_SGN_MOUSE$celltype <- Idents(P25_SGN_MOUSE)
table(Idents(P25_SGN_MOUSE))

# RENAME CLUSTERS WITH "celltype'
Idents(P25_SGN_MOUSE) <- "celltype"
table(P25_SGN_MOUSE$celltype)
new.cluster.ids <- c("TypeII", "TypeI","TypeI","TypeI")
names(new.cluster.ids) <- levels(P25_SGN_MOUSE)
P25_SGN_MOUSE<- RenameIdents(P25_SGN_MOUSE, new.cluster.ids)
P25_SGN_MOUSE$subclass <- Idents(P25_SGN_MOUSE)
table(Idents(P25_SGN_MOUSE))





###############################################################################
#DEG calculate by tt-test in scanpy

macaque_genes <- read.csv(file = "macaque_sgn_subclass_marker_v2_final.csv",row.names = 1)
mouse_genes <- read.csv(file = "mouse_sgn_subclass_marker_v2_final.csv",row.names = 1)
macaque_genes$genes=rownames(macaque_genes)
mouse_genes$genes=rownames(mouse_genes)
all.genes <- data.frame(genes = unique(c(macaque_genes$gene, mouse_genes$gene)))


all.genes$Macaque <- as.character(match(all.genes$genes, macaque_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))

all.genes$Macaque[which(is.na(all.genes$Macaque))] <- FALSE 
all.genes$Macaque[which(all.genes$Macaque != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE

all.genes$Macaque <- as.logical(all.genes$Macaque)
all.genes$Mouse <- as.logical(all.genes$Mouse)
write.csv(all.genes,file="DEGS_sgn_type2_macaque VS mouse by tt-test_scanpy.csv")
library(eulerr)
plot(euler(
  all.genes[ ,2:3]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0("macaque vs. mouse"),
  fills = c("royalblue1", "sienna2")
)


##### plot GO TERMS RELATED TO PBX3
library(ggplot2)
df <- read.csv("GO_TERMS_RELATED_TO_PBX3_V2.csv", header = T)
df$LogP <- -df$negative_log10_of_adjusted_p_value

df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df, 
       aes(LogP, reorder(term_name,LogP))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#FE8D3C",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = term_name),
            size=6, 
            hjust =1.02)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(qvalue)")+
  ggtitle("Enrichment")+
  scale_x_continuous(expand = c(0,0))
###

VlnPlot(sce, features = "PBX3", group.by = "celltype",split.by = "species",
        assay = "SCT")

DotPlot(sce, features = "PBX3", dot.scale = 8,group.by = "celltype_species",assay = "SCT",
        cols  =c("white", "#2141B5")) + RotatedAxis()
####===================type II==========================================
macaque_genes <- read.csv(file = "macaque_sgn_subclass_marker_v2_final.csv",row.names = 1)
mouse_genes <- read.csv(file = "mouse_sgn_subclass_marker_v2_final.csv",row.names = 1)
macaque_genes$genes=rownames(macaque_genes)
mouse_genes$genes=rownames(mouse_genes)
macaque_genes <- macaque_genes[which(macaque_genes$cluster == "TypeII"), ]#chang to 'Type I' 
mouse_genes <- mouse_genes[which(mouse_genes$cluster == "TypeII"), ]

all.genes <- data.frame(genes = unique(c( macaque_genes$genes, mouse_genes$genes)))

all.genes$Macaque <- as.character(match(all.genes$genes, macaque_genes$genes))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$genes))

all.genes$Macaque[which(is.na(all.genes$Macaque))] <- FALSE 
all.genes$Macaque[which(all.genes$Macaque != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE

all.genes$Macaque <- as.logical(all.genes$Macaque)
all.genes$Mouse <- as.logical(all.genes$Mouse)
write.csv(all.genes,file="DEGS_TYPEI_ACROSS_MACAQUE_MOUSE_TT-TEST.csv")
library(eulerr)
plot(euler(
  all.genes[ ,2:3]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0("TypeII:Macaque vs Mouse"),
  fills = c( "maroon4", "sienna2")
)




########################## scatter plots               #####################################################################
############################################################################################################################

#use this code to change the default log expression to log2 expression
macaque=sgn_macaque
mouse=P25_SGN_MOUSE
macaque@assays$RNA@data <- macaque@assays$RNA@counts
macaque@assays$RNA@data@x <- log2(macaque@assays$RNA@data@x + 1)
mouse@assays$RNA@data <- mouse@assays$RNA@counts
mouse@assays$RNA@data@x <- log2(mouse@assays$RNA@data@x + 1)


macaque_mat <- macaque@assays$RNA@data
colnames(macaque_mat) <- macaque$subclass
mouse_mat <- mouse@assays$RNA@data
colnames(mouse_mat) <- mouse$subclass

#find avg expression for TYPEII and TYPEI (log(umi + 1))

macaque_avg <- data.frame(typeII = rowMeans(macaque_mat[ , 
                                which(colnames(macaque_mat) == "TypeII")]))
macaque_avg$typeI <- rowMeans(macaque_mat[ , 
                   which(colnames(macaque_mat) == "TypeI")])
mouse_avg <- data.frame(typeII = rowMeans(mouse_mat[ , 
    which(colnames(mouse_mat) == "TypeII")]))
mouse_avg$typeI <- rowMeans(mouse_mat[ , 
            which(colnames(mouse_mat) == "TypeI")])

#plot scatter
conserved_genes <- all.genes$genes[which( 
                                         all.genes$Macaque ==T & 
                                           all.genes$Mouse == T)]

macaque_avg$conserved <- F
macaque_avg$conserved[which(rownames(macaque_avg) %in% conserved_genes)] <- T

mouse_avg$conserved <- F
mouse_avg$conserved[which(rownames(mouse_avg) %in% conserved_genes)] <- T

macaque_type2_specified <- all.genes$genes[which(all.genes$Macaque==T &
                                                   all.genes$Mouse==F)]

mouse_type2_specified <- all.genes$genes[which(all.genes$Macaque==F &
                                                 all.genes$Mouse==T)]

tf_genelist <- read.csv("human_tf_gene_list.csv")
write.csv(conserved_genes,file = "sgn_type2_conserved_genes_v2.csv")

type2_conserved_genes <- read.csv(file="sgn_type2_conserved_genes_v2.csv",
                                  row.names = 1)
type2_tf_conserved_both <- merge(x=tf_genelist,y=type2_conserved_genes,
                                 by="HGNC.symbol")

write.csv(macaque_type2_specified,file = "sgn_type2_macaque_enriched_genes_v2.csv")
type2_macaque_genes <- read.csv(file="sgn_type2_macaque_enriched_genes_v2.csv",
                                row.names = 1)

type2_tf_macaque <- merge(x=tf_genelist,y=type2_macaque_genes,by="HGNC.symbol")
write.csv(mouse_type2_specified,file = "sgn_type2_mouse_enriched_genes_v2.csv")

type2_mouse_genes <- read.csv(file="sgn_type2_mouse_enriched_genes.csv",
                              row.names = 1)
type2_tf_mouse <- merge(x=tf_genelist,y=type2_mouse_genes,by="HGNC.symbol")




macaque_avg$label_tf_conserved <- NA
macaque_avg$label_tf_conserved[which(rownames(macaque_avg) %in% type2_tf_conserved_both$HGNC.symbol)]<- rownames(macaque_avg)[which(rownames(macaque_avg) %in% type2_tf_conserved_both$HGNC.symbol)]
mouse_avg$label_tf_conserved <- NA
mouse_avg$label_tf_conserved[which(rownames(mouse_avg) %in% type2_tf_conserved_both$HGNC.symbol)]<- rownames(mouse_avg)[which(rownames(mouse_avg) %in% type2_tf_conserved_both$HGNC.symbol)]

macaque_avg$label_tf_macaque <- NA
macaque_avg$label_tf_macaque[which(rownames(macaque_avg) %in% type2_tf_macaque$HGNC.symbol)]<- rownames(macaque_avg)[which(rownames(macaque_avg) %in% type2_tf_macaque$HGNC.symbol)]
mouse_avg$label_tf_macaque <- NA
mouse_avg$label_tf_macaque[which(rownames(mouse_avg) %in% type2_tf_macaque$HGNC.symbol)]<- rownames(mouse_avg)[which(rownames(mouse_avg) %in% type2_tf_macaque$HGNC.symbol)]

macaque_avg$label_tf_mouse <- NA
macaque_avg$label_tf_mouse[which(rownames(macaque_avg) %in% type2_tf_mouse$HGNC.symbol)]<- rownames(macaque_avg)[which(rownames(macaque_avg) %in% type2_tf_mouse$HGNC.symbol)]
mouse_avg$label_tf_mouse <- NA
mouse_avg$label_tf_mouse[which(rownames(mouse_avg) %in% type2_tf_mouse$HGNC.symbol)]<- rownames(mouse_avg)[which(rownames(mouse_avg) %in% type2_tf_mouse$HGNC.symbol)]

library(ggrepel)
library(readxl)

macaque_avg$tf <- F
macaque_avg$tf[which(rownames(macaque_avg) %in% tf_genelist$HGNC.symbol)] <- T
mouse_avg$tf <- F
mouse_avg$tf[which(rownames(mouse_avg) %in% tf_genelist$HGNC.symbol)] <- T

colnames(macaque_avg) <- paste("macaque", colnames(macaque_avg), sep = "_")
colnames(mouse_avg) <- paste("mouse", colnames(mouse_avg), sep = "_")

macaque_avg$genes <- rownames(macaque_avg)
mouse_avg$genes <- rownames(mouse_avg)

to_plot <- merge(x=macaque_avg, y=mouse_avg,by="genes",all.x=TRUE)

# conserved_tfs for type2 vs type 1
p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == T), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & is.na(to_plot$macaque_label_tf_conserved) == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("macaque_typeII (log2 expression)") +
  ylab("macaque_typeI (log2 expression)") +
  xlim(c(0, 2)) +
  ylim(c(0, 2)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$macaque_conserved == T & is.na(to_plot$macaque_label_tf_conserved) == F), ], aes(x = macaque_typeII, y = macaque_typeI, label = macaque_label_tf_conserved),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p2
dev.off()
p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == T), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label_tf_conserved) == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("mouse_typeII (log2 expression)") +
  ylab("mouse_typeI (log2 expression)") +
  xlim(c(0, 5)) +
  ylim(c(0, 5)) +
  theme(aspect.ratio = 1)
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label_tf_conserved) == F), ], aes(x = mouse_typeII, y = mouse_typeI, label = mouse_label_tf_conserved),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p3
### plot macaque specific tf

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == T), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("macaque_typeII (log2 expression)") +
  ylab("macaque_typeI (log2 expression)") +
  xlim(c(0, 2)) +
  ylim(c(0, 2)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], aes(x = macaque_typeII, y = macaque_typeI, label = macaque_label_tf_macaque),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p2
dev.off()
p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == T), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & is.na(to_plot$mouse_label_tf_mouse) == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("mouse_typeII (log2 expression)") +
  ylab("mouse_typeI (log2 expression)") +
  xlim(c(0, 5)) +
  ylim(c(0, 5)) +
  theme(aspect.ratio = 1)
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$mouse_conserved == F & is.na(to_plot$mouse_label_tf_mouse) == F), ], aes(x = mouse_typeII, y = mouse_typeI, label = mouse_label_tf_mouse),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p3


### TYPE2 

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = mouse_typeII), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == T), ], mapping = aes(x = macaque_typeII, y = mouse_typeII), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = mouse_typeII), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], mapping = aes(x = macaque_typeII, y = mouse_typeII), color = "magenta", size = 1.5, shape = 15) +
  xlab("macaque_typeII (log2 expression)") +
  ylab("mouse_typeII (log2 expression)") +
  xlim(c(0, 2)) +
  ylim(c(0, 5)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], aes(x = macaque_typeII, y = mouse_typeII, label = macaque_label_tf_macaque),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 2)
p2
dev.off()
# WE CONSIDER IT IS NOT REASONABLE AND DISGARD THIS PLOT.



####===================type I==========================================
macaque_genes <- read.csv(file = "macaque_sgn_subclass_marker_v2_final.csv",row.names = 1)
mouse_genes <- read.csv(file = "mouse_sgn_subclass_marker_v2_final.csv",row.names = 1)
macaque_genes$genes=rownames(macaque_genes)
mouse_genes$genes=rownames(mouse_genes)
macaque_genes <- macaque_genes[which(macaque_genes$cluster == "TypeI"), ]#chang to 'Type I' 
mouse_genes <- mouse_genes[which(mouse_genes$cluster == "TypeI"), ]

all.genes <- data.frame(genes = unique(c( macaque_genes$genes, mouse_genes$genes)))

all.genes$Macaque <- as.character(match(all.genes$genes, macaque_genes$genes))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$genes))

all.genes$Macaque[which(is.na(all.genes$Macaque))] <- FALSE 
all.genes$Macaque[which(all.genes$Macaque != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE

all.genes$Macaque <- as.logical(all.genes$Macaque)
all.genes$Mouse <- as.logical(all.genes$Mouse)
write.csv(all.genes,file="DEGS_TYPEI_ACROSS_MACAQUE_MOUSE_TT-TEST.csv")
library(eulerr)
plot(euler(
  all.genes[ ,2:3]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0("TypeI:Macaque vs Mouse"),
  fills = c( "maroon4", "sienna2")
)
##### plot GO TERMS RELATED TO PBX3
library(ggplot2)
df <- read.csv("GO_TERMS_RELATED_TO_PBX3.csv", header = T)
df$LogP <- -df$negative_log10_of_adjusted_p_value

df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)

ggplot(data = df, 
       aes(LogP, reorder(term_name,LogP))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#FE8D3C",
           width = 0.8) + 
  geom_text(aes(x=labelx,
                y=labely,
                label = term_name),
            size=6, 
            hjust =1.02)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(qvalue)")+
  ggtitle("Enrichment")+
  scale_x_continuous(expand = c(0,0))
###

VlnPlot(sce, features = "PBX3", group.by = "celltype",split.by = "species",
        assay = "SCT")

DotPlot(sce, features = "PBX3", dot.scale = 8,group.by = "celltype.species",assay = "SCT",
        cols  =c("white", "#2141B5")) + RotatedAxis()

########################## scatter plots               #####################################################################
############################################################################################################################

#use this code to change the default log expression to log2 expression

macaque@assays$RNA@data <- macaque@assays$RNA@counts
macaque@assays$RNA@data@x <- log2(macaque@assays$RNA@data@x + 1)
mouse@assays$RNA@data <- mouse@assays$RNA@counts
mouse@assays$RNA@data@x <- log2(mouse@assays$RNA@data@x + 1)


macaque_mat <- macaque@assays$RNA@data
colnames(macaque_mat) <- macaque$subclass
mouse_mat <- mouse@assays$RNA@data
colnames(mouse_mat) <- mouse$subclass

#find avg expression for TYPEII and TYPEI (log(umi + 1))

macaque_avg <- data.frame(typeII = rowMeans(macaque_mat[ , which(colnames(macaque_mat) == "TypeII")]))
macaque_avg$typeI <- rowMeans(macaque_mat[ , which(colnames(macaque_mat) == "TypeI")])
mouse_avg <- data.frame(typeII = rowMeans(mouse_mat[ , which(colnames(mouse_mat) == "TypeII")]))
mouse_avg$typeI <- rowMeans(mouse_mat[ , which(colnames(mouse_mat) == "TypeI")])

#plot scatter
conserved_genes <- all.genes$genes[which( 
  all.genes$Macaque ==T & 
    all.genes$Mouse == T)]

macaque_avg$conserved <- F
macaque_avg$conserved[which(rownames(macaque_avg) %in% conserved_genes)] <- T

mouse_avg$conserved <- F
mouse_avg$conserved[which(rownames(mouse_avg) %in% conserved_genes)] <- T

macaque_type1_specified <- all.genes$genes[which(all.genes$Macaque==T &
                                                   all.genes$Mouse==F)]

mouse_type1_specified <- all.genes$genes[which(all.genes$Macaque==F &
                                                 all.genes$Mouse==T)]

tf_genelist <- read.csv("human_tf_gene_list.csv")
write.csv(conserved_genes,file = "sgn_type1_conserved_genes_v2.csv")
type1_conserved_genes <- read.csv(file="sgn_type1_conserved_genes_v2.csv",row.names = 1)
type1_tf_conserved_both <- merge(x=tf_genelist,y=type1_conserved_genes,by="HGNC.symbol")
write.csv(macaque_type1_specified,file = "sgn_type1_macaque_enriched_genes_v2.csv")
type1_macaque_genes <- read.csv(file="sgn_type1_macaque_enriched_genes_v2.csv",row.names = 1)
type1_tf_macaque <- merge(x=tf_genelist,y=type1_macaque_genes,by="HGNC.symbol")
write.csv(mouse_type1_specified,file = "sgn_type1_mouse_enriched_genes_v2.csv")
type1_mouse_genes <- read.csv(file="sgn_type1_mouse_enriched_genes_v2.csv",row.names = 1)
type1_tf_mouse <- merge(x=tf_genelist,y=type1_mouse_genes,by="HGNC.symbol")

macaque_avg$label_tf_conserved <- NA
macaque_avg$label_tf_conserved[which(rownames(macaque_avg) %in% type1_tf_conserved_both$HGNC.symbol)]<- rownames(macaque_avg)[which(rownames(macaque_avg) %in% type1_tf_conserved_both$HGNC.symbol)]
mouse_avg$label_tf_conserved <- NA
mouse_avg$label_tf_conserved[which(rownames(mouse_avg) %in% type1_tf_conserved_both$HGNC.symbol)]<- rownames(mouse_avg)[which(rownames(mouse_avg) %in% type1_tf_conserved_both$HGNC.symbol)]

macaque_avg$label_tf_macaque <- NA
macaque_avg$label_tf_macaque[which(rownames(macaque_avg) %in% type1_tf_macaque$HGNC.symbol)]<- rownames(macaque_avg)[which(rownames(macaque_avg) %in% type1_tf_macaque$HGNC.symbol)]
mouse_avg$label_tf_macaque <- NA
mouse_avg$label_tf_macaque[which(rownames(mouse_avg) %in% type1_tf_macaque$HGNC.symbol)]<- rownames(mouse_avg)[which(rownames(mouse_avg) %in% type1_tf_macaque$HGNC.symbol)]

macaque_avg$label_tf_mouse <- NA
macaque_avg$label_tf_mouse[which(rownames(macaque_avg) %in% type1_tf_mouse$HGNC.symbol)]<- rownames(macaque_avg)[which(rownames(macaque_avg) %in% type1_tf_mouse$HGNC.symbol)]
mouse_avg$label_tf_mouse <- NA
mouse_avg$label_tf_mouse[which(rownames(mouse_avg) %in% type1_tf_mouse$HGNC.symbol)]<- rownames(mouse_avg)[which(rownames(mouse_avg) %in% type1_tf_mouse$HGNC.symbol)]

library(ggrepel)
library(readxl)

macaque_avg$tf <- F
macaque_avg$tf[which(rownames(macaque_avg) %in% tf_genelist$HGNC.symbol)] <- T
mouse_avg$tf <- F
mouse_avg$tf[which(rownames(mouse_avg) %in% tf_genelist$HGNC.symbol)] <- T

colnames(macaque_avg) <- paste("macaque", colnames(macaque_avg), sep = "_")
colnames(mouse_avg) <- paste("mouse", colnames(mouse_avg), sep = "_")

macaque_avg$genes <- rownames(macaque_avg)
mouse_avg$genes <- rownames(mouse_avg)

to_plot <- merge(x=macaque_avg, y=mouse_avg,by="genes",all.x=TRUE)

# conserved_tfs for type2 vs type 1
p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == T), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & is.na(to_plot$macaque_label_tf_conserved) == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("macaque_typeII (log2 expression)") +
  ylab("macaque_typeI (log2 expression)") +
  xlim(c(0, 2)) +
  ylim(c(0, 2)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$macaque_conserved == T & is.na(to_plot$macaque_label_tf_conserved) == F), ], aes(x = macaque_typeII, y = macaque_typeI, label = macaque_label_tf_conserved),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p2
dev.off()
p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == T), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label_tf_conserved) == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("mouse_typeII (log2 expression)") +
  ylab("mouse_typeI (log2 expression)") +
  xlim(c(0, 5)) +
  ylim(c(0, 5)) +
  theme(aspect.ratio = 1)
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label_tf_conserved) == F), ], aes(x = mouse_typeII, y = mouse_typeI, label = mouse_label_tf_conserved),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p3
### plot macaque specific tf

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "grey80", size = 2) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == T), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "cyan2", size = 2.5, shape = 15) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "red", size = 2.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], mapping = aes(x = macaque_typeII, y = macaque_typeI), color = "magenta", size = 4, shape = 15) +
  xlab("macaque_typeII (log2 expression)") +
  ylab("macaque_typeI (log2 expression)") +
  xlim(c(0, 2)) +
  ylim(c(0, 2)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], aes(x = macaque_typeII, y = macaque_typeI, label = macaque_label_tf_macaque),
                           nudge_x = 0.5,
                           nudge_y = 0.2,
                           size = 3)
p2
dev.off()
p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == T), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & is.na(to_plot$mouse_label_tf_mouse) == F), ], mapping = aes(x = mouse_typeII, y = mouse_typeI), color = "magenta", size = 1.5, shape = 15) +
  xlab("mouse_typeII (log2 expression)") +
  ylab("mouse_typeI (log2 expression)") +
  xlim(c(0, 5)) +
  ylim(c(0, 5)) +
  theme(aspect.ratio = 1)
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$mouse_conserved == F & is.na(to_plot$mouse_label_tf_mouse) == F), ], aes(x = mouse_typeII, y = mouse_typeI, label = mouse_label_tf_mouse),
                           nudge_x = 0.5,
                           nudge_y = 0.05,
                           size = 3)
p3


### TYPE1 

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeI, y = mouse_typeI), color = "grey80", size = 1.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & to_plot$macaque_tf == T), ], mapping = aes(x = macaque_typeI, y = mouse_typeI), color = "cyan2", size = 2.5, shape = 15) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == T & to_plot$macaque_tf == F), ], mapping = aes(x = macaque_typeI, y = mouse_typeI), color = "red", size = 1.5) +
  geom_point(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], mapping = aes(x = macaque_typeI, y = mouse_typeI), color = "magenta", size = 2, shape = 15) +
  xlab("macaque_typeI (log2 expression)") +
  ylab("mouse_typeI (log2 expression)") +
  xlim(c(0, 2)) +
  ylim(c(0, 6)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$macaque_conserved == F & is.na(to_plot$macaque_label_tf_macaque) == F), ], aes(x = macaque_typeI, y = mouse_typeI, label = macaque_label_tf_macaque),
                           nudge_x = 0.5,
                           nudge_y = 0.1,
                           size = 2)
p2
dev.off()

