library(Seurat)
library(Matrix)
library(readr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(limma)

#In order to integrate species data, they must have the same names for 'nFeatures' before merging
#load both species data that is subsetted to TFs and renamed to human ortholog gene IDs.
species1_EDNRcells_TFs <- readRDS("path/to/species1/EDNRcells_TF.RDS")
species1_EDNRcells_TFs

#PmTF_EDNRB_besthit
species2_EDNRcells_TFs <- readRDS("path/to/species2/EDNRcells_TF.RDS")
species2_EDNRcells_TFs

######
##merge objects
#objects must be in same order as identifiers PmTransc/pmNC 
#cell IDs are how you identify where the cells came from
merged_EDNRcells <- merge(species1_EDNRcells_TFs, y = species2_EDNRcells_TFs, add.cell.ids = c("species1 ID", "species2 ID"))
merged_EDNRcells


#####
##preprocess merged Data
merged_EDNRcells <- NormalizeData(merged_EDNRcells)
#601/4 = 150
merged_EDNRcells <- FindVariableFeatures(merged_EDNRcells, selection.method = 'vst', nfeatures = 150)
all.genes.merged <- rownames(merged_EDNRcells)
merged_EDNRcells <- ScaleData(merged_EDNRcells, features = all.genes.merged)
merged_EDNRcells <- RunPCA(merged_EDNRcells, features = VariableFeatures(object = merged_EDNRcells))
ElbowPlot(merged_EDNRcells)
merged_EDNRcells <- FindNeighbors(merged_EDNRcells, dims = 1:10)
merged_EDNRcells <- FindClusters(merged_EDNRcells)
merged_EDNRcells <- RunUMAP(merged_EDNRcells, dims = 1:10)
DimPlot(merged_EDNRcells, reduction = 'umap', group.by = 'orig.ident')

#####
##Batch Correction video tutorial 11/20/2023
# Create a list to store the modified objects
obj.list <- SplitObject(merged_EDNRcells, split.by = 'orig.ident')

for (i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

#find anchors (CCA)
features <- SelectIntegrationFeatures(object.list = obj.list)


# Find anchors (CCA) with reduced features

anchors_merged <- FindIntegrationAnchors(object.list = obj.list,
                                         anchor.features = features)
#integrate anchors
EDNRcells.integrated <- IntegrateData(anchorset = anchors_merged)

##Scale / PCA / UMAP integrated data
EDNRcells.integrated <- ScaleData(object = EDNRcells.integrated)
EDNRcells.integrated <- RunPCA(object = EDNRcells.integrated)
EDNRcells.integrated 
EDNRcells.integrated <- FindVariableFeatures(EDNRcells.integrated)
EDNRcells.integrated <- FindNeighbors(EDNRcells.integrated, dims = 1:10)
EDNRcells.integrated <- FindClusters(EDNRcells.integrated, resolution = 0.5)
EDNRcells.integrated <- RunUMAP(object = EDNRcells.integrated , dims = 1:10)


#visualize by species
species_plot <- DimPlot(EDNRcells.integrated, reduction = 'umap', group.by = 'orig.ident')
species_plot

#visualize by cluster
cluster_plot <- DimPlot(EDNRcells.integrated, reduction = 'umap', group.by = 'seurat_clusters')
cluster_plot

grid.arrange(species_plot, cluster_plot, ncol = 2, nrow = 1)
View(EDNRcells.integrated@meta.data)

obj.list <- SplitObject(EDNRcells.integrated, split.by = 'orig.ident')
obj.list


saveRDS(EDNRcells.integrated, file = "EDNRcells.integrated.RDS")