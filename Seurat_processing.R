library(Seurat)
library(Matrix)
library(readr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(limma)

# Load the Data Matrices (.tsv & .mtx)
#directory path to files from CellRanger (barcodes.tsv/genes.tsv/matrix.mtx)
allgenes.data <- Read10X(data.dir = "path/to/cellranger/output/")

#create seurat object
allgenes <- CreateSeuratObject(counts = allgenes.data, project = "species_sample", min.cells = 3, min.features = 50)

#Subset to only include cells expressing EDNR (NCCs)
#Subset on the expression level of a gene/feature
#in lamprey reference "EDNRB-COJTA" "EDNRB-MOUSE" = only ednr genes > 0 means at any expression level
EDNRcells <- subset(allgenes, subset = `EDNRB-MOUSE` > 0 | `EDNRB-COTJA`)

#Filter object
VlnPlot(EDNRcells, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
EDNRcells <- subset(EDNRcells, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & nCount_RNA < 20000)
EDNRcells


######
#####SUBSET object EDNRcells to only contain transcription factor gene features with one-to-many orthologs between
# lamprey/zebrafish and human
#convert object into a matrix so 'rownames' are genes within 
EDNRcells.counts <- GetAssayData(object = EDNRcells, assay = "RNA")
rownames(EDNRcells.counts)

# load gene mapping .TSV or .CSV that has the lamprey/zebrafish and human TF orthologs and accession IDs
gene_map <- read.csv(file = "Path/to/mappingfile.csv")

# Select only old and new names and remove duplicates
#lamprey.gene.ID = original gene annotation of transcript
#human.ortholog = human ortholog of lamprey/zebrafish gene determined by OrthoFinder
gene_map <- select(gene_map,lamprey.gene.ID, human.ortholog) %>% unique()

# retain only 1-to-1 orthologs in look-up table (LUT)
gene_map_uniq_original <- group_by(gene_map, lamprey.gene.ID) %>% filter(n() == 1)
gene_map_1to1 <- group_by(gene_map_uniq_original, human.ortholog) %>% filter(n() == 1)

# Identify genes not present in ENDRcells
genes_not_present <- setdiff(gene_map_1to1$lamprey.gene.ID, rownames(EDNRcells.counts))

# Subset gene_map_1to1 to exclude genes_not_present
gene_map_filtered <- subset(gene_map_1to1, !(lamprey.gene.ID %in% genes_not_present))
EDNRcells_TFs <- EDNRcells.counts[gene_map_filtered$lamprey.gene.ID, ]

# Confirm that rownames(gene names) are in same order as LUT
all.equal(gene_map_filtered$lamprey.gene.ID, rownames(EDNRcells_TFs))
#TRUE= yes row names equal original names on list.
#if FALSE = check gene names match between mapping file and R script

# Rename genes using LUT
EDNRcells_TFs.data <- EDNRcells_TFs
rownames(EDNRcells_TFs.data) <- gene_map_filtered$human.ortholog
rownames(EDNRcells_TFs.data)


#create seruat object of EDNRB expressing cells with ONLY TF orthos within the data
EDNRcells_TFs <- CreateSeuratObject(counts = EDNRcells_TFs.data, project = "species_sample_EDNR_TFs", min.cells = 3)
EDNRcells_TFs

##PREPROCESSING
EDNRcells_TFs <- NormalizeData(EDNRcells_TFs)
saveRDS(EDNRcells_TFs, file = "EDNRcells_TFs_object.RDS")
#An object of class Seurat 

######
# SUBSET object EDNRcells to only contain transcription factor gene features with 1-to-1 orthologs between
# lamprey/zebrafish and human\

# Select only old and new names and remove duplicates
gene_map <- select(gene_map, original_name, new_name) %>% unique()

# View genes with 2+ mappings to reassign later
group_by(gene_map, original_name) %>% filter(n() > 1) %>% arrange(original_name) %>% View()

# retain only 1-to-1 orthologs in lookup table (LUT)
gene_map_uniq_original <- group_by(gene_map, original_name) %>% filter(n() == 1)
gene_map_1to1 <- group_by(gene_map_uniq_original, new_name) %>% filter(n() == 1)


# Identify genes not present in PmTransc_EDNRB.counts
genes_not_present <- setdiff(gene_map_1to1$original_name, rownames(EDNRcells.counts))

# Subset gene_map_1to1 to exclude genes_not_present
gene_map_filtered <- subset(gene_map_1to1, !(original_name %in% genes_not_present))

EDNRcells_1to1_TFs.data <- EDNRcells.counts[gene_map_filtered$original_name, ]

# Confirm that rownames(gene names) are in same order as LUT
all.equal(gene_map_filtered$original_name, rownames(EDNRcells_1to1_TFs.data))
#TRUE= yes row names equal original names on list. 

#rename 
rownames(EDNRcells_1to1_TFs.data) <- gene_map_filtered$new_name
rownames(EDNRcells_1to1_TFs.data)

#create seruat object of EDNRB expressing cells with ONLY TF orthos within the data
EDNRcells_1to1_TFs <- CreateSeuratObject(counts = EDNRcells_1to1_TFs.data, project = "species_sample_EDNR_1to1_TFs", min.cells = 3)
EDNRcells_1to1_TFs

#PREPROCESSING
EDNRcells_1to1_TFs <- NormalizeData(EDNRcells_1to1_TFs)
saveRDS(EDNRcells_1to1_TFs, file = "EDNRcells_1to1_TFs.RDS")


#####
#Same commands can be run on the original, 'allgenes' object or 'EDNRcells_1to1_TFs'
#if nFeatures >10,000 then set vst = 2000
#when TFs only, set vst = (#rownames/4) or (1/4th nfeature) 
EDNRcells_TFs <- FindVariableFeatures(EDNRcells_TFs, selection.method = "vst", nfeatures = XXX)

#scale data
all.genes.EDNR <- rownames(EDNRcells_TFs)
EDNRcells_TFs <- ScaleData(EDNRcells_TFs, features = all.genes.EDNR)

#PCA
EDNRcells_TFs <- RunPCA(EDNRcells_TFs, features = VariableFeatures(object = EDNRcells_TFs))

#cluster/neighbor
EDNRcells_TFs <- FindNeighbors(EDNRcells_TFs, dims = 1:10)
# If processing an object featuring ALL GENES (>10,000) use resolution of 0.5
EDNRcells_TFs <- FindClusters(EDNRcells_TFs, resolution = 0.8)

###UMAP
EDNRcells_TFs <- RunUMAP(EDNRcells_TFs, dims = 1:10)
DimPlot(EDNRcells_TFs, reduction = "umap")

saveRDS(EDNRcells_TFs, file = "EDNRcells_TFs_PCA_UMAP.RDS")

