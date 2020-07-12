
rm(list=ls())
gc()
options(stringsAsFactors=FALSE)
# load packages
library(dplyr)
library(export)
library(Seurat)
library(cowplot)
library(tidyverse)
library(reticulate)
library(scater)
library(scran)

# load data
if(!file.exists("GSE102130_gSet.Rdata")){
  gSet <- read.table("GSE102130_rawdata.txt", header = T)
  row.names(gSet) <- gSet[ , 1]
  gSet <- gSet[ , -1]
  save(gSet, file = "GSE102130_gSet.Rdata")
}
load("GSE102130_gSet.Rdata")

# create a Seurat objective
sce <- CreateSeuratObject(
  counts = as.matrix(gSet),
  min.cells = 3,
  min.features = 2500,
  assay = "scRNA",
  project="DIPG")
dim(sce)

# # cell QC
# VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
# plot1 <- FeatureScatter(scset, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(scset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# scset <- subset(scset, subset = nFeature_RNA > 5000 & nFeature_RNA < 12500 & percent.mt < 5)
# dim(scset)

# Normalizing the data
scset <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
scset <- FindVariableFeatures(scset, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scset), 10)
# Plot varible features with and without labels
plot1 <- VariableFeaturePlot(scset)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

# Scaling the data
all.genes <- rownames(scset)
scset <- ScaleData(scset, features = all.genes)

# Perform linear dimensional reduction
scset <- RunPCA(scset, features = VariableFeatures(object = scset))
print(scset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scset, dims = 1:2, reduction = "pca")
DimPlot(scset, reduction = "pca")
DimHeatmap(scset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scset, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the dimensionality of the dataset
scset <- JackStraw(scset, num.replicate = 100)
scset <- ScoreJackStraw(scset, dims = 1:20)
JackStrawPlot(scset, dims = 1:15)
ElbowPlot(scset)

# Cluster the cells
scset <- FindNeighbors(scset, dims = 1:15)
scset <- FindClusters(scset, resolution = 0.5)
head(Idents(scset), 6)

# Run non-linear dimensional reduction(UMAP/tSNE)
scset <- RunUMAP(scset, dims = 1:10)
DimPlot(scset, reduction = "umap", pt.size = 1)
scset <- RunTSNE(scset, dims = 1:10)
DimPlot(scset, reduction = "tsne")

# Find cluster biomarker
scset.markers <- FindAllMarkers(scset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
FeaturePlot(scset, features = c("Krt8", "Ngp", "Ighm", "Pld4", "Hba-a1"))

# Save the data
save(scset, scset.markers, file = "GSE102130_gSet_saving.Rdata")
load("GSE102130_saving.Rdata")
