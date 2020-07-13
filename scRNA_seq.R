## 0. Prepare software environment and load libraries
rm(list = ls())
gc()
set.seed(12345)
graphics.off()
options(stringsAsFactors = FALSE)
# load packages
pkgs <- c("dplyr", "Seurat", "cowplot", "tidyverse", "stringr", "reticulate", "ggpubr",
          "scran", "pheatmap", "ggfortify", "readxl", "gplots", "export")
# installpkgs <- function(pkgs) {
#   new.pkgs <- pkgs[!(pkgs %in% installed.packages()[ , "Package"])]
#   if (length(new.pkgs))
#     BiocManager::install(new.pkgs, ask = F, update = F)
#   sapply(pkgs, require, character.only = T)
# }
# installpkgs(pkgs)
lapply(pkgs, library, character.only = T)

## 1. Prepare rawdata to proceed
if(!file.exists("GSE102130_gSet.Rdata")){
  gSet <- read.table("GSE102130_rawdata.txt", header = T)
  row.names(gSet) <- gSet[ , 1]
  gSet <- gSet[ , -1]
  save(gSet, file = "GSE102130_gSet.Rdata")
}
load("GSE102130_gSet.Rdata")
counts <- gSet
# create annotationdata
cellName <- as.character(colnames(counts))
id <- strsplit(cellName, "[.]")
names(id) <- cellName
anno <- matrix(nrow = length(id), ncol = 3)
for (i in 1:length(id)) {
  anno[i, 1] <- cellName[i]
  anno[i, 2] <- id[[i]][1]
  anno[i, 3] <- id[[i]][2]
}
colnames(anno) <- c("cellName", "patientID", "plateID")
rownames(anno) <- anno[ , 1]
save(counts, anno, file = "GSE102130.Rdata")

## 2. Load counts to analyse 
# load data
load("GSE102130.Rdata")
{
  Anno1 <- anno[which(anno[ , 2] == "MUV1"), ]
  Anno5 <- anno[which(anno[ , 2] == "MUV5"), ]
  Anno10 <- anno[which(anno[ , 2] == "MUV10"), ]
  Anno836 <- anno[which(anno[ , 2] == "BCH836"), ]
  Anno869 <- anno[which(anno[ , 2] == "BCH869"), ]
  Anno126 <- anno[which(anno[ , 2] == "BCH1126"), ]
}
Anno <- rbind(Anno1, Anno5, Anno10, Anno836, Anno869, Anno126)
exprset <- counts[ , which(colnames(counts) %in% rownames(Anno))]
save(Anno, exprset, file = "tmp.Rdata")
# filter genes by the paper source: Ea = log2(average(TPMgene)+1) < 4 
ave_row <- rowSums(exprset)/ncol(exprset)
log_TPM <- as.matrix(log2(ave_row + 1))
filter_by_Ea <- as.matrix(log_TPM[which(log_TPM >= 4), ])
exprSet <- exprset[rownames(filter_by_Ea), ]
# create a Seurat object and quick QC
sce <- CreateSeuratObject(counts = as.matrix(exprSet),
                          min.cells = 3, min.features = 200,
                          assay = "scRNA", project="DIPG")
dim(sce)
save(sce, file = "sce_GSE102130.Rdata")

## 3. Single cell data analysis 
# normalizing the data
load("sce_GSE102130.Rdata")
scset <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
# identification of highly variable features
scset <- FindVariableFeatures(scset, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scset), 10)
# plot varible features with and without labels
plot1 <- VariableFeaturePlot(scset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
# scaling the data
all.genes <- rownames(scset)
scset <- ScaleData(scset, features = all.genes)
# perform linear dimensional reduction
scset <- RunPCA(scset, features = VariableFeatures(object = scset))
# examine and visualize PCA results by a few different ways:
# 1). print the first 5 PCs with 5 genes
print(scset[["pca"]], dims = 1:5, nfeatures = 5)
# 2). plot the vizplot of each PC according to features
VizDimLoadings(scset, dims = 1:2, reduction = "pca")
# 3). PCA scatter plot
DimPlot(scset, reduction = "pca")
# 4). heatmap to decide those PCs according to PCA scores
DimHeatmap(scset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scset, dims = 1:20, cells = 500, balanced = TRUE)
# determine the dimensionality of the dataset by a few different ways:
# 1). JackStraw for resampling test to enrich the lowest p-value
scset <- JackStraw(scset, num.replicate = 100)
scset <- ScoreJackStraw(scset, dims = 1:20)
# find the sharp drop-off in significance after the first PCs.
JackStrawPlot(scset, dims = 1:20)
# 2). ElbowPlot to rank of principle components based on the percentage of variance,
# and find the  ¡®elbow¡¯ around PCs that capture the majority of true signals
ElbowPlot(scset)
save(scset, file = "scset_GSE102130.Rdata")

## 4. Downstream analysis
# cluster the cells
scset <- FindNeighbors(scset, dims = 1:10)
scset <- FindClusters(scset, resolution = 0.15)
head(Idents(scset), 6)
# run non-linear dimensional reduction(tSNE/UMAP):
# tSNE plot
scset <- RunTSNE(scset, dims = 1:20)
DimPlot(scset, reduction = "tsne", label = T, label.size = 5)
# UMAP plot
scset <- RunUMAP(scset, dims = 1:10)
DimPlot(scset, reduction = "umap", pt.size = 1)
#-----------------Start from here---------------------------
# Find cluster biomarker
scset.markers <- FindAllMarkers(scset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
FeaturePlot(scset, features = c("Krt8", "Ngp", "Ighm", "Pld4", "Hba-a1"))

# Save the data
save(scset, scset.markers, file = "GSE102130_gSet_saving.Rdata")
load("GSE102130_saving.Rdata")

## load rawdata
if(F){
  load("0114.Rdata")
  scset <- RunTSNE(pbmc_tsne, dims = 1:10)
  save(scset, file = "0222_raw.Rdata")
}

## Run non-linear dimensional reduction(UMAP/tSNE)
load("0222.Rdata")
DimPlot(scset, 
        pt.size = 1,
        reduction = "tsne",
        label = T,
        label.size = 5)

## Preprocess and Identify the cluster
cluster_id <- data.frame(scset@active.ident)
cluster <- cbind(rownames(cluster_id), cluster_id)
sample <- data.frame(matrix(nrow = nrow(cluster), ncol = 1))
for (i in 1: nrow(cluster)) {
  sample[i, 1] <- substr(cluster[i, 1], 1, 7)
}
sample_matrix <- as.matrix(sample)
if(F){
  sample <- gsub(".P0", "", sample_matrix[ , 1])
  sample <- gsub(".P1", "", sample)
  sample <- gsub(".P", "", sample)
  sample <- gsub("BCH869.", "BCH869", sample)
  sample <- gsub("BCH836.", "BCH836", sample)
}
sample_cluster <- as.matrix(sample)
cluster <- cbind(sample_cluster, cluster)
colnames(cluster) <- c("sample", "name", "id")
if(F){save(scset, cluster, file = "0222.Rdata")}

## Define the cluster
MUV1 <- cluster[ which(cluster[ , 1] == "MUV1"), ]
MUV5 <- cluster[ which(cluster[ , 1] == "MUV5"), ]
MUV10 <- cluster[ which(cluster[ , 1] == "MUV10"), ]
BCH836 <- cluster[ which(cluster[ , 1] == "BCH836"), ]
BCH869 <- cluster[ which(cluster[ , 1] == "BCH869"), ]
BCH1126 <- cluster[ which(cluster[ , 1] == "BCH1126"), ]
Oligo <- cluster[ which(cluster[ , 1] == "Oligo"), ]

## Assign markers to clusters
if(F){
  CellCyle <- c("PCNA", "CDK1")
  AC_like <- c("GFAP", "APOE")
  OC_like <- c("MBP", "PLP1")
  OPC_like <- c("PDGFRA", "CSPG4")
  OPC_variable <- c("ITM2C", "SCG3")
}
# par(mfrow = c(2, 2))
VlnPlot(object = scset, sort = T,
        features = c(CellCyle[2], AC_like[2], OC_like[2],
                     OPC_like[1], OPC_variable[2], OC_like[1]))

FeaturePlot(object = scset, label = T, reduction = "tsne", cols = c("grey", "blue"),
            features = c(CellCyle[2], AC_like[2], OC_like[2],
                         OPC_like[1], OPC_variable[2], OC_like[1]))


# Assign the celltype identity to the cluster
DimPlot(scset, pt.size = 1, reduction = "umap", label = T, label.size = 5)
FeaturePlot(object = scset, label = T, reduction = "umap", cols = c("grey", "blue"),
            features = c(CellCyle[2], AC_like[2], OC_like[2],
                         OPC_like[1], OPC_variable[2], OC_like[1]))
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
new.cluster.ids <- c("OPC-like1", "OPC-like2", "Cell Cycle1", "Cell Cycle2", "Normal", "OC-like1",
                     "OPC-variable", "AC-like1", "AC-like2", "OC-like2")
scset@active.ident <- plyr::mapvalues(x = scset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(scset, pt.size = 1, reduction = "umap", label = T, label.size = 3)


## Pseudotime Analysis
load("0222.Rdata")
dipgdata <- as.matrix(scset@assays[["scRNA"]]@counts)
# data preprocess
procdata <- preprocess(dipgdata)
# extract subpopulation
subpopulation <- data.frame(cell = colnames(procdata), 
                            sub = cluster[ , 3])
# perform model-based clustering on expression values
dipgmclust <- exprmclustv(scset)
# dipgmclust[["clusterid"]] <- as.integer(subpopulation[ , 2])
# names(dipgmclust[["clusterid"]]) <- subpopulation[ , 1]
# dipgmclust[["pcareduceres"]][ , 1:2] <- scset@reductions[["umap"]]@cell.embeddings
plotmclustv(dipgmclust, 
            show_full_tree = F, 
            cell_point_size = 1, 
            linetype = 2, arrow_length = 0)
# perform cell-level ordering
dipgorder <- TSCANorder(dipgmclust)
# single-cell pseudotime gene expression
singlegene <- "MBP"
singlegene_expr <- log2(dipgdata[singlegene, ] + 1)
singlegeneplotv(singlegene_expr, dipgorder, cell_size = 1.5, genename = singlegene)




sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.936 
# [2] LC_CTYPE=Chinese (Simplified)_China.936   
# [3] LC_MONETARY=Chinese (Simplified)_China.936
# [4] LC_NUMERIC=C                              
# [5] LC_TIME=Chinese (Simplified)_China.936    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
# [1] httr_1.4.1          tidyr_1.0.0         viridisLite_0.3.0  
# [4] jsonlite_1.6.1      splines_3.6.3       lsei_1.2-0         
# [7] R.utils_2.9.0       leiden_0.3.1        gtools_3.8.1       
# [10] RcppParallel_4.4.4  Rdpack_0.11-0       assertthat_0.2.1   
# [13] ggrepel_0.8.1       globals_0.12.4      pillar_1.4.2       
# [16] lattice_0.20-38     glue_1.3.1          reticulate_1.13    
# [19] digest_0.6.25       RColorBrewer_1.1-2  SDMTools_1.1-221.1 
# [22] colorspace_1.4-1    cowplot_1.0.0       htmltools_0.4.0    
# [25] Matrix_1.2-18       R.oo_1.23.0         plyr_1.8.4         
# [28] pkgconfig_2.0.3     bibtex_0.4.2        tsne_0.1-3         
# [31] listenv_0.7.0       purrr_0.3.3         scales_1.1.0       
# [34] RANN_2.6.1          gdata_2.18.0        Rtsne_0.15         
# [37] tibble_2.1.3        ggplot2_3.2.1       ROCR_1.0-7         
# [40] pbapply_1.4-2       lazyeval_0.2.2      survival_3.1-8     
# [43] magrittr_1.5        crayon_1.3.4        R.methodsS3_1.7.1  
# [46] future_1.15.1       nlme_3.1-142        MASS_7.3-51.4      
# [49] gplots_3.0.1.1      ica_1.0-2           tools_3.6.3        
# [52] fitdistrplus_1.0-14 data.table_1.12.6   gbRd_0.4-11        
# [55] lifecycle_0.1.0     stringr_1.4.0       plotly_4.9.1       
# [58] munsell_0.5.0       cluster_2.1.0       irlba_2.3.3        
# [61] compiler_3.6.3      rsvd_1.0.2          caTools_1.17.1.2   
# [64] rlang_0.4.6         grid_3.6.3          ggridges_0.5.1     
# [67] rstudioapi_0.10     RcppAnnoy_0.0.14    htmlwidgets_1.5.1  
# [70] igraph_1.2.4.1      bitops_1.0-6        Seurat_3.1.1       
# [73] npsurv_0.4-0        gtable_0.3.0        codetools_0.2-16   
# [76] reshape2_1.4.3      R6_2.4.1            gridExtra_2.3      
# [79] zoo_1.8-6           dplyr_0.8.3         uwot_0.1.4         
# [82] future.apply_1.3.0  KernSmooth_2.23-16  metap_1.1          
# [85] ape_5.3             stringi_1.4.3       parallel_3.6.3     
# [88] Rcpp_1.0.3          sctransform_0.2.0   png_0.1-7          
# [91] vctrs_0.3.0         tidyselect_1.1.0    lmtest_0.9-37 


#============================#
#       Musician: Resonance  #
#           Date: 2020/07/12 #
# Revised author: Resonance  #
#     1st revise: 2020/07/12 #
#     2nd revise: 2020/07/13 #
#============================#
