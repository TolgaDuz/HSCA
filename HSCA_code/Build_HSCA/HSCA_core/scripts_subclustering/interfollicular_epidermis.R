set.seed(42)

setwd("~/HSCA_data/core/processed_seurat_objects/")

# Load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(BPCells)
library(reticulate)

# Configure reticulate to use scvi conda environment
reticulate::use_condaenv(
  "scvi",
  conda = "/opt/conda/condabin/conda",
  required = TRUE
)

# Import helper functions
source("~/HSCA_code/Build_HSCA/HSCA_core/helpers.R")

# ======================
# Load raw HSCA core object
# ======================

core_scvi <- readRDS("./HSCA_core_scvi.rds")

# UMAP visualization by different metadata
DimPlot(core_scvi, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(core_scvi, group.by = "orig_celltype_lvl_3", label = TRUE, raster = FALSE, repel = FALSE)

# Backup original cluster assignments
core_scvi$seurat_clusters_full <- core_scvi$seurat_clusters

# Feature visualization
FeaturePlot(core_scvi, c("percent_mito", "nFeature_RNA"), raster = FALSE)

# ==== Marker gene discovery ====
core_scvi <- JoinLayers(core_scvi)

# Identify markers (only positive markers, limited cells per ident)
markers <- FindMarkers(
  object = core_scvi,
  ident.1 = 80,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

head(markers, 20)
FeaturePlot(core_scvi, rownames(markers)[1:12], raster = FALSE)

# ===============================================
# Subclustering: Interfollicular epidermis (IFE)
# ===============================================

core_scvi

# UMAP overview to find IFE cells
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(core_scvi, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

core_scvi$seurat_clusters_full <- core_scvi$seurat_clusters
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters_full", label = T, raster = F)

# subset to IFE
core_ife <- subset(core_scvi, subset = seurat_clusters %in% c(
  49,59,44,27,32,22,31,8,60,63,43,3,54,26,19,50,45,28,34,39,5,14,57,36,11,35,
  21,18,13,20,68,56
))

# Quick checks
core_ife
core_scvi
core_ife$orig_celltype_lvl_3 %>% table()
core_scvi$orig_celltype_lvl_3 %>% table()
core_ife$seurat_clusters %>% table()

DimPlot(core_ife, reduction = "umap", group.by = "seurat_clusters", label = T)

# Cluster refinement
core_ife <- FindClusters(core_ife, resolution = 3)
DimPlot(core_ife, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(core_ife, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE)

core_ife

core_ife$Dataset %>% table()
# Remove dataset with very few cells
core_ife <- subset(
  core_ife,
  subset = Dataset == "Ganier_Lynch_2024_2",
  invert = T
)
core_ife$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_ife[["RNA"]] <- split(core_ife[["RNA"]], f = core_ife$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_ife_scvi <- process_with_scvi(core_ife)
core_ife_scvi <- FindClusters(core_ife_scvi, resolution = 3.7)

# Visualize clusters and metadata
DimPlot(core_ife_scvi, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(core_ife_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(core_ife_scvi, group.by = "Dataset", label = TRUE, raster = FALSE)
DimPlot(core_ife_scvi, group.by = "anatomical_region_level2", label = TRUE, raster = FALSE)

FeaturePlot(core_ife_scvi, c("percent_mito", "nFeature_RNA"), raster = F)

# Feature plots
FeaturePlot(core_ife_scvi, c("percent_mito", "nFeature_RNA"), raster = FALSE)
FeaturePlot(core_ife_scvi, c("MKI67", "LCE1A", "SPRR1B", "FLG"), raster = FALSE)

core_ife_scvi <- JoinLayers(core_ife_scvi)

# ==== Feature Plots and Cell type annotation ====

FeaturePlot(core_ife_scvi, c("SERPINB4", "KRT1", "PTN", "C1QTNF12"))
FeaturePlot(core_ife_scvi, c("KRT15", "KRT1", "MKI67", "FLG"))
FeaturePlot(core_ife_scvi, c("COL17A1", "KRT14", "POSTN", "KRT5"))
FeaturePlot(core_ife_scvi, c("DIO2", "WIF1", "FRZB", "SHROOM3"))
FeaturePlot(core_ife_scvi, c("S100A8", "S100A7", "S100A9", "KRT6A"))
FeaturePlot(core_ife_scvi, c("PTN", "GATA6", "KRT79", "LGR6"))
FeaturePlot(core_ife_scvi, c("C1QTNF12", "GATA6", "KRT79", "PTN"))
FeaturePlot(core_ife_scvi, c("CALB2", "PI3", "KRT79", "IVL"))
FeaturePlot(core_ife_scvi, c("KLK6", "LCN2", "RHCG", "CDA"))
FeaturePlot(core_ife_scvi, c("CRAT", "SEC14L6", "PLIN5", "PNPLA3"))
FeaturePlot(core_ife_scvi, c("NNAT", "IL1R2", "KRT7", "WFDC2"))
FeaturePlot(core_ife_scvi, c("MKI67", "TOP2A"))
FeaturePlot(core_ife_scvi, c("LGR5", "RARB", "TBX2", "LHX2"))
FeaturePlot(core_ife_scvi, c("RXRA", "RARG"))
FeaturePlot(core_ife_scvi, "KRT1")

# DimPlots with/without repel
DimPlot(core_ife_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = TRUE)
DimPlot(core_ife_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = FALSE)

core_ife_scvi$inherited_celltype_lvl_4 <- ifelse(
  core_ife_scv$seurat_clusters %in% c(
    14,70,7,51,68,71,59,40,48,22,28,49,89,74,0,75,82,80,
    53,43,63,5,31,47,57,50
  ),
  "Basal KC_3",
  core_ife_scv$seurat_clusters
)
core_ife_scvi$inherited_celltype_lvl_4 <- ifelse(
  core_ife_scv$seurat_clusters %in% c(
    61,25,23,78,19,36,30,13,46,39,44,57,10,66,67,60,2,42,4,8,3,29,6,90,16,18,
    26,9,79,58,88,54,52,81,33,72,20,35,37,65,85,56,27,24,17,12,1,91
  ),
  "Spinous KC_3",
  core_ife_scv$inherited_celltype_lvl_4
)
core_ife_scvi$inherited_celltype_lvl_4 <- ifelse(core_ife_scv$seurat_clusters %in% c(11,84,45,41,69,83), "Prolif. KC_3", core_ife_scv$inherited_celltype_lvl_4)
core_ife_scvi$inherited_celltype_lvl_4 <- ifelse(core_ife_scv$seurat_clusters %in% c(76,86,64,73), "Granular KC_3", core_ife_scv$inherited_celltype_lvl_4)
core_ife_scvi$inherited_celltype_lvl_4 <- ifelse(core_ife_scv$seurat_clusters %in% c(87), "Cornified KC_3", core_ife_scv$inherited_celltype_lvl_4)
core_ife_scvi$inherited_celltype_lvl_4 <- ifelse(core_ife_scv$seurat_clusters %in% c(77,21,32,38,15,55, 62,34), "LQ", core_ife_scv$inherited_celltype_lvl_4)

# Final DimPlots
DimPlot(core_ife_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE, raster = FALSE)
DimPlot(core_ife_scvi, group.by = "celltype_lvl_4", label = TRUE, raster = FALSE)

# Save annotated object
saveRDS(core_ife_scvi, "./Core_subsets/raw/ife_core.rds")
