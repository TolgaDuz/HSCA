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

# Visualization
DimPlot(core_scvi, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(core_scvi, group.by = "orig_celltype_lvl_3", label = TRUE, raster = FALSE, repel = FALSE)

# Backup original cluster assignments
core_scvi$seurat_clusters_full <- core_scvi$seurat_clusters

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
# Subclustering: "Other" cell lineages
# ===============================================

# subset
core_other <- subset(
  core_scvi,
  subset = seurat_clusters %in% c(76,80,75,55,78,12),
  invert = F
)

core_other

core_other$orig_celltype_lvl_3 %>% table()
core_other$seurat_clusters %>% table()

DimPlot(core_other, group.by = "seurat_clusters", label = T)
DimPlot(core_other, group.by = "orig_celltype_lvl_3", label = T)

core_other

core_other$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_other[["RNA"]] <- split(core_other[["RNA"]], f = core_other$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_other_scvi <- process_with_scvi(core_other)

core_other_scvi <- FindClusters(core_other_scvi, resolution = 2.2)

# Inspect
DimPlot(core_other_scvi, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(core_other_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(core_other_scvi, group.by = "Dataset", label = TRUE, raster = FALSE)
DimPlot(core_other_scvi, split.by = "anatomical_region_level2", label = TRUE, raster = FALSE)

FeaturePlot(core_other_scvi, c("percent_mito", "nFeature_RNA"), raster = F)

core_other_scvi <- JoinLayers(core_other_scvi)


# ==== Marker gene discovery  ====
Idents(core_other_scvi)
markers <- FindMarkers(core_other_scvi, ident.1 = 24, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 20)

# FeaturePlots for top markers
FeaturePlot(core_other_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(core_other_scvi, rownames(markers)[11:20], raster = FALSE)
FeaturePlot(core_other_scvi, rownames(markers)[21:30], raster = FALSE)

# Melanocytes
FeaturePlot(core_other_scvi, features = c("MLANA", "QPCT", "TYRP1", "PMEL"))
# Merkel cells
FeaturePlot(core_other_scvi, features = c("CCER2", "CCK", "KCNMB2", "MIAT", "KRT8", "KRT18", "KRT20"))
# Erythrocyte
FeaturePlot(core_other_scvi, features = c("ALAS2", "HBA1", "HBB", "HBA2"))
# Adipocyte
FeaturePlot(core_other_scvi, features = c("FABP4", "ADIPOQ", "RBP4", "PLIN4"))
# Schwann cells
FeaturePlot(core_other_scvi, features = c("PRX", "GLDN", "MPZ", "MLIP"))
# Neuron
FeaturePlot(core_other_scvi, features = c("NRXN1", "STARD13", "XKR4", "NTM"))
# Skeletal muscle
FeaturePlot(core_other_scvi, features = c("MYL1", "ACTA1", "CKM", "COX6A2"))
# CNTNAP2+ neuron
FeaturePlot(core_other_scvi, features = c("CNTNAP2", "CSMD1", "CNTN5", "PTPRD"))

DimPlot(core_other_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = TRUE)

core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(28), "Adipocyte_2", core_other_scvi$seurat_clusters)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(16,0,9,2,3,6,10,18,5,23,31,29,32,13,8,22,11,7,21,14,1), "Melanocyte_2", core_other_scvi$inherited_celltype_lvl_4)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(33), "Merkel cell_3", core_other_scvi$inherited_celltype_lvl_4)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(12,17,4), "Neuron_2", core_other_scvi$inherited_celltype_lvl_4)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(25), "Skeletal muscle_3", core_other_scvi$inherited_celltype_lvl_4)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(19), "Schwann cell_2", core_other_scvi$inherited_celltype_lvl_4)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(24), "CNTNAP2+ Neuron", core_other_scvi$inherited_celltype_lvl_4)
core_other_scvi$inherited_celltype_lvl_4 <- ifelse(core_other_scvi$seurat_clusters %in% c(26,20,15,30,27,34), "LQ", core_other_scvi$inherited_celltype_lvl_4)

# Plot annotations
DimPlot(core_other_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE, raster = FALSE)
DimPlot(core_other_scvi, group.by = "celltype_lvl_4", label = TRUE, raster = FALSE)

saveRDS(core_other_scvi, "./Core_subsets/raw/ohter_cells_core.rds")
