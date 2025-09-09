# setting working dir
setwd("~/HSCA_data/core/processed_seurat_objects/")

set.seed(42)

# load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(BPCells)
library(harmony)
library(reticulate)

# Configure reticulate to use scvi conda environment
reticulate::use_condaenv(
  "scvi",
  conda = "/opt/conda/condabin/conda",
  required = TRUE
)

# Increase allowed memory for large datasets
options(future.globals.maxSize = 3e+09)

# Load raw merged Seurat object
core_raw <- readRDS("./HSCA_core_raw.rds")

core_raw

# check metadata
core_raw$Dataset %>% table()
core_raw$Accession_source %>% table()
core_raw$orig_celltype_lvl_3 %>% table()

core_raw$genes <- rownames(core_raw)


# ==== Standard Seurat preprocessing ====

core_raw <- NormalizeData(core_raw)
core_raw[["RNA"]] <- split(core_raw[["RNA"]], f = core_raw$Dataset)
core_raw <- FindVariableFeatures(core_raw)
core_raw <- ScaleData(core_raw)
core_raw <- RunPCA(core_raw, assay = "RNA", reduction.name = "pca", npcs = 30)

core_raw

# ======================
# Integration with scVI
# ======================

set.seed(42)
core_scvi <- IntegrateLayers(
  object = core_raw,
  method = scVIIntegration,
  assay = "RNA",
  orig.reduction = "pca",
  new.reduction = "pca.svi",
  conda_env = "/opt/conda/envs/scvi",
  verbose = TRUE,
  seed = 42
)

# ==== Dimensionality reduction and clustering ====

# UMAP embedding on scVI-reduced dimensions
set.seed(42)
core_scvi <- RunUMAP(
  core_scvi,
  reduction = "pca.svi",
  reduction.name = "umap",
  return.model = TRUE,
  dims = 1:30
)

# Construct neighborhood graph
core_scvi <- FindNeighbors(
  core_scvi,
  assay = "sketch",
  dims = 1:30,
  reduction = "pca.svi"
)

# Cluster cells
core_scvi <- FindClusters(
  core_scvi,
  resolution = 1.6
)

# Visualiztation
DimPlot(core_scvi, label = TRUE)
DimPlot(core_scvi, group.by = "seurat_clusters")
DimPlot(core_scvi, group.by = "Dataset")
DimPlot(core_scvi, group.by = "orig_celltype_lvl_3", label = TRUE, raster = FALSE)
DimPlot(core_scvi, group.by = "Accession_source")

# Feature plots
FeaturePlot(core_scvi, c("CCER2", "FLG"), raster = FALSE)
FeaturePlot(seurat_object, feature = "nFeature_Spatial")

saveRDS(
  object = core_scvi,
  file = "./HSCA_core_scvi.rds"
)
