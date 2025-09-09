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
# Subclustering: Eccrine sweat glands
# ===============================================

DefaultAssay(core_scvi) <- "RNA"

# UMAP overview to find eccrine sweat gland cells
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(core_scvi, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

# subset eccrine sweat gland cells
core_esg <- subset(core_scvi, subset = seurat_clusters %in% c(67,74,7,61,51,38,30))

# Visualization of ESG subset
DimPlot(core_esg, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

core_esg

core_esg$Dataset %>% table()

# Remove specific datasets with very few cells
core_esg_sub <- subset(
  core_esg,
  subset = Dataset == "Ganier_Lynch_2024_2",
  invert = T
)

core_esg_sub$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_esg_sub[["RNA"]] <- split(core_esg_sub[["RNA"]], f = core_esg_sub$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_esg_scvi <- process_with_scvi(core_esg_sub)
core_esg_scvi <- FindClusters(core_esg_scvi, resolution = 1.6)

DimPlot(core_esg_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_esg_scvi, group.by = "Dataset", label = TRUE)

FeaturePlot(core_esg_scvi, c("percent_mito", "nFeature_RNA"))

core_esg_scvi <- JoinLayers(core_esg_scvi)

# ==== Marker gene discovery  ====
markers <- FindMarkers(core_esg_scvi,
                       ident.1 = 23,
                       max.cells.per.ident = 500,
                       only.pos = TRUE)
head(markers, 20)

FeaturePlot(core_esg_scvi, rownames(markers)[1:10])
FeaturePlot(core_esg_scvi, rownames(markers)[11:20])
FeaturePlot(core_esg_scvi, rownames(markers)[21:30])

# ==== Feature Plots and Cell type annotation ====

# ESG consists of coil and duct parts
# Coil: dark, clear, and myoepithelial cells
# Duct: luminal and basal cells

# Coil cells
FeaturePlot(core_esg_scvi, c("DCD", "SCGB2A1", "SCGB2A2", "LIPH"))   # Dark & clear cells
FeaturePlot(core_esg_scvi, c("LCN2", "WFDC2", "ANKRD36C"))            # Clear cell 1
FeaturePlot(core_esg_scvi, c("LCN2", "CHRM3", "NRG3", "CLDN10"))      # Clear cell 2
FeaturePlot(core_esg_scvi, c("S100A1", "KRT19", "KRT8", "KRT18"))     # Secretory coil cells

# Duct cells
FeaturePlot(core_esg_scvi, c("S100A2", "KRT5", "KRT14", "CCL2", "CCR4"))   # Basal luminal
FeaturePlot(core_esg_scvi, c("KRT6A", "KRT16", "KRT77"))                  # Luminal C1, C3
FeaturePlot(core_esg_scvi, c("IVL", "SPINK5", "IFI27"))                   # C1 marker
FeaturePlot(core_esg_scvi, c("S100P"))                                    # C3 marker

# Myoepithelial
FeaturePlot(core_esg_scvi, c("MYH11", "ACTG2", "TAGLN", "ACTA2"))
FeaturePlot(object, c("MYH11", "ACTG2", "TAGLN", "ACTA2", "TPM2"))

DimPlot(core_esg_scvi, group.by = "seurat_clusters", label = T)

core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(13,6,10,2,16), "Dark cell", core_esg_scvi$seurat_clusters)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(11), "Myoepithelial", core_esg_scvi$inherited_celltype_lvl_4)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(31,19,23), "LQ", core_esg_scvi$inherited_celltype_lvl_4)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(18,14), "S100P+ luminal", core_esg_scvi$inherited_celltype_lvl_4)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(29), "S100P- luminal", core_esg_scvi$inherited_celltype_lvl_4)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(8,22,25,0,9,28,3,20,1,4,7,17,30), "Basal luminal", core_esg_scvi$inherited_celltype_lvl_4)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(5,21), "Clear cell 1", core_esg_scvi$inherited_celltype_lvl_4)
core_esg_scvi$inherited_celltype_lvl_4 <- ifelse(core_esg_scvi$seurat_clusters %in% c(26,24,27,12,15), "Clear cell 2", core_esg_scvi$inherited_celltype_lvl_4)

# Plot final annotations
DimPlot(core_esg_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_esg_scvi, group.by = "Dataset", label = TRUE)

saveRDS(core_esg_scvi,"./Core_subsets/raw/esg_core.rds")
