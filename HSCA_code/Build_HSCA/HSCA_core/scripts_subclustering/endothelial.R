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
FeaturePlot(core_scvi, "PAX7", raster = FALSE)


# ===============================================
# Subclustering: Endothelial cells
# ===============================================

DefaultAssay(core_scvi) <- "RNA"

# UMAP overview to find endothelial cells
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(core_scvi, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

# Subset endothelial clusters
core_ec <- subset(core_scvi, subset = seurat_clusters %in% c(0,64,29,6,52,37))

DimPlot(core_ec, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
core_ec$Dataset %>% table()

# Exclude selected datasets with very few cells
core_ec_sub <- subset(
  core_ec,
  subset = Dataset %in% c("Cheng_Cho_2018_1", "Takahashi_Lowry_2019"),
  invert = TRUE
)

core_ec_sub$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_ec_sub[["RNA"]] <- split(core_ec_sub[["RNA"]], f = core_ec_sub$Dataset)

core_ec_sub

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_ec_scvi <- process_with_scvi(core_ec_sub)

set.seed(42)
core_ec_scvi <- FindClusters(core_ec_scvi, resolution = 1)

core_ec_scvi

# Inspect clustering
DimPlot(core_ec_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_ec_scvi, group.by = "Dataset", label = TRUE)
DimPlot(core_ec_scvi, group.by = "celltype_lvl_4", label = TRUE)
FeaturePlot(core_ec_scvi, c("percent_mito", "nFeature_RNA"))

core_ec_scvi <- JoinLayers(core_ec_scvi)

# ==== Marker gene discovery  ====
markers <- FindMarkers(core_ec_scvi, ident.1 = 12, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 20)
FeaturePlot(core_ec_scvi, rownames(markers)[1:12])

# ==== Feature Plots and Cell type annotation ====

#### VEC ####
# Arterial EC
FeaturePlot(core_ec_sub.scvi, c("IGFBP3", "ARL15", "CXCL12", "HEY1"))
FeaturePlot(core_ec_sub.scvi, c("IGFBP3", "PCSK5", "NEBL", "SEM3AG", "FBLN5"))

# Capillary EC
FeaturePlot(core_ec_sub.scvi, c("APLN", "VWA1", "RGCC", "H19"))
FeaturePlot(core_ec_sub.scvi, c("FABP4", "CD36", "BTNL9", "RBP7"))

# Venule subsets
# Venous 1
FeaturePlot(core_ec_sub.scvi, c("G0S2", "SELE", "ACKR1", "CSF3", "STC1"))
FeaturePlot(core_ec_sub.scvi, c("CSF3", "VCAM1", "SOD2", "IL6"))
# Venous 2
FeaturePlot(core_ec_sub.scvi, c("CCL14", "AQP1", "ID1", "SOX18"))
FeaturePlot(core_ec_sub.scvi, c("NAV3", "LDB2", "TLL1", "CCSER1"))

# RGS5+ cells
FeaturePlot(core_ec_sub.scvi, c("RGS5", "PDFGRB", "ANPEP"))

### Lymphatic EC (LEC) ###
FeaturePlot(core_ec_sub.scvi, c("SCG3", "NRXN3", "LYPD6", "RADIL", "HGF", "ADM"))
FeaturePlot(core_ec_sub.scvi, c("NEO1", "SCN3A", "DPP4", "SCG3"))
FeaturePlot(core_ec_sub.scvi, c("CCL21", "LYVE1", "FXYD6", "NRP2", "PTX3"))

### Low-quality (LQ) cells ###
FeaturePlot(core_ec_sub.scvi, c("PTPRC", "SAMSN1", "CD53", "IKZF1", "RHOH", "BCL11B"))

core_ec_scvi$orig_celltype_lvl_3 %>% table()

core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(7) , "RSG5+ EC", core_ec_scvi$seurat_clusters)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(10,2,8,4,1), "Venous 1 EC", core_ec_scvi$inherited_celltype_lvl_4)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(0,11,16), "Venous 2 EC", core_ec_scvi$inherited_celltype_lvl_4)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(17), "SCG3+ LEC", core_ec_scvi$inherited_celltype_lvl_4)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(5), "LYVE1+ LEC", core_ec_scvi$inherited_celltype_lvl_4)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(9,6,15), "Capillary EC_3", core_ec_scvi$inherited_celltype_lvl_4)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(3,14), "Arterial EC_3", core_ec_scvi$inherited_celltype_lvl_4)
core_ec_scvi$inherited_celltype_lvl_4 <- ifelse(core_ec_scvi$seurat_clusters %in% c(13,12), "LQ", core_ec_scvi$inherited_celltype_lvl_4)

# Final visualization
DimPlot(core_ec_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_ec_scvi, group.by = "Dataset", label = TRUE)

# Save annotated object
saveRDS(
  core_ec_scvi,
  "./Core_subsets/raw/ec_core.rds"
)