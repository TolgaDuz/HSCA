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
# Subclustering: Muscle cells
# ===============================================

core_scvi

# UMAP overview to find muscle cells
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(core_scvi, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

# subset muscle cells
core_muscle <- subset(object, subset = seurat_clusters %in% c(79,62,1,15))

core_muscle
object

# Visualize muscle cluster
DimPlot(core_muscle, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(core_muscle, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE)
core_muscle$Dataset %>% table()

# Exclude selected datasets with very few cells
core_muscle_sub <- subset(
  core_muscle,
  subset = Dataset %in% c("Cheng_Cho_2018_1", "Takahashi_Lowry_2019"),
  invert = TRUE
)

core_muscle_sub$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_muscle_sub[["RNA"]] <- split(core_muscle_sub[["RNA"]], f = core_muscle_sub$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_muscle_scvi <- process_with_scvi(core_muscle_sub)
core_muscle_scvi <- FindClusters(core_muscle_scvi, resolution = 1.9)

# Cluster inspection
DimPlot(core_muscle_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_muscle_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(core_muscle_scvi, group.by = "Dataset", label = TRUE)
DimPlot(core_muscle_scvi, group.by = "celltype_lvl_4", label = TRUE)

core_muscle_scvi <- JoinLayers(core_muscle_scvi)

# ==== Marker gene discovery  ====
markers <- FindMarkers(
  object = core_muscle_scvi,
  ident.1 = 8,
  max.cells.per.ident = 500,
  only.pos = TRUE
)
head(markers, 20)

# Marker visualization
FeaturePlot(core_muscle_scvi, rownames(markers)[1:12])
FeaturePlot(core_muscle_scvi, rownames(markers)[13:24])

# ==== Feature Plots and Cell type annotation ====

# Muscle stem cell
FeaturePlot(core_muscle_scvi, c("MYF5", "RBP1", "CHRNA1", "DLK1", "CADM2", "PAX7"))

# Skeletal muscle
FeaturePlot(core_muscle_scvi, c("COX6A2", "KLHL41", "MYL1", "CSRP3", "SPMX", "STAC3", "MB", "TCAP", "CKM", "MYLPF"))

# DES+ SMC
FeaturePlot(core_muscle_scvi, c("DES", "PCP4", "ACTG2", "PDE4D", "SLC8A1", "SORBS1", "PRUNE2", "PCDH7", "ROBO2", "ITGA8"))

# HMCN2 + SMC
FeaturePlot(core_muscle_scvi, c("CYFIP2", "ALDH1B1", "MYH10", "HMCN2"))

# STEAP4+ SMC
FeaturePlot(core_muscle_scvi, c("STEAP4", "COL6A3"))

# FRMD3+
FeaturePlot(core_muscle_scvi, c("FRMD3", "PDZD2", "TSHZ2", "CD36"))

# TM4SF1+ SMC
FeaturePlot(core_muscle_scvi, c("TM4SF1", "LMCD1", "C2orf40", "GGT5", "CFD"))

# BASP1+ SMC
FeaturePlot(core_muscle_scvi, c("CYP26B1", "CFH", "P2RY14", "BASP1", "ADAMTS5", "APCDD1", "COL23A1", "GPM6B")) #

# RERGL+ SMC
FeaturePlot(core_muscle_scvi, c("RERGL", "SORBS2", "PLN", "SBCG", "MYH11", "BCAM", "PGAM2", "NRGN"))

# General smooth muscle markers
FeaturePlot(core_muscle_scvi, c("ACTA2", "TAGLN"))
FeaturePlot(core_muscle_scvi, c("KRT2", "CALML5"))

# core_muscle_scvi$inherited_celltype_lvl_4 <- ifelse(core_muscle_scvi@meta.data$seurat_clusters == 24 , "Skeletal Muscle_3", core_muscle_scvi@meta.data$seurat_clusters)
core_muscle_scvi$inherited_celltype_lvl_4 <- ifelse(core_muscle_scvi$seurat_clusters == 27, "Muscle progenitor", core_muscle_scvi$seurat_clusters)
core_muscle_scvi$inherited_celltype_lvl_4 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(20,13), "DES+ SMC", core_muscle_scvi$inherited_celltype_lvl_4)
core_muscle_scvi$inherited_celltype_lvl_4 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(2,25,14,0,15,17,26,4), "RERGL+ SMC", core_muscle_scvi$inherited_celltype_lvl_4)
core_muscle_scvi$inherited_celltype_lvl_4 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(12,21,5,16,9,7,8,10,11,1,19,6,24,3), "STEAP4+ SMC", core_muscle_scvi$inherited_celltype_lvl_4)
core_muscle_scvi$inherited_celltype_lvl_4 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(22,23,18), "LQ", core_muscle_scvi$inherited_celltype_lvl_4)
# Subtype of STEAP4+ cells, hence assgined to lvl_5
core_muscle_scvi$inherited_celltype_lvl_5 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(12,21,5,16,19,9,11), "CYP26B1+ SMC", core_muscle_scvi$seurat_clusters)
core_muscle_scvi$inherited_celltype_lvl_5 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(10,8,7), "FRMD3+ SMC", core_muscle_scvi$inherited_celltype_lvl_5)
core_muscle_scvi$inherited_celltype_lvl_5 <- ifelse(core_muscle_scvi$seurat_clusters %in% c(6,1,24,3), "TM4SF1+ SMC", core_muscle_scvi$inherited_celltype_lvl_5)

# Plot annotation
DimPlot(core_muscle_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_muscle_scvi, group.by = "inherited_celltype_lvl_5", label = TRUE)

saveRDS(core_muscle_scvi,"./Core_subsets/raw/muscle_core.rds")
