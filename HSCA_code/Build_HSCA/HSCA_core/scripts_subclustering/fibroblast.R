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
# Subclustering: Fibroblasts
# ===============================================


core_scvi

# UMAP overview to find fibroblasts
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(core_scvi, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

# subset fibroblasts
core_fb <- subset(
  core_scvi,
  subset = seurat_clusters %in% c(9,46,65,4,53,58,23,66,47,24,70)
)

DimPlot(core_fb, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
core_fb

core_fb <- JoinLayers(core_fb)

core_fb$Dataset %>% table()

# Remove specific datasets with very few cells
core_fb_sub <- subset(
  core_fb,
  subset = Dataset %in% c("Cheng_Cho_2018_1", "Takahashi_Lowry_2019"),
  invert = TRUE
)

core_fb_sub$Dataset %>% table()

DefaultAssay(core_fb_sub) <- "RNA"
# Split RNA assay by dataset for scVI processing
core_fb_sub[["RNA"]] <- split(core_fb_sub[["RNA"]], f = core_fb_sub$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_fb_scvi <- process_with_scvi(core_fb_sub)

core_fb_scvi <- FindClusters(core_fb_scvi, resolution = 2.3)

# Cluster visualization
DimPlot(core_fb_scvi, label = TRUE, raster = FALSE)
DimPlot(core_fb_scvi, group.by = "Dataset", label = TRUE, raster = FALSE)

FeaturePlot(core_fb_scvi, c("percent_mito", "nFeature_RNA"))

core_fb_scvi <- JoinLayers(core_fb_scvi)
DimPlot(core_fb_scvi, group.by = "seurat_clusters", label = TRUE)

# ==== Marker gene discovery ====
markers <- FindMarkers(
  object = core_fb_scvi,
  ident.1 = 23,
  max.cells.per.ident = 500,
  only.pos = TRUE
)
head(markers, 20)

FeaturePlot(core_fb_scvi, c(rownames(markers)[1:12]))
FeaturePlot(core_fb_scvi, c(rownames(markers)[13:24]))

# Metadata plots
DimPlot(core_fb_scvi, group.by = "Dataset", label = TRUE)
DimPlot(core_fb_scvi, group.by = "anatomical_region_level3", label = TRUE)
DimPlot(core_fb_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(core_fb_scvi, group.by = "age_range", label = TRUE)
DimPlot(core_fb_scvi, group.by = "sample", label = TRUE)
DimPlot(core_fb_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)

object.fb

# ==== Feature Plots and Cell type annotation ====
FeaturePlot(core_fb_scvi, c("IFI27"))
FeaturePlot(core_fb_scvi, c("COL6A5", "COL18A1", "CCL2", "CCL19"))

# -------------------------
# Marker sets by Ascension
# Paper 1: Human Dermal Fibroblast Subpopulations Are
# Conserved across Single-Cell RNA Sequencing Studies.
# J Invest Dermatol, 141(7), 1735-1744 e1735. https://doi.org/10.1016/j.jid.2020.11.028
# Paper 2: Ascensión, A. M., & Izeta, A. (2024).
# A consensus single-cell transcriptomic atlas of dermal fibroblast heterogeneity.
# bioRxiv, 2024.2009. 2005.611379.
# -------------------------

# A
FeaturePlot(core_fb_scvi, features = c("SFRP2", "ELN", "MMP2", "QPCT"))
# A1
FeaturePlot(core_fb_scvi, features = c("IGFBP6", "PI116", "SLPI", "CCN5"))   # Paper 1
FeaturePlot(core_fb_scvi, features = c("WISP2", "SEMA3B", "LGR5", "ANGPTL5"))  # Paper 2
# A2
FeaturePlot(core_fb_scvi, features = c("APCDD1", "COL18A1", "COMP", "NKD2"))  # Paper 1
FeaturePlot(core_fb_scvi, features = c("HSPB3", "COL18A1", "COL6A5", "NKD2")) # Paper 2
# A3
FeaturePlot(core_fb_scvi, features = c("ELN", "RGCC", "SGCA", "WIF1"))        # Paper 1
FeaturePlot(core_fb_scvi, features = c("SOSTDC1", "CORIN", "SGCA", "WIF1"))  # Paper 2
# A4
FeaturePlot(core_fb_scvi, features = c("FBN1", "PCOLCE2", "PRG4", "SFRP4"))  # Paper 1
FeaturePlot(core_fb_scvi, features = c("C1QTNF3", "SCARA5", "PRG4", "TRAC")) # Paper 2

# B
FeaturePlot(core_fb_scvi, c("APOE", "C7", "CYGB", "IGFBP7"))
# B1
FeaturePlot(core_fb_scvi, c("CCL2", "ITM2A", "SPSB1", "TNFAIP6")) # Paper 1
FeaturePlot(core_fb_scvi, c("GEM", "CXCL2", "CXCL1", "TNFAIP6")) # Paper 2
# B2
FeaturePlot(core_fb_scvi, c("CCDC146", "CCL19", "CD74", "TNFSF13B")) # Paper 1
FeaturePlot(core_fb_scvi, c("GGT5", "IL33", "C7", "SCN4B")) # Paper 2
# B3
FeaturePlot(core_fb_scvi, c("CCL19", "CTSH", "RBP5", "ACHE")) # only in paper 2
# B4
FeaturePlot(core_fb_scvi, c("EFEMP1", "ITM2A", "MYOC", "GDF10")) # only in paper 2

# C
FeaturePlot(core_fb_scvi, c("SFRP1", "TNMD", "DKK3", "TNN"))
# C1 (DS)
FeaturePlot(core_fb_scvi, c("COL11A1", "DPEP1", "TNMD", "WFDC1")) # only in paper 1
FeaturePlot(core_fb_scvi, c("COL11A1", "MEF2C", "DPEP1", "WFDC1")) # only in paper 2
# C2 (Upper Bulge DP)
FeaturePlot(core_fb_scvi, c("COCH", "CRABP1", "FIBIN", "RSPO4")) # only in paper 1
FeaturePlot(core_fb_scvi, c("COCH", "CRABP1", "NDNF", "SLITRK6")) # only in paper 2
# C3
FeaturePlot(core_fb_scvi, c("ASPN", "F2R", "GPM6B", "POSTN")) # only in paper 1
FeaturePlot(core_fb_scvi, c("LTBP2", "LRRC15", "BGN", "POSTN")) # only in paper 2
# C4
FeaturePlot(core_fb_scvi, c("ANGPTL7", "APOD", "ECRG4", "TM4SF1")) # only in paper 1
FeaturePlot(core_fb_scvi, c("ANGPTL7", "APOD", "ECRG4", "TM4SF1")) # only in paper 2
# C5
FeaturePlot(core_fb_scvi, c("IGFBP3", "SLC5A3", "WNT5A", "LUZP2")) # only in paper 2
FeaturePlot(core_fb_scvi, features = c("PGM2L1", "PTCH1", "KIF26B", "GRIK1"))
FeaturePlot(core_fb_scvi, c("IGFBP3", "DKK2", "RSPO4","WIF1","EDN3","DIO3"))

# D
# D1
FeaturePlot(core_fb_scvi, c("ANGPTL7", "ENTPD2", "CDH19","ATP1A2"))
# D2
FeaturePlot(core_fb_scvi, c("BNC2", "ITGA6", "ITGB4","TNNC1"))
# E1
FeaturePlot(core_fb_scvi, features = c("IGFBP2", "COL26A1", "WNT2", "RAMP1"))

# Interesting
# COL8A1
FeaturePlot(core_fb_scvi, features = c("COL8A1", "PTX3"))
# PDZRN4
FeaturePlot(core_fb_scvi, features = c("DKK1","PDZRN4", "SGCZ"))

core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(4,23,38,14,39), "A1", core_fb_scvi$seurat_clusters)
core_fb_scvi$inherited_celltype_lvl_5 <- ifelse(core_fb_scvi$seurat_clusters %in% c(23), "PDZRN4+ FB", core_fb_scvi$seurat_clusters)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(28,45,12,18,0,41,32,25,33,7,26,21,40), "A2", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(11,2), "A3", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(19), "A4", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(15), "B1", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(17,13), "B2", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(3,10,16,8), "B3", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(44,24,20,29,36), "B4", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(34), "D1", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(30), "D2", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(42,9), "Dermal sheath (C1)", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(1,27), "C3", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(43,5), "Outer bulge DP (C2)", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(6), "RAMP1+ Fibro (E1)", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(31,35), "Anagen DP (C5)", core_fb_scvi$inherited_celltype_lvl_4)
core_fb_scvi$inherited_celltype_lvl_4 <- ifelse(core_fb_scvi$seurat_clusters %in% c(37,22), "LQ", core_fb_scvi$inherited_celltype_lvl_4)

DimPlot(core_fb_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_fb_scvi, group.by = "inherited_celltype_lvl_5", label = TRUE)

saveRDS(core_fb_scvi,"./Core_subsets/raw/fb_core.rds")
