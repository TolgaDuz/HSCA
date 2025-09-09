set.seed(42)

setwd("~/HSCA_data/extended/processed_seurat_objects/")

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
source("~/HSCA_code/Build_HSCA/HSCA_extended/helpers.R")

# ======================
# Load raw HSCA extended
# ======================

hsca_extended <- readRDS("./HSCA_extended_raw.rds")

DimPlot(hsca_extended, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(hsca_extended, group.by = "Dataset",  label = T, raster = F)
DimPlot(hsca_extended, group.by = "orig_celltype_lvl_3",  label = T, raster = F, repel = T)
DimPlot(hsca_extended, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)

FeaturePlot(hsca_extended, "S100A8", raster = F)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"), raster = F)


# ==== Marker gene discovery ====
hsca_extended <- JoinLayers(hsca_extended)

# Identify markers (only positive markers, limited cells per ident)
markers <- FindMarkers(
  object = hsca_extended,
  ident.1 = 80,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

head(markers, 20)
FeaturePlot(hsca_extended, rownames(markers)[1:12], raster = FALSE)

# ===============================================
# Subclustering: muscle cells
# ===============================================

# UMAP overview to find muscle cells
DimPlot(hsca_extended, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(hsca_extended, group.by = "sample", label = TRUE, raster = FALSE)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"))

# Low quality removal
hsca_extended <- subset(hsca_extended, subset = seurat_clusters == 26, invert = T)
# Remove samples with strong wound healing signatures
hsca_extended <- subset(hsca_extended, subset = sample %in% c("GSM8238438", "GSM8238437"), invert = T)
hsca_extended

# subset muscle cells
extended_muscle <- subset(hsca_extended, subset = seurat_clusters %in% c(84,36,44,13,7,37,103,109))

extended_muscle

DimPlot(extended_muscle, group.by = "seurat_clusters", label = T,raster = F)
FeaturePlot(extended_muscle, "ACTA", raster = F)

extended_muscle$Dataset %>% table()

# Remove datasets with very few cells
extended_muscle <- subset(
  extended_muscle,
  subset = Dataset %in% c(
    "Chu_Midha_2022", "Wang_Atwood_2020",
    "Cheng_Cho_2018_1", "Kim_Krüger_2022",
    "Cheng_Cho_2018_2", "Sun_Liu_2022"
  ),
  invert = T
)

extended_muscle$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_muscle[["RNA"]] <- split(extended_muscle[["RNA"]], f = extended_muscle$Dataset)

extended_muscle

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_muscle_scvi <- process_with_scvi(extended_muscle)
extended_muscle_scvi <- FindClusters(extended_muscle_scvi, resolution = 1.3)

DimPlot(extended_muscle_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_muscle_scvi, group.by = "seurat_clusters_full", label = T)
DimPlot(extended_muscle_scvi, group.by = "Dataset", label = T)
DimPlot(extended_muscle_scvi, group.by = "inherited_celltype_lvl_4", label = T, repel = T)

extended_muscle_scvi$seurat_clusters %>% table()

FeaturePlot(extended_muscle_scvi, c("percent_mito", "nFeature_RNA"))
FeaturePlot(extended_muscle_scvi, c("inherited_celltype_lvl_5_transfer_uncert"), raster = F)
extended_muscle_scvi

extended_muscle_scvi <- JoinLayers(extended_muscle_scvi)

extended_muscle_scvi

# ==== Marker gene discovery  ====
DefaultAssay(extended_muscle_scvi)
Idents(extended_muscle_scvi) <- extended_muscle_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_muscle_scvi,
  ident.1 = 64,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(extended_muscle_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(extended_muscle_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_muscle_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(extended_muscle_scvi, c("percent_mito", "nFeature_RNA"))

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

extended_muscle_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_muscle_scvi$seurat_clusters == 41 , "Skeletal muscle_3", extended_muscle_scvi$seurat_clusters)
extended_muscle_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_muscle_scvi$seurat_clusters == 34, "Muscle progenitor", extended_muscle_scvi$inherited_celltype_lvl_4_extended)
extended_muscle_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(8,36,15), "DES+ SMC", extended_muscle_scvi$inherited_celltype_lvl_4_extended)
extended_muscle_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(26,23,27,32,43,12,37,31,0,18,4,1,3,11,35), "RERGL+ SMC", extended_muscle_scvi$inherited_celltype_lvl_4_extended)
extended_muscle_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(33,28,29,7,6,16,25,21,9,30,22,38,40,17,2,13,24,5,10,20), "STEAP4+ SMC", extended_muscle_scvi$inherited_celltype_lvl_4_extended)
extended_muscle_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(14,19,39,42,20,40), "LQ", extended_muscle_scvi$inherited_celltype_lvl_4_extended)
extended_muscle_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(10,5,13,24,7,33,28), "CYP26B1+ SMC", extended_muscle_scvi$seurat_clusters)
extended_muscle_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(29,6,16), "FRMD3+ SMC", extended_muscle_scvi$inherited_celltype_lvl_5_extended)
extended_muscle_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_muscle_scvi$seurat_clusters %in% c(2,17,38,22,30,9,21,25), "TM4SF1+ SMC", extended_muscle_scvi$inherited_celltype_lvl_5_extended)

DimPlot(extended_muscle_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_muscle_scvi, group.by = "seurat_clusters", label = T)

saveRDS(extended_muscle_scvi, "./Extended_subsets/raw/muscle_extended.rds")
