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
# Subclustering: Immune cells
# ===============================================

DefaultAssay(core_scvi) <- "RNA"

# UMAP overview to find immune cells
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(core_scvi, reduction = "umap", group.by = "orig_celltype_lvl_3", label = TRUE, raster = FALSE)

# subset immune cells
core_immu <- subset(
  core_scvi,
  subset = seurat_clusters %in% c(77,48,72,73,25,41,2,42)
)

# Quick inspection
DimPlot(core_immu, group.by = "seurat_clusters", label = TRUE)

core_immu$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_immu[["RNA"]] <- split(core_immu[["RNA"]], f = core_immu$Dataset)

core_immu

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_immu_scvi <- process_with_scvi(core_immu)

set.seed(42)
core_immu_scvi <- FindClusters(core_immu_scvi, resolution = 5.9)
core_immu_scvi

# Inspect
DimPlot(core_immu_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_immu_scvi, group.by = "seurat_clusters_full", label = TRUE)
DimPlot(core_immu_scvi, group.by = "Dataset", label = TRUE)
DimPlot(core_immu_scvi, group.by = "orig_celltype_lvl_3", label = TRUE, repel = TRUE)

FeaturePlot(core_immu_scvi, c("percent_mito", "nFeature_RNA"))

core_immu_scvi <- JoinLayers(core_immu_scvi)

# ==== Marker gene discovery  ====
Idents(core_immu_scvi)
markers <- FindMarkers(
  object = core_immu_scvi,
  ident.1 = 52,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
core_immu_scvi$seurat_clusters %>% table()

# Feature plots for top marker genes
FeaturePlot(core_immu_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(core_immu_scvi, rownames(markers)[13:24], raster = FALSE)

DimPlot(core_immu_scvi, group.by = "seurat_clusters", label = T)

# ==== Feature Plots and Cell type annotation ====

# CD4+ T cells (including Th17)
FeaturePlot(core_immu_scvi, c(
  "CCR6", "NTRK2", "IL26", "TOX", "CAMK4", "BATF", "SPOCK2", "DUSP4"
))
FeaturePlot(core_immu_scvi, c("BACH2", "CNOT6L", "PBX4", "PGAP1"))

# CD8+ T cells (cytotoxic)
FeaturePlot(core_immu_scvi, c(
  "GZMK", "GZMA", "CD8A", "CRTAM", "CD8B", "ZNF683", "NKG7", "CCL5", "IFNG", "CST7", "GZMH"
))

# Regulatory T cells
FeaturePlot(core_immu_scvi, c("IKZF2", "TIGIT", "CTLA4", "CD27", "FOXP3", "IL2RA"))

# Naive B cell & pDC
FeaturePlot(core_immu_scvi, c("CD79A", "IGKC", "MS4A1", "IGHM", "PAX5"))
FeaturePlot(core_immu_scvi, c("MS4A1", "CD79A", "BANK1", "CD79B", "RALGPS2", "IGKC", "IGHM", "GZMB"))

# Langerhans cells (LC)
FeaturePlot(core_immu_scvi, c("CD207", "FCGBP", "CD1A", "HPGDS", "PLEK2"))

# Proliferating dendritic cells
FeaturePlot(core_immu_scvi, c("CD1C", "FCER1A", "CLEC10A", "MKI67", "PCLAF", "TOP2A", "CDK1"))

# Conventional DCs (cDC1 / cDC2)
FeaturePlot(core_immu_scvi, c("CD1C", "FCER1A", "CLEC10A", "IL1B", "CXCL8", "IL1R2", "LGALS2", "CD1E"))
FeaturePlot(core_immu_scvi, c("WDFY4", "DNASE1L3", "CADM1", "C1orf54", "CLNK", "CLEC9A"))

# Mature DC
FeaturePlot(core_immu_scvi, c("LAMP3", "CCR7", "FSCN1", "WNT5B", "MARCKSL1", "CERS6", "SLCO5A1", "NAV1", "CCL2"))
FeaturePlot(core_immu_scvi, c("CCL22", "POGLUT1", "LAMP3", "CCR7"))

# Macrophages
FeaturePlot(core_immu_scvi, c("SELENOP", "F13A1", "C1QA", "STAB1", "DAB2", "FOLR2")) # Anti-inflammatory
FeaturePlot(core_immu_scvi, c("CXCL2", "CXCL3", "RNASE1", "SELENOP", "CTSL", "C1QA")) # Inflammatory
FeaturePlot(core_immu_scvi, c("OLR1", "TREM2", "C3", "C1QB")) # TREM2+ C3+
FeaturePlot(core_immu_scvi, c("TREM2", "LPL", "SPP1", "CYP27A1")) # TREM2+ LPL+

# Monocytes & neutrophils
FeaturePlot(core_immu_scvi, c("S100A9", "FCN1", "AQP9", "VCAN")) # Monocytes
FeaturePlot(core_immu_scvi, c("S100A8", "S100A9", "FCGR3B", "CMTM2", "AQP9")) # Neutrophils

# Plasma cells
FeaturePlot(core_immu_scvi, c("JCHAIN", "LILRA4", "LRRC26", "LAMP5", "PTCRA", "SCT", "SMIM5", "GZMB"))

# Mast cells
FeaturePlot(core_immu_scvi, c(
  "TPSB2", "TPSAB1", "CTSG", "HPGD", "GATA2",
  "CPA3", "KIT", "CMA1", "TMOD1", "MEIS2", "SLC24A3"
))

# nonclassical Monocyte
FeaturePlot(core_immu_scvi, c("SMIM25", "FCGR3A", "LILRA5", "LILRB2")) # 35

# plasmacytoid DC (pDC) not present
FeaturePlot(core_immu_scvi, c("TCF4", "JCHAIN", "LILRA4", "PTCRA"))
FeaturePlot(core_immu_scvi, c("CLEC4C", "IRF7", "LILRA4", "TCF4"))

# Erythrocytes
FeaturePlot(core_immu_scvi, c("HBA2", "HBB", "HBA1", "HBG2"))

# Low quality
FeaturePlot(core_immu_scvi, c("CALD1", "MYLK", "KRT17", "ACTA2", "KRT1", "KRT14"))

# Naive T cells not separated as a cluster
FeaturePlot(core_immu_scvi, c("SELL", "LEF1", "CD8A", "NELL2"))

# gamma-delta T cells (GD-T cells)
FeaturePlot(core_immu_scvi, c("FXYD2", "XCL1", "TRGC2", "ZNF683"))
FeaturePlot(core_immu_scvi, c("KLRC3", "KLRC2"))

# NK cells
FeaturePlot(core_immu_scvi, c("SPON2", "PRF1", "FGFBP2", "FCGR3A")) # FGFBP2+
FeaturePlot(core_immu_scvi, c("XCL1", "XCL2", "KLRC1", "GZMK")) # XCL2+
FeaturePlot(core_immu_scvi, c("SPINK2", "MB", "TNPRSS11E", "KLRC1")) # SPINK2+

FeaturePlot(core_immu_scvi, c(rownames(markers)[1:12]))
FeaturePlot(core_immu_scvi, c(rownames(markers)[13:24]))

core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(6,2,12), "CD8+ T cell", core_immu_scvi$seurat_clusters)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(21,63,56), "GD-T cell", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(
  core_immu_scvi$seurat_clusters %in% c(
    32,19,43,20,24,58,0,4,18,31,17,54,14,1,7,27,3,13
  ),
  "CD4+ T cell",
  core_immu_scvi$inherited_celltype_lvl_4
)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters  %in% c(28,8,9), "Reg. T cell", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters == 62, "Plasma cell", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters == 68, "Naive B cell & pDC", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(36,53,45,15,22,69), "LC", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(65,50), "Prolif. DC", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(59,30,38,26,16,25,48,51), "cDC", core_immu_scvi$inherited_celltype_lvl_4)
# Subtype of cDC, hence assgined to lvl_5
core_immu_scvi$inherited_celltype_lvl_5 <- ifelse(core_immu_scvi$seurat_clusters %in% c(59,30,38,26,16,25), "cDC2", core_immu_scvi$seurat_clusters)
core_immu_scvi$inherited_celltype_lvl_5 <- ifelse(core_immu_scvi$seurat_clusters %in% c(51,48), "cDC1", core_immu_scvi$inherited_celltype_lvl_5)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(55) , "Mature DC", core_immu_scvi$inherited_celltype_lvl_4)
# Low quality cluster
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(29,47,61,40,44,52), "LQ", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(37,23,42), "Anti-inflammatory Mph", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(57,49,39), "Inflammatory Mph", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(64,67), "TREM2+ Mph", core_immu_scvi$inherited_celltype_lvl_4)
# Subtype of TREM2+ Mph, hence assgined to lvl_5
core_immu_scvi$inherited_celltype_lvl_5 <- ifelse(core_immu_scvi$seurat_clusters == 64, "C3+ Mph", core_immu_scvi$inherited_celltype_lvl_5)
core_immu_scvi$inherited_celltype_lvl_5 <- ifelse(core_immu_scvi$seurat_clusters == 67, "LPL+ Mph", core_immu_scvi$inherited_celltype_lvl_5)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters == 46, "Monocyte_3", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters == 60, "Neutrophil_3", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(33,34,5,10,11,70), "Mast cell_3", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters == 35, "FGFBP2+ NK", core_immu_scvi$inherited_celltype_lvl_4)
core_immu_scvi$inherited_celltype_lvl_4 <- ifelse(core_immu_scvi$seurat_clusters %in% c(66,41), "XCL2+ NK", core_immu_scvi$inherited_celltype_lvl_4)
# Subtype of XCL2+ NK, hence assgined to lvl_5
core_immu_scvi$inherited_celltype_lvl_5 <- ifelse(core_immu_scvi$seurat_clusters == 41, "SPINK2+ NK", core_immu_scvi$inherited_celltype_lvl_5)
core_immu_scvi$inherited_celltype_lvl_5 <- ifelse(core_immu_scvi$seurat_clusters == 66, "GZMK+ NK", core_immu_scvi$inherited_celltype_lvl_5)

# Plot annotations
DimPlot(core_immu_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_immu_scvi, group.by = "inherited_celltype_lvl_5", label = TRUE)

saveRDS(core_immu_scvi, "./Core_subsets/raw/immune_core.rds")
