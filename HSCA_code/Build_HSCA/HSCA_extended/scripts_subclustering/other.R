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
# Subclustering: other cells
# ===============================================

# UMAP overview to find other_cells lineage cells
DimPlot(hsca_extended, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(hsca_extended, group.by = "sample", label = TRUE, raster = FALSE)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"))

# Subset to other_cells lineage
other_cells <- subset(hsca_extended, subset = seurat_clusters %in% c(101,64,17,88,105,75,83,108,91,61,49,104))

other_cells
atlas

DimPlot(other_cells, group.by = "seurat_clusters", label = T, raster = T)

other_cells$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
other_cells[["RNA"]] <- split(other_cells[["RNA"]], f = other_cells$Dataset)

other_cells
# ==== Re-process immune subset with scVI (Subclustering step) ====
other_cells_scvi <- process_with_scvi(other_cells)
other_cells_scvi <- FindClusters(other_cells_scvi, resolution = 3)

DimPlot(other_cells_scvi, group.by = "seurat_clusters", label = T)
DimPlot(other_cells_scvi, group.by = "seurat_clusters_full", label = T)
DimPlot(other_cells_scvi, group.by = "Dataset", label = T)
DimPlot(other_cells_scvi, group.by = "inherited_celltype_lvl_4", label = T, repel = T)

other_cells_scvi$seurat_clusters %>% table()

FeaturePlot(other_cells_scvi, c("percent_mito", "nFeature_RNA"))
FeaturePlot(other_cells_scvi, c("inherited_celltype_lvl_5_transfer_uncert"), raster = F)
other_cells_scvi

other_cells_scvi <- JoinLayers(other_cells_scvi)

other_cells_scvi

# ==== Marker gene discovery  ====
DefaultAssay(other_cells_scvi)
Idents(other_cells_scvi) <- other_cells_scvi$seurat_clusters
markers <- FindMarkers(
  object = other_cells_scvi,
  ident.1 = 64,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(other_cells_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(other_cells_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(other_cells_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(other_cells_scvi, c("percent_mito", "nFeature_RNA"))

# ==== Feature Plots and Cell type annotation ====

# Melanocytes
FeaturePlot(other_cells_scvi, features = c("MLANA", "QPCT", "TYRP1", "PMEL"))
# Merkel cells
FeaturePlot(other_cells_scvi, features = c("CCER2", "CCK", "KCNMB2", "MIAT", "KRT8", "KRT18", "KRT20"))
# Erythrocyte
FeaturePlot(other_cells_scvi, features = c("ALAS2", "HBA1", "HBB", "HBA2"))
# Adipocyte
FeaturePlot(other_cells_scvi, features = c("FABP4", "ADIPOQ", "RBP4", "PLIN4"))
# Schwann cells
FeaturePlot(other_cells_scvi, features = c("PRX", "GLDN", "MPZ", "MLIP"))
# Neuron
FeaturePlot(other_cells_scvi, features = c("NRXN1", "STARD13", "XKR4", "NTM"))
# Skeletal muscle
FeaturePlot(other_cells_scvi, features = c("MYL1", "ACTA1", "CKM", "COX6A2"))
# CNTNAP2+ neuron
FeaturePlot(other_cells_scvi, features = c("CNTNAP2", "CSMD1", "CNTN5", "PTPRD"))

# Chondrocytes
FeaturePlot(other_cells_scvi, features = c("ACAN", "COL2A1"))
FeaturePlot(other_cells_scvi, features = c("KRT15", "KRT1"))
FeaturePlot(other_cells_scvi, features = c("DCX", "STMN2"))
FeaturePlot(other_cells_scvi, features = c("TOP2A", "MKI67"))
FeaturePlot(other_cells_scvi, features = c("CD74", "HLA-DRA"))
FeaturePlot(other_cells_scvi, features = c("S100A4", "HLA-DRA"))

# Schwann cell markers from Chu_Midha paper
FeaturePlot(other_cells_scvi, features=c("S100B", "NES", "JUN", "POU3F1", "SOX10", "MPZ", "GAP43", "CDH19"))

other_cells_scvi

other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(41,37,44,11,64,39,13,32,38,45,54,57,35,22,12,52,10,42,3,62,19,
                                                                                                   14,15,26,16,43,2,9,55,34,27,59,36,50,58), "Melanocyte_2", other_cells_scvi$seurat_clusters)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(29,56,9,7,60,63), "Neuron_2", other_cells_scvi$inherited_celltype_lvl_4_extended)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(28), "Schwann cell_2", other_cells_scvi$inherited_celltype_lvl_4_extended)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(33,51), "Neural progenitor_2", other_cells_scvi$inherited_celltype_lvl_4_extended)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(6,20), "Chondrocyte_2", other_cells_scvi$inherited_celltype_lvl_4_extended)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(49,48,40), "Erythrocyte_2", other_cells_scvi$inherited_celltype_lvl_4_extended)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(17,5,31,18,8,0,30,1,24,25), "Schwann cell (Ex vivo)_2", other_cells_scvi$inherited_celltype_lvl_4_extended)
other_cells_scvi$inherited_celltype_lvl_4_extended <- ifelse(other_cells_scvi$seurat_clusters %in% c(21,23,4,53,47,61,46), "LQ", other_cells_scvi$inherited_celltype_lvl_4_extended)

DimPlot(other_cells_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T)
DimPlot(other_cells_scvi, group.by = "seurat_clusters", label = T)

saveRDS(other_cells_scvi, "./Extended_subsets/raw/other_cells_extended.rds")
