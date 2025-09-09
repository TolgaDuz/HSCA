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

# ===============================================
# Subclustering: Endothelial cells
# ===============================================

# UMAP overview to find endothelial cells
DimPlot(hsca_extended, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(hsca_extended, group.by = "sample", label = TRUE, raster = FALSE)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"))

# Low quality removal
hsca_extended <- subset(hsca_extended, subset = seurat_clusters == 26, invert = T)
# Remove samples with strong wound healing signatures
hsca_extended <- subset(hsca_extended, subset = sample %in% c("GSM8238438", "GSM8238437"), invert = T)
hsca_extended

# subset endothelial cells
extended_ec <- subset(
  hsca_extended,
  subset = seurat_clusters %in% c(27,39,10,67,41,42,14)
)


DimPlot(extended_ec, group.by = "seurat_clusters", label = T, raster = F)

extended_ec$Dataset %>% table()

# Remove datasets with very few cells
extended_ec <- subset(
  extended_ec,
  subset = Dataset %in% c(
    "Cheng_Cho_2018_1",
    "Cheng_Cho_2018_2",
    "Kim_Krüger_2022",
    "Sun_Liu_2022",
    "Takahashi_Lowry_2019"
  ),
  invert = T
)

extended_ec$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_ec[["RNA"]] <- split(extended_ec[["RNA"]], f = extended_ec$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_ec_scvi <- process_with_scvi(extended_ec)
extended_ec_scvi <- FindClusters(extended_ec_scvi, resolution = 1.7)

DimPlot(extended_ec_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_ec_scvi, group.by = "Dataset", label = T)
FeaturePlot(extended_ec_scvi, c("percent_mito", "nFeature_RNA"))
FeaturePlot(extended_ec_scvi, c("ACT"), raster = F)

extended_ec_scvi <- JoinLayers(extended_ec_scvi)

DimPlot(extended_ec_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_ec_scvi, group.by = "seurat_clusters_full", label = T)
DimPlot(extended_ec_scvi, group.by = "Dataset", label = T)

# ==== Marker gene discovery  ====
DefaultAssay(extended_ec_scvi)
Idents(extended_ec_scvi) <- extended_ec_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_ec_scvi,
  ident.1 = 41,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(extended_ec_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(extended_ec_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_ec_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(extended_ec_scvi, c("percent_mito", "nFeature_RNA"))

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

extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(13) , "RSG5+ EC", extended_ec_scvi$seurat_clusters)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(18,17,26,9,30,6,5,31,10,37,43,20,24,22,11,1,19,12,40,39,16,45,32), "Venous 1 EC", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(35,29,27,48,46,25,41,2,3,28), "Venous 2 EC", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(33), "SCG3+ LEC", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(44), "NEO1+ LEC", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(0,7,42,36,34), "LYVE1+ LEC", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(15,14,32,23), "Capillary EC_3", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(21,4,47,8), "Arterial EC_3", extended_ec_scvi$inherited_celltype_lvl_4_extended)
extended_ec_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_ec_scvi$seurat_clusters %in% c(49,38), "LQ", extended_ec_scvi$inherited_celltype_lvl_4_extended)

DimPlot(extended_ec_scvi, group.by = "seurat_clusters", label = T)

saveRDS(extended_ec_scvi, "./Extended_subsets/raw/ec_extended.rds")
