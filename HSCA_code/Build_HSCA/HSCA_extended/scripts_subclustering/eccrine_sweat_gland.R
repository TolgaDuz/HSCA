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
# Subclustering: Eccrine sweat gland cells
# ===============================================

# UMAP overview to find esg cells
DimPlot(hsca_extended, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(hsca_extended, group.by = "sample", label = TRUE, raster = FALSE)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"))

# Low quality removal
hsca_extended <- subset(hsca_extended, subset = seurat_clusters == 26, invert = T)
# Remove samples with strong wound healing signatures
hsca_extended <- subset(hsca_extended, subset = sample %in% c("GSM8238438", "GSM8238437"), invert = T)
hsca_extended

# ecne subset
extended_esg <- subset(
  hsca_extended,
  subset = seurat_clusters %in% c(57,58,102,90,22,99)
)

extended_esg

DimPlot(extended_esg, group.by = "seurat_clusters", label = T, raster = F)

FeaturePlot(extended_esg, "CNTNAP2")

extended_esg$Dataset %>% table()

extended_esg <- subset(
  extended_esg,
  subset = Dataset %in% c(
    "Kim_Krüger_2022",
    "Chu_Midha_2022",
    "Ahlers_Siracusa_2022",
    "Ganier_Lynch_2024_2"
  ),
  invert = T
)

extended_esg$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_esg[["RNA"]] <- split(extended_esg[["RNA"]], f = extended_esg$Dataset)

extended_esg

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_esg_scvi <- process_with_scvi(extended_esg)
extended_esg_scvi <- FindClusters(extended_esg_scvi, resolution = 2)

DimPlot(extended_esg_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_esg_scvi, group.by = "Dataset", label = T)

FeaturePlot(extended_esg_scvi, c("percent_mito", "nFeature_RNA"))

extended_esg_scvi

extended_esg_scvi <- JoinLayers(extended_esg_scvi)

extended_esg_scvi

# ==== Marker gene discovery  ====
DefaultAssay(extended_esg_scvi)
Idents(extended_esg_scvi) <- extended_esg_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_esg_scvi,
  ident.1 = 64,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(extended_esg_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(extended_esg_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_esg_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(extended_esg_scvi, c("percent_mito", "nFeature_RNA"))

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

# CNTNAP2+ Neuron
FeaturePlot(extended_esg_scvi, c("CNTNAP2", "CNTN5", "CSMD1", "PTPRD"))

DimPlot(extended_esg_scvi, group.by = "seurat_clusters", label = T)

extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(34,4,20,21,17,12,19,22,27), "Dark cell", extended_esg_scvi$seurat_clusters)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(13,31), "Myoepithelial", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(7,14,36), "LQ", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(15,23,29), "S100P+ luminal", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(37,28), "S100P- luminal", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(32,2,18,0,35,6,11,8,3,5,10,16,33,25), "Basal luminal", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(1), "Clear cell 1", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(38), "CNTNAP2+ Neuron", extended_esg_scvi$inherited_celltype_lvl_4_extended)
extended_esg_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_esg_scvi$seurat_clusters %in% c(26,24,9,30), "Clear cell 2", extended_esg_scvi$inherited_celltype_lvl_4_extended)

DimPlot(extended_esg_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T)

saveRDS(extended_esg_scvi, "./Extended_subsets/raw/esg_extended.rds")
