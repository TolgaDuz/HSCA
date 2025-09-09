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
# Subclustering: Immune cells
# ===============================================

# UMAP overview to find immune cells
DimPlot(hsca_extended, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(hsca_extended, group.by = "sample", label = TRUE, raster = FALSE)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"))

# Low quality removal
hsca_extended <- subset(hsca_extended, subset = seurat_clusters == 26, invert = T)
# Remove samples with strong wound healing signatures
hsca_extended <- subset(hsca_extended, subset = sample %in% c("GSM8238438", "GSM8238437"), invert = T)
hsca_extended

# subset immune cells
extended_immu <- subset(
  hsca_extended,
  subset = seurat_clusters %in% c(
    40,11,50,15,19,32,29,79,60,89,96,72,69,98,95,48,20,92,23,46,70,94,30
  ),
  invert = F
)

DimPlot(extended_immu, group.by = "seurat_clusters", label = T, raster = F)

extended_immu$Dataset %>% table()

# Remove datasets with very few cells
extended_immu <- subset(extended_immu, subset = Dataset %in% c(
  "Chu_Midha_2022"
), invert = T)

extended_immu$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_immu[["RNA"]] <- split(extended_immu[["RNA"]], f = extended_immu$Dataset)

extended_immu

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_immu_scvi <- process_with_scvi(extended_immu)
extended_immu_scvi <- FindClusters(extended_immu_scvi, resolution = 2.4)

DimPlot(extended_immu_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_immu_scvi, group.by = "seurat_clusters_full", label = T)
DimPlot(extended_immu_scvi, group.by = "Dataset", label = T)
DimPlot(extended_immu_scvi, group.by = "inherited_celltype_lvl_4", label = T, repel = T)

extended_immu_scvi$seurat_clusters %>% table()

FeaturePlot(extended_immu_scvi, c("percent_mito", "nFeature_RNA"))
FeaturePlot(extended_immu_scvi, c("inherited_celltype_lvl_5_transfer_uncert"), raster = F)
extended_immu_scvi

extended_immu_scvi <- JoinLayers(extended_immu_scvi)

extended_immu_scvi

# ==== Marker gene discovery  ====
DefaultAssay(extended_immu_scvi)
Idents(extended_immu_scvi) <- extended_immu_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_immu_scvi,
  ident.1 = 64,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(extended_immu_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(extended_immu_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_immu_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(extended_immu_scvi, c("percent_mito", "nFeature_RNA"))

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
FeaturePlot(core_immu_scvi, c("S100A9", "FCN1", "AQP9", "VCAN", "EREG")) # Monocytes
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

# plasmacytoid DC (pDC)
FeaturePlot(core_immu_scvi, c("TCF4", "JCHAIN", "LILRA4", "PTCRA"))
FeaturePlot(core_immu_scvi, c("CLEC4C", "IRF7", "LILRA4", "TCF4"))

# Erythrocytes
FeaturePlot(core_immu_scvi, c("HBA2", "HBB", "HBA1", "HBG2"))

# Naive T cells not separated as a cluster
FeaturePlot(core_immu_scvi, c("SELL", "LEF1", "CD8A", "NELL2"))

# gamma-delta T cells (GD-T cells)
FeaturePlot(core_immu_scvi, c("FXYD2", "XCL1", "TRGC2", "ZNF683"))
FeaturePlot(core_immu_scvi, c("KLRC3", "KLRC2"))

# NK cells
FeaturePlot(core_immu_scvi, c("SPON2", "PRF1", "FGFBP2", "FCGR3A")) # FGFBP2+
FeaturePlot(core_immu_scvi, c("XCL1", "XCL2", "KLRC1", "GZMK")) # XCL2+
FeaturePlot(core_immu_scvi, c("SPINK2", "MB", "TNPRSS11E", "KLRC1")) # SPINK2+

# Low quality
FeaturePlot(extended_immu_scvi, c("CALD1", "MYLK", "KRT17", "ACTA2", "KRT1", "KRT14"))

extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(62,43) , "Mature DC", extended_immu_scvi$seurat_clusters)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(51), "Prolif. DC", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(17,25,42), "LQ", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(39,0,19), "Inflammatory DC", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(13,12,10), "Inflammatory Mph", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(36,16), "Anti-inflammatory Mph", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(47,32,46), "LC", extended_immu_scvi$inherited_celltype_lvl_4_extended)

extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(8,7,34,38,18,11,64,19,55,66), "cDC", extended_immu_scvi$inherited_celltype_lvl_4_extended)

# subtype of cDC, therefore assignment to lvl_5
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(8,7,34,38,18,11,64), "cDC2", extended_immu_scvi$seurat_clusters)
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(19,55), "cDC1", extended_immu_scvi$inherited_celltype_lvl_5_extended)
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(66), "pDC", extended_immu_scvi$inherited_celltype_lvl_5_extended)


extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(54,29), "TREM2+ Mph", extended_immu_scvi$inherited_celltype_lvl_4_extended)
# subtype of trem2+, therefore assignment to lvl_5
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters == 29, "C3+ Mph", extended_immu_scvi$inherited_celltype_lvl_5_extended)
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters == 54, "LPL+ Mph", extended_immu_scvi$inherited_celltype_lvl_5_extended)


extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(33,14,15,53), "Mast cell_3", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters  %in% c(50,41), "Naive B cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters  %in% c(69,67,48), "Plasma cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(23,39,58,63), "cMono", extended_immu_scvi$inherited_celltype_lvl_4_extended)

extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(23), "EREG+ Mono", extended_immu_scvi$inherited_celltype_lvl_5_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters == 59, "ncMono", extended_immu_scvi$inherited_celltype_lvl_4_extended)

extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters == 60, "Neutrophil_3", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters  %in% c(4,24), "Reg. T cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(40,27,44), "FGFBP2+ NK", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(35,30), "XCL2+ NK", extended_immu_scvi$inherited_celltype_lvl_4_extended)
# subtype of XCL2+ NK, therefore assignment to lvl_5
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters == 30, "SPINK2+ NK", extended_immu_scvi$inherited_celltype_lvl_5_extended)
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters == 35, "GZMK+ NK", extended_immu_scvi$inherited_celltype_lvl_5_extended)

extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(6,21,20,45, 61), "CD8+ T cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(5,3,37,2,49,31,68,1,9,0,52,26,70,28), "CD4+ T cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)

extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(28), "naive CD4+ T cell", extended_immu_scvi$inherited_celltype_lvl_5_extended)
extended_immu_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(61), "naive CD8+ T cell", extended_immu_scvi$inherited_celltype_lvl_5_extended)


extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(22,65,57), "GD-T cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)
extended_immu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_immu_scvi$seurat_clusters %in% c(56), "Prolif. T cell", extended_immu_scvi$inherited_celltype_lvl_4_extended)

DimPlot(extended_immu_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, repel = T)
DimPlot(extended_immu_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T, repel = T)
DimPlot(extended_immu_scvi, group.by = "seurat_clusters", label =T)

saveRDS(extended_immu_scvi, "./Extended_subsets/raw/immune_extended.rds")
