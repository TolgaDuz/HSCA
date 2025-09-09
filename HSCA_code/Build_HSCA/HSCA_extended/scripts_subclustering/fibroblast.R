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
# Subclustering: Fibroblasts
# ===============================================

# UMAP overview to find fibroblast cells
DimPlot(hsca_extended, group.by = "seurat_clusters", label = TRUE, raster = FALSE)
DimPlot(hsca_extended, group.by = "sample", label = TRUE, raster = FALSE)
FeaturePlot(hsca_extended, c("percent_mito", "nFeature_RNA"))

# Low quality removal
hsca_extended <- subset(hsca_extended, subset = seurat_clusters == 26, invert = T)
# Remove samples with strong wound healing signatures
hsca_extended <- subset(hsca_extended, subset = sample %in% c("GSM8238438", "GSM8238437"), invert = T)
hsca_extended

# Firbo subset
extended_fb <- subset(
  hsca_extended,
  subset = seurat_clusters %in% c(4,82,28,71,73,66,0,31,9,12,5,24,106,76,33,55,86,93)
)

DimPlot(extended_fb, group.by = "seurat_clusters", label = T, raster = T)
FeaturePlot(extended_fb, "CCER2", raster = F)

extended_fb$Dataset %>% table()

# Remove datasets with very few cells
extended_fb <- subset(extended_fb, subset = Dataset %in% c(
  "Cheng_Cho_2018_1", "Cheng_Cho_2018_2", "Chu_Midha_2022",
  "Ganier_Lynch_2024_2", "Kim_Krüger_2022", "Sun_Liu_2022"
  ),
  invert = T
)

extended_fb$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_fb[["RNA"]] <- split(extended_fb[["RNA"]], f = extended_fb$Dataset)

extended_fb

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_fb_scvi <- process_with_scvi(extended_fb)
extended_fb_scvi <- FindClusters(extended_fb_scvi, resolution = 3)

extended_fb_scvi
DimPlot(extended_fb_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_fb_scvi, group.by = "seurat_clusters_full", label = T)
DimPlot(extended_fb_scvi, group.by = "Dataset", label = T)
DimPlot(extended_fb_scvi, group.by = "inherited_celltype_lvl_4", label = T, repel = T)

extended_fb_scvi$seurat_clusters %>% table()

FeaturePlot(extended_fb_scvi, c("percent_mito", "nFeature_RNA"))
FeaturePlot(extended_fb_scvi, c("inherited_celltype_lvl_5_transfer_uncert"), raster = F)
extended_fb_scvi

extended_fb_scvi <- JoinLayers(extended_fb_scvi)

extended_fb_scvi

# ==== Marker gene discovery  ====
DefaultAssay(extended_fb_scvi)
Idents(extended_fb_scvi) <- extended_fb_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_fb_scvi,
  ident.1 = 64,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(extended_fb_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(extended_fb_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_fb_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(extended_fb_scvi, c("percent_mito", "nFeature_RNA"))

# ==== Feature Plots and Cell type annotation ====
FeaturePlot(extended_fb_scvi, c("IFI27"))
FeaturePlot(extended_fb_scvi, c("COL6A5", "COL18A1", "CCL2", "CCL19"))

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
FeaturePlot(extended_fb_scvi, features = c("SFRP2", "ELN", "MMP2", "QPCT"))
# A1
FeaturePlot(extended_fb_scvi, features = c("IGFBP6", "PI116", "SLPI", "CCN5"))   # Paper 1
FeaturePlot(extended_fb_scvi, features = c("WISP2", "SEMA3B", "LGR5", "ANGPTL5"))  # Paper 2
# A2
FeaturePlot(extended_fb_scvi, features = c("APCDD1", "COL18A1", "COMP", "NKD2"))  # Paper 1
FeaturePlot(extended_fb_scvi, features = c("HSPB3", "COL18A1", "COL6A5", "NKD2")) # Paper 2
# A3
FeaturePlot(extended_fb_scvi, features = c("ELN", "RGCC", "SGCA", "WIF1"))        # Paper 1
FeaturePlot(extended_fb_scvi, features = c("SOSTDC1", "CORIN", "SGCA", "WIF1"))  # Paper 2
# A4
FeaturePlot(extended_fb_scvi, features = c("FBN1", "PCOLCE2", "PRG4", "SFRP4"))  # Paper 1
FeaturePlot(extended_fb_scvi, features = c("C1QTNF3", "SCARA5", "PRG4", "TRAC")) # Paper 2

# B
FeaturePlot(extended_fb_scvi, c("APOE", "C7", "CYGB", "IGFBP7"))
# B1
FeaturePlot(extended_fb_scvi, c("CCL2", "ITM2A", "SPSB1", "TNFAIP6")) # Paper 1
FeaturePlot(extended_fb_scvi, c("GEM", "CXCL2", "CXCL1", "TNFAIP6")) # Paper 2
# B2
FeaturePlot(extended_fb_scvi, c("CCDC146", "CCL19", "CD74", "TNFSF13B")) # Paper 1
FeaturePlot(extended_fb_scvi, c("GGT5", "IL33", "C7", "SCN4B")) # Paper 2
# B3
FeaturePlot(extended_fb_scvi, c("CCL19", "CTSH", "RBP5", "ACHE")) # only in paper 2
# B4
FeaturePlot(extended_fb_scvi, c("EFEMP1", "ITM2A", "MYOC", "GDF10")) # only in paper 2

# C
FeaturePlot(extended_fb_scvi, c("SFRP1", "TNMD", "DKK3", "TNN"))
# C1 (DS)
FeaturePlot(extended_fb_scvi, c("COL11A1", "DPEP1", "TNMD", "WFDC1")) # only in paper 1
FeaturePlot(extended_fb_scvi, c("COL11A1", "MEF2C", "DPEP1", "WFDC1")) # only in paper 2
# C2 (Upper Bulge DP)
FeaturePlot(extended_fb_scvi, c("COCH", "CRABP1", "FIBIN", "RSPO4")) # only in paper 1
FeaturePlot(extended_fb_scvi, c("COCH", "CRABP1", "NDNF", "SLITRK6")) # only in paper 2
# C3
FeaturePlot(extended_fb_scvi, c("ASPN", "F2R", "GPM6B", "POSTN")) # only in paper 1
FeaturePlot(extended_fb_scvi, c("LTBP2", "LRRC15", "BGN", "POSTN")) # only in paper 2
# C4
FeaturePlot(extended_fb_scvi, c("ANGPTL7", "APOD", "ECRG4", "TM4SF1")) # only in paper 1
FeaturePlot(extended_fb_scvi, c("ANGPTL7", "APOD", "ECRG4", "TM4SF1")) # only in paper 2
# C5
FeaturePlot(extended_fb_scvi, c("IGFBP3", "SLC5A3", "WNT5A", "LUZP2")) # only in paper 2
FeaturePlot(extended_fb_scvi, features = c("PGM2L1", "PTCH1", "KIF26B", "GRIK1"))
FeaturePlot(extended_fb_scvi, c("IGFBP3", "DKK2", "RSPO4", "WIF1", "EDN3", "DIO3"))

# D
# D1
FeaturePlot(extended_fb_scvi, c("ANGPTL7", "ENTPD2", "CDH19", "ATP1A2"))
# D2
FeaturePlot(extended_fb_scvi, c("BNC2", "ITGA6", "ITGB4", "TNNC1"))
# E1
FeaturePlot(extended_fb_scvi, features = c("IGFBP2", "COL26A1", "WNT2", "RAMP1"))

# Interesting
# COL8A1
FeaturePlot(extended_fb_scvi, features = c("COL8A1", "PTX3"))
# PDZRN4
FeaturePlot(extended_fb_scvi, features = c("DKK1", "PDZRN4", "SGCZ"))

# CAFs
FeaturePlot(extended_fb_scvi, c("NDUFA4L2", "GMFG", "MCAM", "PDGFA", "TAGLN", "COL4A1", "TPM2", "LOXL2", "POSTN")) # 21 kommt nur in einem datensatz vor!
FeaturePlot(extended_fb_scvi, c("RRAD", "CHI3L1", "CA12", "WISP1", "WNT5A", "BCAT1", "CD82", "TGFBI", "SLC16A3")) # 17 nur ein patient!



# Cluster 77 consists keratinocytes but also includes Merkel cells,
# which we separate here for downstream annotation and visualization
merkel_cells_and_kera <- subset(extended_fb_scvi, subset = seurat_clusters %in% c(77))

# Check where the annotated merkel cells from the HSCA core cluster
DimPlot(merkel_cells_and_kera, group.by = "Core")
merkel_cells_and_kera <- FindClusters(merkel_cells_and_kera, resolution = 1)
DimPlot(merkel_cells_and_kera)
merkel_cells_and_kera$seurat_clusters %>% table()
FeaturePlot(atlas, "CCER2")

# Annotate cluster identity (Merkel cells)
# and later transfer updated labels back into the fibroblast object
merkel_cells_and_kera$inherited_celltype_lvl_4_extended <- ifelse(merkel_cells_and_kera$seurat_clusters %in% c(5), "Merkel cell_3", merkel_cells_and_kera$seurat_clusters)
merkel_cells_and_kera$inherited_celltype_lvl_4_extended <- ifelse(merkel_cells_and_kera$seurat_clusters %in% c(0,1,2,3,4), "LQ", merkel_cells_and_kera$inherited_celltype_lvl_4_extended)
DimPlot(merkel_cells_and_kera, group.by = "inherited_celltype_lvl_4_extended")


### for classifier ###
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(1,49,68,16,7,38,83,22,20,60,90,78,87,54), "A1", extended_fb_scvi$seurat_clusters)
# lvl 5
extended_fb_scvi$inherited_celltype_lvl_5_extended  <- ifelse(extended_fb_scvi$seurat_clusters %in% c(78), "PDZRN4+ FB", extended_fb_scvi$seurat_clusters)

extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(23,66,8,18,56,81,21,62,63,29,26,52), "A2", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(71,46,27,12,35,82,84,25,31,69), "A3", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(64,11,0,4,43,44,42,86,57), "A4", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(65,58,34,72,33,19,3,59,37,39,17,55,67), "B1", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(24,2,10,53,28), "B3", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(5,85,50), "B2", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(32,88,45,80,73,51,40), "B4", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(48), "D1", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(70,61), "D2", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(41,47,15), "Dermal sheath (C1)", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(13,36,79,14), "C3", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(9), "Outer bulge DP (C2)", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(6), "RAMP1+ extended_fb (E1)", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(30), "Anagen DP (C5)", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(74), "CAF1", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(75), "CAF2", extended_fb_scvi$inherited_celltype_lvl_4_extended)
extended_fb_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_fb_scvi$seurat_clusters %in% c(76,40,89), "LQ", extended_fb_scvi$inherited_celltype_lvl_4_extended)

# Transfer merkel cell labels
common_cells <- intersect(Cells(extended_fb_scvi), Cells(merkel_cells_and_kera))

extended_fb_scvi$inherited_celltype_lvl_4_extended[common_cells] <-
  merkel_cells_and_kera$inherited_celltype_lvl_4_extended[common_cells]

DimPlot(extended_fb_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T)
DimPlot(extended_fb_scvi, group.by = "seurat_clusters", label = T)

saveRDS(extended_fb_scvi, "./Extended_subsets/raw/fb_extended.rds")
