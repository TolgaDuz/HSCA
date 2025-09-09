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
# Subclustering:IFE + PSU cells together
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


DimPlot(hsca_extended, group.by = "seurat_clusters", label = T, raster = F, repel = F)

extended_epid <- subset(
  hsca_extended,
  subset = seurat_clusters %in% c(
    100,25,38,68,18,2,51,54,21,87,16,62,1,59,80,
    65,8,63,3,47,34,6, 5,78,53,77,97,35,43,74
  )
)

DimPlot(extended_epid, group.by = "seurat_clusters", label = T,raster = F)
DimPlot(extended_epid, group.by = "sample", label = F,raster = F)

extended_epid$Dataset %>% table()

# Remove datasets with very few cells
extended_epid <- subset(
  extended_epid,
  subset = Dataset %in% c(
    "Chu_Midha_2022"
  ),
  invert = T
)

extended_epid$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_epid[["RNA"]] <- split(extended_epid[["RNA"]], f = extended_epid$Dataset)

extended_epid

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_epid_scvi <- process_with_scvi(extended_epid)
extended_epid_scvi <- FindClusters(extended_epid_scvi, resolution = 1.5)

DimPlot(extended_epid_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_epid_scvi, group.by = "seurat_clusters_full", label = T)
DimPlot(extended_epid_scvi, group.by = "Dataset", label = T)
DimPlot(extended_epid_scvi, group.by = "inherited_celltype_lvl_4", label = T, repel = T)

extended_epid_scvi$seurat_clusters %>% table()

FeaturePlot(extended_epid_scvi, c("percent_mito", "nFeature_RNA"))
FeaturePlot(extended_epid_scvi, c("inherited_celltype_lvl_5_transfer_uncert"), raster = F)
extended_epid_scvi

extended_epid_scvi <- JoinLayers(extended_epid_scvi)

extended_epid_scvi

# ==== Marker gene discovery  ====
DefaultAssay(extended_epid_scvi)
Idents(extended_epid_scvi) <- extended_epid_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_epid_scvi,
  ident.1 = 64,
  max.cells.per.ident = 500,
  only.pos = TRUE
)

# Inspect top markers
head(markers, 40)
DimPlot(extended_epid_scvi, group.by = "seurat_clusters", label = T)

FeaturePlot(extended_epid_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_epid_scvi, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(extended_epid_scvi, c("percent_mito", "nFeature_RNA"))


# ===============================================
# Subclustering: Pilosebaceous Unit (PSU)
# ===============================================

# UMAP overview to find PSU cells
DimPlot(extended_epid_scvi, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(extended_epid_scvi, group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

FeaturePlot(extended_epid_scvi, c("percent_mito", "nFeature_RNA"), reduction = "umap")
FeaturePlot(extended_epid_scvi, "KRT6A", raster = F)

# Track clusters
extended_epid_scvi$seurat_clusters_full <- object$seurat_clusters
DimPlot(extended_epid_scvi, reduction = "umap", group.by = "seurat_clusters_full", label = T, raster = F)

extended_psu <- subset(extended_epid_scvi, subset = seurat_clusters %in% c(72,6,61,1,25,45,20,60,70))

extended_psu$orig_celltype_lvl_3 %>% table()
extended_psu$seurat_clusters %>% table()

DimPlot(extended_psu, group.by = "seurat_clusters", label = T)
DimPlot(extended_psu, group.by = "orig_celltype_lvl_3", label = T)

# Remove datasets with very few cells
extended_psu <- subset(
  extended_psu,
  subset = Dataset %in% c(
    "Kim_Nagao_2020_five_prime", "Rustagi_Singh_2022",
    "Ahlers_Siracusa_2022", "Sole-boldo_Lyko_2020",
    "Kim_Krüger_2022"
  ),
  invert = T
)

extended_psu$Dataset %>% table()
extended_psu
# Split RNA assay by dataset for scVI processing
extended_psu[["RNA"]] <- split(extended_psu[["RNA"]], f = extended_psu$Dataset)

extended_psu

# ==== Re-process immune subset with scVI (Subclustering step) ====
extended_psu_scvi <- process_with_scvi(extended_psu)
extended_psu_scvi <- FindClusters(extended_psu_scvi, resolution = 3.5)

DimPlot(extended_psu_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_psu_scvi, group.by = "sample", label = F)
DimPlot(extended_psu_scvi, group.by = "Dataset", label = T)
DimPlot(extended_psu_scvi, group.by = "anatomical_region_level2", label = T)
DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_4", label = T)
FeaturePlot(extended_psu_scvi, "C1QTNF12")

extended_psu_scvi$Dataset %>% table()

extended_psu_scvi <- JoinLayers(extended_psu_scvi)

FeaturePlot(extended_psu_scvi, c("percent_mito", "nFeature_RNA"))

# ==== Marker gene discovery  ====
DefaultAssay(extended_psu_scvi)
Idents(extended_psu_scvi) <- extended_psu_scvi$seurat_clusters
markers <- FindMarkers(
  object = extended_psu_scvi,
  ident.1 = 15, # Change
  only.pos = TRUE
)

head(markers, 50)

# Inspect top markers
FeaturePlot(extended_psu_scvi, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(extended_psu_scvi, c(rownames(markers)[13:24]), raster = F)
DimPlot(extended_psu_scvi, group.by = "seurat_clusters", label = T, raster = F)

# Remove LQ
extended_psu_scvi <- subset(extended_psu_scvi, subset = seurat_clusters %in% c(28,44), invert = T)
DimPlot(extended_psu_scvi, group.by = "seurat_clusters", label = T, raster = F)

# ==== Feature Plots and Cell type annotation ====
FeaturePlot(extended_psu_scvi, c("SERPINB4", "KRT1", "PTN", "C1QTNF12"))
FeaturePlot(extended_psu_scvi, c("KRT15", "KRT1", "MKI67", "FLG"))
FeaturePlot(extended_psu_scvi, c("COL17A1", "KRT14", "POSTN", "KRT5"))
FeaturePlot(extended_psu_scvi, c("DIO2", "WIF1", "FRZB", "SHROOM3"))
FeaturePlot(extended_psu_scvi, c("S100A8", "S100A7", "S100A9", "KRT6A"))
FeaturePlot(extended_psu_scvi, c("PTN", "GATA6", "KRT79", "LGR6"))
FeaturePlot(extended_psu_scvi, c("C1QTNF12", "GATA6", "KRT79", "PTN"))
FeaturePlot(extended_psu_scvi, c("CALB2", "PI3", "KRT79", "IVL"))
FeaturePlot(extended_psu_scvi, c("KLK6", "LCN2", "RHCG", "CDA"))
FeaturePlot(extended_psu_scvi, c("CRAT", "SEC14L6", "PLIN5", "PNPLA3"))
FeaturePlot(extended_psu_scvi, c("NNAT", "IL1R2", "KRT7", "WFDC2"))
FeaturePlot(extended_psu_scvi, c("MKI67", "TOP2A"))
FeaturePlot(extended_psu_scvi, c("LGR5", "RARB", "TBX2", "LHX2"))
FeaturePlot(extended_psu_scvi, c("RXRA", "RARG"))
FeaturePlot(extended_psu_scvi, c("FST", "LGR6"))

# DimPlots with/without repel
DimPlot(extended_psu_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = TRUE)
DimPlot(extended_psu_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = FALSE)

extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(1,58,49), "SG progenitors", extended_psu_scvi$seurat_clusters)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(39,25,57), "Supr. SG", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(51,40,3,11,52), "Basal SG", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(9,5,14,43), "Basal KC_3", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(15), "Basal SD", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(5,56,60,0,18,17,21,45,2,46,10), "Basal JZ", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(19), "Supr. JZ & SD", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(36,62), "Supr. Infund. KC", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(26,8,23,38,47,7,37,24,13,35), "Infund. KC", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(34,55,22,64,63,53,29,31,33,42,30,50,16,27,12,41,32,20,54,6,4,59,48), "lower HF", extended_psu_scvi$inherited_celltype_lvl_4_extended)
extended_psu_scvi$inherited_celltype_lvl_4_extended <- ifelse(extended_psu_scvi$seurat_clusters %in% c(61), "LQ", extended_psu_scvi$inherited_celltype_lvl_4_extended)

DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)

saveRDS(extended_psu_scvi, "./Extended_subsets/raw/psu_extended.rds")


# ===============================================
# Subclustering: Lower hairfollicle cluster from PSU
# ===============================================

DefaultAssay(extended_psu_scvi) <- "RNA"
extended_psu_scvi[["RNA"]] <- JoinLayers(extended_psu_scvi[["RNA"]])
extended_psu_scvi

extended_psu_scvi

# Quick overview
DimPlot(extended_psu_scvi, group.by = "seurat_clusters", label = TRUE, raster = FALSE)

# subset to lower hairfollicle
extended_lower_hf <- subset(
  extended_psu_scvi,
  subset = seurat_clusters %in% c(
    34,55,22,64,63,53,29,31,33,42,30,50,16,27,12,41,32,20,54,6,4,59,48
  ),
  invert = F
)

extended_lower_hf

DimPlot(extended_lower_hf, group.by = "seurat_clusters", label = T)
DimPlot(extended_lower_hf, group.by = "inherited_celltype_lvl_3", label = T)
extended_lower_hf$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_lower_hf[["RNA"]] <- split(extended_lower_hf[["RNA"]], f = extended_lower_hf$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
set.seed(42)
extended_lower_hf_scvi <- process_with_scvi(extended_lower_hf)
extended_lower_hf_scvi <- FindClusters(extended_lower_hf_scvi, resolution = 8)

# Visualize clusters and metadata
DimPlot(extended_lower_hf_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(extended_lower_hf_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(extended_lower_hf_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(extended_lower_hf_scvi, group.by = "anatomical_region_level2", label = TRUE)
DimPlot(extended_lower_hf_scvi, group.by = "Dataset", label = FALSE)

extended_lower_hf_scvi$inherited_celltype_lvl_3 %>% table()

FeaturePlot(extended_lower_hf_scvi, c("percent_mito", "nFeature_RNA"))

extended_lower_hf_scvi <- JoinLayers(extended_lower_hf_scvi)

# ==== Marker gene discovery ====
markers <- FindMarkers(extended_lower_hf_scvi, ident.1 = 40, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 20)

# FeaturePlots of top markers
FeaturePlot(extended_lower_hf_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(extended_lower_hf_scvi, rownames(markers)[13:24], raster = FALSE)
FeaturePlot(extended_lower_hf_scvi, rownames(markers)[21:30], raster = FALSE)

# ==== Feature Plots and Cell type annotation ====

# Huxley's layer
FeaturePlot(extended_lower_hf_scvi, c("TGFA", "CDSN", "SPINK5", "CLDN4"))
# Henle's layer (lower)
FeaturePlot(extended_lower_hf_scvi, c("SLPI", "RHCG", "KRT71", "CAPN8"))
# Henle's layer (middle)
FeaturePlot(extended_lower_hf_scvi, c("CRCT1", "CPA4", "CNFN", "ASPRV1"))
# Henle's layer (upper)
# Caution: comp and selenop also in ORS
FeaturePlot(extended_lower_hf_scvi, c("CHI3L1", "KRT75", "COMP", "SELENOP"))
# Late Matrix
FeaturePlot(extended_lower_hf_scvi, c("LPAR6", "MYCN", "FABP5", "SCEL"))
# Early Matrix
FeaturePlot(extended_lower_hf_scvi, c("HIST1H1B", "TOP2A", "MYCN", "HIST1H1D"))
# Medulla
FeaturePlot(extended_lower_hf_scvi, c("KRTAP10-2", "SCYGR4", "KRTAP12-1", "SCYGR8"))
FeaturePlot(extended_lower_hf_scvi, c("VSIG8", "KRT86", "KRT83", "KRTAP5-8"))
# Cortex
FeaturePlot(extended_lower_hf_scvi, c("SELENBP1", "KRT85", "LY6G6D", "GPNMB")) # LEF1+ Cortex
FeaturePlot(extended_lower_hf_scvi, c("LEF1")) # LEF1+ Cortex
FeaturePlot(extended_lower_hf_scvi, c("KRT81", "KRT83", "KRT86")) # KRT83+ Cortex
# Melanocytes
FeaturePlot(extended_lower_hf_scvi, c("TYRP1", "MLANA", "TYR", "QPCT"))
# DPC
FeaturePlot(extended_lower_hf_scvi, c("IGFBP3", "DKK2", "RSPO4", "WIF1"))
# Cuticle
# IRS cuticle
FeaturePlot(extended_lower_hf_scvi, c("KRT25", "KRT26", "KRT27", "KRT28", "KRT71", "KRT72", "KRT73"))
# Cuticle hair shaft
FeaturePlot(extended_lower_hf_scvi, c("KRT82", "PROCR", "VSIG8", "ACTBL2"))
FeaturePlot(extended_lower_hf_scvi, c("KRT32", "SOX21", "CYP26B1", "PPP2R1B", "PADI3"))
FeaturePlot(extended_lower_hf_scvi, c("KRT32", "KRT35", "KRT39", "KRT40", "KRT82", "KRT85"))
# CL
FeaturePlot(extended_lower_hf_scvi, c("GAL", "GABRP", "SERPINA1", "CALB2"))
FeaturePlot(extended_lower_hf_scvi, c("KRT6A", "SH3D21", "NES", "ADAMTS1"))
FeaturePlot(extended_lower_hf_scvi, c("VEGFA", "CA9", "MEST", "ADAMTS1"))
# Basal ORS
FeaturePlot(extended_lower_hf_scvi, c("LGR5", "COMP", "SELENOP", "ANGPTL7"))
FeaturePlot(extended_lower_hf_scvi, c("FGF18", "CRLF1", "ALAD", "RNF152"))
# Bulge
FeaturePlot(extended_lower_hf_scvi, c("KRT6A", "KRT16", "S100A2", "SERPINA3"))
FeaturePlot(extended_lower_hf_scvi, c("DSG3", "DSP", "PKP1", "FABP5"))
FeaturePlot(extended_lower_hf_scvi, c("CST6", "KRT6A"))
FeaturePlot(extended_lower_hf_scvi, c("MUCL1"))
DFeaturePlot(extended_lower_hf_scvi, c("FABP5"))
FeaturePlot(extended_lower_hf_scvi, c("DIO2", "PDK4", "CXCL14"))
FeaturePlot(extended_lower_hf_scvi, c("COL6A1"))
FeaturePlot(extended_lower_hf_scvi, c("ANGPTL7"))
# lower catagen
FeaturePlot(extended_lower_hf_scvi, c("SULT1E1", "NMB", "THBS1", "CLU", "PAPLN"))
# SHG
FeaturePlot(extended_lower_hf_scvi, c(
  "KREMEN2", "SLC2A3", "EPCAM", "TNF", "SOX11", "CASC15",
  "SPRY4", "SPX5", "ROBO1", "CUX2", "RARB"
))

DimPlot(extended_lower_hf_scvi, group.by = "seurat_clusters", label = T)
DimPlot(extended_lower_hf_scvi, group.by = "seurat_clusters_full", label = T)

extended_lower_hf_scvi <- JoinLayers(extended_lower_hf_scvi)

extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(67,48,59,26), "Matrix (Early)", extended_lower_hf_scvi$seurat_clusters)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(33,44,65), "Matrix (Late) & Henle's layer (lower)", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(62,27,78), "LEF1+ Cortex", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(60,81), "KRT83+ Cortex", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(71), "Cuticle hair shaft", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(47), "Cuticle IRS", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(73), "Companion layer", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(68,69), "Henle's layer (upper)", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(58,9,51,1), "Basal ORS", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(0,61,64), "WFDC3+ Inner bulge", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(14,49,12), "LQ", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(21,34), "Basal JZ_4", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(15,38,66,28,5,6,63,13), "SHG", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(41,25,29,32,39,22,46,79,11,45,54), "Supr. ORS", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(72,20,36,53), "CST6+ Inner bulge", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(52), "MUCL1+ Inner bulge", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)
extended_lower_hf_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_lower_hf_scvi$seurat_clusters %in% c(40,77,57,10,2,55,24,37,18,16,17,2,42,4,75,3,70,19,76,74,43,31,7,50,56,23,80,8,35,30), "Outer bulge", extended_lower_hf_scvi$inherited_celltype_lvl_5_extended)

DimPlot(extended_lower_hf_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)

saveRDS(extended_lower_hf_scvi,"./Extended_subsets/raw/lower_hf_extended.rds")

# Transfer lvl_5 cell type annotations from lower hair follicle subset to final PSU dataset
common_cells <- intersect(rownames(extended_psu_scvi@meta.data), rownames(extended_lower_hf_scvi@meta.data))

extended_psu_scvi@meta.data[common_cells, "inherited_celltype_lvl_5_extended"] <- extended_lower_hf_scvi@meta.data[common_cells, "inherited_celltype_lvl_5_extended"]

DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_4", label = T, raster = F )
DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F )


# ===============================================
# Subclustering: Sebaceous gland (SG) cluster from PSU
# ===============================================

DefaultAssay(extended_psu_scvi) <- "RNA"
extended_psu_scvi
extended_psu_scvi[["RNA"]] <- JoinLayers(extended_psu_scvi[["RNA"]])


# subset SG
extended_sg <- subset(extended_psu_scvi, subset = seurat_clusters %in% c(39,25,57,3,40,51,11,52,49,58,1))

DimPlot(extended_sg, group.by = "seurat_clusters", label = T)
DimPlot(extended_sg, group.by = "orig_celltype_lvl_3", label = T)

extended_sg$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
extended_sg[["RNA"]] <- split(extended_sg[["RNA"]], f = extended_sg$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
set.seed(42)
extended_sg_scvi <- process_with_scvi(extended_sg)
extended_sg_scvi <- FindClusters(extended_sg_scvi, resolution = 2.7)

DimPlot(extended_sg_scvi, group.by = "seurat_clusters", label = T)

# Visualizations
DimPlot(extended_sg_scvi, group.by = "seurat_clusters", label = TRUE, reduction = "umap")
DimPlot(extended_sg_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(extended_sg_scvi, group.by = "Dataset", label = FALSE)
FeaturePlot(extended_sg_scvi, c("percent_mito", "nFeature_RNA"), reduction = "umap")

extended_sg_scvi <- JoinLayers(extended_sg_scvi)

# ==== Marker gene discovery  ====
Idents(extended_sg_scvi) <- extended_sg_scvi$seurat_clusters
markers <- FindMarkers(extended_sg_scvi, ident.1 = 15, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 50)

# Feature plots for markers
FeaturePlot(extended_sg_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(extended_sg_scvi, rownames(markers)[13:24], raster = FALSE)

# ==== Feature Plots and Cell type annotation ====

# Basal SD
FeaturePlot(extended_sg_scvi, c("KRT15", "PTN", "COL17A1", "C1QTNF12"))
FeaturePlot(extended_sg_scvi, c("KRT5", "SERPINB4", "COL17A1", "POSTN"))

# Suprabasal SD
FeaturePlot(extended_sg_scvi, c("KRT6A", "IVL", "KRT79", "KRT1"))

# SEB-B, SEB-T, SEB-1, SEB-2, SEB-2L
feature_lists <- list(
  SEB_B = c("NNAT", "IL1R2", "MKI67", "TINAGL1"),
  SEB_T = c("WFDC2", "KRT7", "TBC1D4", "MKI67"),
  SEB_1 = c("HAO2", "ACO1", "HSD11B1", "FASN"),
  SEB_2 = c("PNPLA3", "PLIN5", "PPARG", "FADS2"),
  SEB_2L = c("CRAT", "SEC14L6", "DNASE1L2")
)

for (genes in feature_lists) {
  FeaturePlot(extended_sg_scvi, genes, reduction = "umap")
}

# Additional markers
FeaturePlot(extended_sg_scvi, c("CALML5", "FABP4", "CLDN4"), reduction = "umap")

extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(0), "Basal JZ_4", extended_sg_scvi$seurat_clusters)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(11,28,7,9,17,19,20,8,6), "SG progenitors_4", extended_sg_scvi$inherited_celltype_lvl_5_extended)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(24,22,4,3,27,16,1,5,10,26), "SEB-B", extended_sg_scvi$inherited_celltype_lvl_5_extended)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(15), "Prolif. SEB-B", extended_sg_scvi$inherited_celltype_lvl_5_extended)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(14,2,25), "SEB-T", extended_sg_scvi$inherited_celltype_lvl_5_extended)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(18,12,23), "SEB-1", extended_sg_scvi$inherited_celltype_lvl_5_extended)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(21), "SEB-2", extended_sg_scvi$inherited_celltype_lvl_5_extended)
extended_sg_scvi$inherited_celltype_lvl_5_extended <- ifelse(extended_sg_scvi$seurat_clusters %in% c(13), "SEB-2L", extended_sg_scvi$inherited_celltype_lvl_5_extended)

DimPlot(extended_sg_scvi, group.by = "inherited_celltype_lvl_5_extended", label =T, raster = F)

extended_sg_scvi$inherited_celltype_lvl_5_extended %>% table()

saveRDS(extended_sg_scvi,"./Extended_subsets/raw/sg_extended.rds")

# Transfer lvl_5 cell type annotations from lower hair follicle subset to final PSU dataset

common_cells <- intersect(rownames(extended_psu_scvi@meta.data), rownames(extended_sg_scvi@meta.data))

extended_psu_scvi@meta.data[common_cells, "inherited_celltype_lvl_5_extended"] <- extended_sg_scvi@meta.data[common_cells, "inherited_celltype_lvl_5_extended"]

DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_4", label = T, raster = F)
DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_5", label = T, raster = F)

saveRDS(extended_psu_scvi,"./Extended_subsets/raw/psu_extended.rds")


# ===============================================
# Subclustering: Interfollicular epidermis (IFE)
# ===============================================

DimPlot(extended_epid_scvi,  group.by = "seurat_clusters", label = T, raster = F, repel = F)

# subset to IFE
extended_IFE <- subset(
  extended_epid_scvi,
  subset = seurat_clusters %in% c(72,6,61,1,25,45,20,60,70),
  invert = T
)

DimPlot(extended_IFE, group.by  = "inherited_celltype_lvl_3", label = T, raster = F)
DimPlot(extended_IFE, group.by = "seurat_clusters", label = T, raster = F, repel = T)

extended_IFE$Dataset %>% table()

extended_IFE <- FindClusters(extended_IFE, resolution = 1.5)

extended_IFE

# ==== Feature Plots and Cell type annotation ====

FeaturePlot(extended_IFE, c("SERPINB4", "KRT1", "PTN", "C1QTNF12"))
FeaturePlot(extended_IFE, c("KRT15", "KRT1", "MKI67", "FLG"))
FeaturePlot(extended_IFE, c("COL17A1", "KRT14", "POSTN", "KRT5"))
FeaturePlot(extended_IFE, c("DIO2", "WIF1", "FRZB", "SHROOM3"))
FeaturePlot(extended_IFE, c("S100A8", "S100A7", "S100A9", "KRT6A"))
FeaturePlot(extended_IFE, c("PTN", "GATA6", "KRT79", "LGR6"))
FeaturePlot(extended_IFE, c("C1QTNF12", "GATA6", "KRT79", "PTN"))
FeaturePlot(extended_IFE, c("CALB2", "PI3", "KRT79", "IVL"))
FeaturePlot(extended_IFE, c("KLK6", "LCN2", "RHCG", "CDA"))
FeaturePlot(extended_IFE, c("CRAT", "SEC14L6", "PLIN5", "PNPLA3"))
FeaturePlot(extended_IFE, c("NNAT", "IL1R2", "KRT7", "WFDC2"))
FeaturePlot(extended_IFE, c("MKI67", "TOP2A"))
FeaturePlot(extended_IFE, c("LGR5", "RARB", "TBX2", "LHX2"))
FeaturePlot(extended_IFE, c("RXRA", "RARG"))
FeaturePlot(extended_IFE, "KRT1")

DimPlot(extended_IFE, group.by = "seurat_clusters", label = T)

extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(33,49,23,48,22,57,59,18,8,5,28,35,78,73,7,3,67,66,65,40,30), "Basal KC_3", extended_IFE$seurat_clusters)
extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(2), "Prolif. KC_3", extended_IFE$inherited_celltype_lvl_4_extended)
extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(24,79), "Granular KC_3", extended_IFE$inherited_celltype_lvl_4_extended)
extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(43), "Cornified KC_3", extended_IFE$inherited_celltype_lvl_4_extended)
extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(80), "Adipocyte_2", extended_IFE$inherited_celltype_lvl_4_extended)
extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(68,21,26,19,14,39,31,34,32,10,16,17,51,63,9,12,54,41,47,37,64,42,56,69,15,4,52,29,50,11,44,17,36,62,76,0,74,46,55,38,58,77,71,27,75,13), "Spinous KC_3", extended_IFE$inherited_celltype_lvl_4_extended)
extended_IFE$inherited_celltype_lvl_4_extended <- ifelse(extended_IFE$seurat_clusters %in% c(10,53), "LQ", extended_IFE$inherited_celltype_lvl_4_extended)

DimPlot(extended_IFE, group.by = "inherited_celltype_lvl_4_extended", label = T)

saveRDS(extended_IFE, "./Extended_subsets/raw/ife_extended.rds")
