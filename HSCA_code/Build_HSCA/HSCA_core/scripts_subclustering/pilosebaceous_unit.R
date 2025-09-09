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
# Subclustering: Pilosebaceous Unit (PSU)
# ===============================================

core_scvi

# UMAP overview to find PSU cells
DimPlot(core_scvi, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(core_scvi, group.by = "orig_celltype_lvl_3", label = TRUE, alpha = 0.2)

FeaturePlot(core_scvi, c("percent_mito", "nFeature_RNA"), reduction = "umap")
FeaturePlot(core_scvi, "KRT6A", raster = F)

# Track clusters
core_scvi$seurat_clusters_full <- object$seurat_clusters
DimPlot(core_scvi, reduction = "umap", group.by = "seurat_clusters_full", label = T, raster = F)

# Subset to PSU cells
core_psu <- subset(
  core_scvi,
  subset = seurat_clusters %in% c(10,17,33,71,40,16,69)
)

core_psu

core_psu$orig_celltype_lvl_3 %>% table()
core_psu$seurat_clusters %>% table()

DimPlot(core_psu, group.by = "seurat_clusters", label = T)
DimPlot(core_psu, group.by = "orig_celltype_lvl_3", label = T)

# Split RNA assay by dataset for scVI processing
core_psu[["RNA"]] <- split(core_psu[["RNA"]], f = core_psu$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
core_psu_scvi <- process_with_scvi(core_psu)
core_psu_scvi <- FindClusters(core_psu_scvi, resolution = 3.9) #3.5

# Inspect
DimPlot(core_psu_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_psu_scvi, group.by = "seurat_clusters_full", label = TRUE)
DimPlot(core_psu_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(core_psu_scvi, group.by = "Dataset", label = TRUE, raster = FALSE)
DimPlot(core_psu_scvi, group.by = "anatomical_region_level3", label = TRUE, raster = FALSE)

FeaturePlot(core_psu_scvi, c("percent_mito", "nFeature_RNA"), raster = FALSE)
core_psu_scvi <- JoinLayers(core_psu_scvi)

# ==== Marker gene discovery  ====
markers <- FindMarkers(core_psu_scvi, ident.1 = 57, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 20)

FeaturePlot(core_psu_scvi, rownames(markers)[1:30], raster = FALSE)

# ==== Feature Plots and Cell type annotation ====

FeaturePlot(core_psu_scvi, c("SERPINB4", "KRT1", "PTN", "C1QTNF12"))
FeaturePlot(core_psu_scvi, c("KRT15", "KRT1", "MKI67", "FLG"))
FeaturePlot(core_psu_scvi, c("COL17A1", "KRT14", "POSTN", "KRT5"))
FeaturePlot(core_psu_scvi, c("DIO2", "WIF1", "FRZB", "SHROOM3"))
FeaturePlot(core_psu_scvi, c("S100A8", "S100A7", "S100A9", "KRT6A"))
FeaturePlot(core_psu_scvi, c("PTN", "GATA6", "KRT79", "LGR6"))
FeaturePlot(core_psu_scvi, c("C1QTNF12", "GATA6", "KRT79", "PTN"))
FeaturePlot(core_psu_scvi, c("CALB2", "PI3", "KRT79", "IVL"))
FeaturePlot(core_psu_scvi, c("KLK6", "LCN2", "RHCG", "CDA"))
FeaturePlot(core_psu_scvi, c("CRAT", "SEC14L6", "PLIN5", "PNPLA3"))
FeaturePlot(core_psu_scvi, c("NNAT", "IL1R2", "KRT7", "WFDC2"))
FeaturePlot(core_psu_scvi, c("MKI67", "TOP2A"))
FeaturePlot(core_psu_scvi, c("LGR5", "RARB", "TBX2", "LHX2"))
FeaturePlot(core_psu_scvi, c("RXRA", "RARG"))
FeaturePlot(core_psu_scvi, c("FST", "LGR6"))

# DimPlots with/without repel
DimPlot(core_psu_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = TRUE)
DimPlot(core_psu_scvi, group.by = "seurat_clusters", reduction = "umap", label = TRUE, raster = FALSE, repel = FALSE)

core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(1), "SG progenitors", core_psu_scvi$seurat_clusters)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(43,31,50), "Supr. SG", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(57,8,13,18), "Basal SG", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(26), "Prolif. basal SG", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(22,41,27), "Basal KC_3", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(20), "Basal SD", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(5,28,17,14,36,37,24,3,51,11,4,59), "Basal JZ", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(60,15), "Supr. JZ & SD", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(54,39), "Supr. Infund. KC", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(core_psu_scvi$seurat_clusters %in% c(30,40,48,7,9,10,56,46,33,21), "Infund. KC", core_psu_scvi$inherited_celltype_lvl_4)
core_psu_scvi$inherited_celltype_lvl_4 <- ifelse(
  core_psu_scvi$seurat_clusters %in% c(
    45,32,52,58,53,23,42,61,34,29,47,44,12,38,19,6,16,35,49,2,0,55,25
  ),
  "lower HF",
  core_psu_scvi$inherited_celltype_lvl_4
)

# Plot annotations
DimPlot(core_psu_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE, raster = FALSE)
DimPlot(core_psu_scvi, group.by = "celltype_lvl_4", label = TRUE, raster = FALSE)

# Save annotated object
saveRDS(core_psu_scvi, "./Core_subsets/raw/psu_core.rds")


# ===============================================
# Subclustering: Lower hairfollicle cluster from PSU
# ===============================================

DefaultAssay(core_psu_scvi) <- "RNA"
core_psu_scvi
core_psu_scvi[["RNA"]] <- JoinLayers(core_psu_scvi[["RNA"]])

core_psu_scvi

# Quick overview
DimPlot(core_psu_scvi, group.by = "seurat_clusters", label = TRUE, raster = FALSE)

# subset to lower hairfollicle
core_lower_hf <- subset(core_psu_scvi, subset = seurat_clusters %in% c(
  45,32,52,58,53,23,42,61,34,29,47,44,12,38,19,6,16,35,49,2,0,55,25
))

# Inspect subset
core_lower_hf
DimPlot(core_lower_hf, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_lower_hf, group.by = "orig_celltype_lvl_3", label = TRUE)
core_lower_hf$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_lower_hf[["RNA"]] <- split(core_lower_hf[["RNA"]], f = core_lower_hf$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
set.seed(42)
core_lower_hf_scvi <- process_with_scvi(core_lower_hf)
core_lower_hf_scvi <- FindClusters(core_lower_hf_scvi, resolution = 5.5)  # 6.5

# Visualize clusters and metadata
DimPlot(core_lower_hf_scvi, group.by = "seurat_clusters", label = TRUE)
DimPlot(core_lower_hf_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_lower_hf_scvi, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(core_lower_hf_scvi, group.by = "anatomical_region_level2", label = TRUE)
DimPlot(core_lower_hf_scvi, group.by = "Dataset", label = FALSE)

core_lower_hf_scvi$inherited_celltype_lvl_3 %>% table()

FeaturePlot(core_lower_hf_scvi, c("percent_mito", "nFeature_RNA"))

core_lower_hf_scvi <- JoinLayers(core_lower_hf_scvi)

# ==== Marker gene discovery ====
markers <- FindMarkers(core_lower_hf_scvi, ident.1 = 40, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 20)

# FeaturePlots of top markers
FeaturePlot(core_lower_hf_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(core_lower_hf_scvi, rownames(markers)[13:24], raster = FALSE)
FeaturePlot(core_lower_hf_scvi, rownames(markers)[21:30], raster = FALSE)

# ==== Feature Plots and Cell type annotation ====
# Huxley's layer
FeaturePlot(core_lower_hf_scvi, c("TGFA", "CDSN", "SPINK5", "CLDN4"))
# Henle's layer (lower)
FeaturePlot(core_lower_hf_scvi, c("SLPI", "RHCG", "KRT71", "CAPN8"))
# Henle's layer (middle)
FeaturePlot(core_lower_hf_scvi, c("CRCT1", "CPA4", "CNFN", "ASPRV1"))
# Henle's layer (upper)
# Caution: comp and selenop also in ORS
FeaturePlot(core_lower_hf_scvi, c("CHI3L1", "KRT75", "COMP", "SELENOP"))
# Late Matrix
FeaturePlot(core_lower_hf_scvi, c("LPAR6", "MYCN", "FABP5", "SCEL"))
# Early Matrix
FeaturePlot(core_lower_hf_scvi, c("HIST1H1B", "TOP2A", "MYCN", "HIST1H1D"))
# Medulla
FeaturePlot(core_lower_hf_scvi, c("KRTAP10-2", "SCYGR4", "KRTAP12-1", "SCYGR8"))
FeaturePlot(core_lower_hf_scvi, c("VSIG8", "KRT86", "KRT83", "KRTAP5-8"))
# Cortex
FeaturePlot(core_lower_hf_scvi, c("SELENBP1", "KRT85", "LY6G6D", "GPNMB")) # LEF1+ Cortex
FeaturePlot(core_lower_hf_scvi, c("LEF1")) # LEF1+ Cortex
FeaturePlot(core_lower_hf_scvi, c("KRT81", "KRT83", "KRT86")) # KRT83+ Cortex
# Melanocytes
FeaturePlot(core_lower_hf_scvi, c("TYRP1", "MLANA", "TYR", "QPCT"))
# DPC
FeaturePlot(core_lower_hf_scvi, c("IGFBP3", "DKK2", "RSPO4", "WIF1"))
# Cuticle
# IRS cuticle
FeaturePlot(core_lower_hf_scvi, c("KRT25", "KRT26", "KRT27", "KRT28", "KRT71", "KRT72", "KRT73"))
# Cuticle hair shaft
FeaturePlot(core_lower_hf_scvi, c("KRT82", "PROCR", "VSIG8", "ACTBL2"))
FeaturePlot(core_lower_hf_scvi, c("KRT32", "SOX21", "CYP26B1", "PPP2R1B", "PADI3"))
FeaturePlot(core_lower_hf_scvi, c("KRT32", "KRT35", "KRT39", "KRT40", "KRT82", "KRT85"))
# CL
FeaturePlot(core_lower_hf_scvi, c("GAL", "GABRP", "SERPINA1", "CALB2"))
FeaturePlot(core_lower_hf_scvi, c("KRT6A", "SH3D21", "NES", "ADAMTS1"))
FeaturePlot(core_lower_hf_scvi, c("VEGFA", "CA9", "MEST", "ADAMTS1"))
# Basal ORS
FeaturePlot(core_lower_hf_scvi, c("LGR5", "COMP", "SELENOP", "ANGPTL7"))
FeaturePlot(core_lower_hf_scvi, c("FGF18", "CRLF1", "ALAD", "RNF152"))
# Bulge
FeaturePlot(core_lower_hf_scvi, c("KRT6A", "KRT16", "S100A2", "SERPINA3"))
FeaturePlot(core_lower_hf_scvi, c("DSG3", "DSP", "PKP1", "FABP5"))
FeaturePlot(core_lower_hf_scvi, c("CST6", "KRT6A"))
FeaturePlot(core_lower_hf_scvi, c("MUCL1"))
DFeaturePlot(core_lower_hf_scvi, c("FABP5"))
FeaturePlot(core_lower_hf_scvi, c("DIO2", "PDK4", "CXCL14"))
FeaturePlot(core_lower_hf_scvi, c("COL6A1"))
FeaturePlot(core_lower_hf_scvi, c("ANGPTL7"))
# lower catagen
FeaturePlot(core_lower_hf_scvi, c("SULT1E1", "NMB", "THBS1", "CLU", "PAPLN"))
# SHG
FeaturePlot(core_lower_hf_scvi, c(
  "KREMEN2", "SLC2A3", "EPCAM", "TNF", "SOX11", "CASC15",
  "SPRY4", "SPX5", "ROBO1", "CUX2", "RARB"
))

DimPlot(core_lower_hf_scvi, group.by = "seurat_clusters", label = T)
DimPlot(core_lower_hf_scvi, group.by = "seurat_clusters_full", label = T)

core_lower_hf_scvi <- JoinLayers(core_lower_hf_scvi)

core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(50,17), "Matrix (Early)", core_lower_hf_scvi$seurat_clusters)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(42,45), "Matrix (Late) & Henle's layer (lower)", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(52,33), "LEF1+ Cortex", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(49), "KRT83+ Cortex", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(62), "Cuticle hair shaft", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(57), "Cuticle IRS", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(35), "Huxley's layer", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(61), "Companion layer", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(55,54), "Henle's layer (upper)", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(20,8,44), "Basal ORS", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(1), "WFDC3+ Inner bulge", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(63, 3, 39, 11, 15, 56, 2, 23, 46, 51), "LQ", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(10), "Lower isthmus", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(21,13,4,31,38), "SHG", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(7,58,36,28,5,59,32), "Supr. ORS", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(19,43,22), "CST6+ Inner bulge", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(53), "MUCL1+ Inner bulge", core_lower_hf_scvi$inherited_celltype_lvl_5)
core_lower_hf_scvi$inherited_celltype_lvl_5 <- ifelse(core_lower_hf_scvi$seurat_clusters %in% c(41,27,0,9,37,24,34,6,16,29,18,25,26,60,40,48,14,47,30,12), "Outer bulge", core_lower_hf_scvi$inherited_celltype_lvl_5)

DimPlot(core_lower_hf_scvi, group.by = "inherited_celltype_lvl_5", label = T, raster = F)

saveRDS(core_lower_hf_scvi,"./Core_subsets/raw/lower_hf_core.rds")

# Transfer lvl_5 cell type annotations from lower hair follicle subset to final PSU dataset
common_cells <- intersect(rownames(core_psu_scvi@meta.data), rownames(core_lower_hf_scvi@meta.data))

core_psu_scvi@meta.data[common_cells, "inherited_celltype_lvl_5"] <- core_lower_hf_scvi@meta.data[common_cells, "inherited_celltype_lvl_5"]

DimPlot(core_psu_scvi, group.by = "inherited_celltype_lvl_4", label = T, raster = F)
DimPlot(core_psu_scvi, group.by = "inherited_celltype_lvl_5", label = T, raster = F)


# ===============================================
# Subclustering: Sebaceous gland (SG) cluster from PSU
# ===============================================

DefaultAssay(core_psu_scvi) <- "RNA"
core_psu_scvi
core_psu_scvi[["RNA"]] <- JoinLayers(core_psu_scvi[["RNA"]])

core_psu_scvi

# UMAP overview to find SG cells
DimPlot(core_psu_scvi, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)

core_psu_scvi <- FindClusters(core_psu_scvi, resolution = 3.3)

# subset SG
core_sg <- subset(core_psu_scvi, subset = seurat_clusters %in% c(50,31,43,13,8,57,18,26,1))

DimPlot(core_sg, group.by = "seurat_clusters", label = T)
DimPlot(core_sg, group.by = "orig_celltype_lvl_3", label = T)

core_sg$Dataset %>% table()

# Split RNA assay by dataset for scVI processing
core_sg[["RNA"]] <- split(core_sg[["RNA"]], f = core_sg$Dataset)

# ==== Re-process immune subset with scVI (Subclustering step) ====
set.seed(42)
core_sg_scvi <- process_with_scvi(core_sg)
core_sg_scvi <- FindClusters(core_sg_scvi, resolution = 2.1)

# Visualizations
DimPlot(core_sg_scvi, group.by = "seurat_clusters", label = TRUE, reduction = "umap")
DimPlot(core_sg_scvi, group.by = "inherited_celltype_lvl_4", label = TRUE)
DimPlot(core_sg_scvi, group.by = "Dataset", label = FALSE)
FeaturePlot(core_sg_scvi, c("percent_mito", "nFeature_RNA"), reduction = "umap")

core_sg_scvi <- JoinLayers(core_sg_scvi)

# ==== Marker gene discovery  ====
Idents(core_sg_scvi) <- core_sg_scvi$seurat_clusters
markers <- FindMarkers(core_sg_scvi, ident.1 = 15, max.cells.per.ident = 500, only.pos = TRUE)
head(markers, 50)

# Feature plots for markers
FeaturePlot(core_sg_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(core_sg_scvi, rownames(markers)[13:24], raster = FALSE)

# ==== PHATE  for the visualization of SG differentiation ====
library(phateR)
library(readr)
library(reticulate)
library(RColorBrewer)
reticulate::use_condaenv("/opt/conda/envs/scvi")
py_discover_config("phate")

oupPhate <- phate(Embeddings(
  core_sg_scvi,
  reduction = "pca.scvi"
),
  seed = 200,
  knn = 30,
  npca = 30
)

oupDR <- oupPhate$embedding
rownames(oupDR) <- colnames(core_sg_scvi)
colnames(oupDR) <- c("PHATE_1","PHATE_2")
core_sg_scvi[["phate1"]] <- CreateDimReducObject(
  embeddings = oupDR, key = "PHATE_",
  assay = DefaultAssay(core_sg_scvi)
)

nClust <- uniqueN(Idents(core_sg_scvi))         # Setup color palette
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)

p1 <- DimPlot(core_sg_scvi, reduction = "umap", label = TRUE,
              label.size = 3, cols = colCls) + coord_fixed()

p2 <- DimPlot(core_sg_scvi, reduction = "phate1", label = TRUE,
              label.size = 3, cols = colCls) + coord_fixed()

p2

# Scale for Seurat FeaturePlot
subset100 <- core_sg_scvi@reductions$phate1@cell.embeddings * 100
core_sg_scvi[["phate"]] <- CreateDimReducObject(embeddings = subset100, key = "PHATE_",
                                                assay = DefaultAssay(core_sg_scvi))

# ==== Feature Plots and Cell type annotation ====

# Basal SD
FeaturePlot(core_sg_scvi, c("KRT15", "PTN", "COL17A1", "C1QTNF12"))
FeaturePlot(core_sg_scvi, c("KRT5", "SERPINB4", "COL17A1", "POSTN"))
FeaturePlot(core_sg_scvi, c("KRT15", "PTN", "COL17A1", "C1QTNF12"), reduction = "phate")

# Suprabasal SD
FeaturePlot(core_sg_scvi, c("KRT6A", "IVL", "KRT79", "KRT1"))
FeaturePlot(core_sg_scvi, c("KRT6A", "IVL", "KRT79", "KRT1"), reduction = "phate")

# SEB-B, SEB-T, SEB-1, SEB-2, SEB-2L
feature_lists <- list(
  SEB_B = c("NNAT", "IL1R2", "MKI67", "TINAGL1"),
  SEB_T = c("WFDC2", "KRT7", "TBC1D4", "MKI67"),
  SEB_1 = c("HAO2", "ACO1", "HSD11B1", "FASN"),
  SEB_2 = c("PNPLA3", "PLIN5", "PPARG", "FADS2"),
  SEB_2L = c("CRAT", "SEC14L6", "DNASE1L2")
)

for (genes in feature_lists) {
  FeaturePlot(core_sg_scvi, genes, reduction = "umap")
  FeaturePlot(core_sg_scvi, genes, reduction = "phate")
}

# Additional markers
FeaturePlot(core_sg_scvi, c("CALML5", "FABP4", "CLDN4"), reduction = "phate")
FeaturePlot(core_sg_scvi, c("CALML5", "FABP4", "CLDN4"), reduction = "umap")

core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(5), "Basal JZ_4", core_sg_scvi$seurat_clusters)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(10,9,16,13,0), "SG progenitors_4", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(4,7,19,1,8,20), "SEB-B", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(11), "Prolif. SEB-B", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(14,3,18), "SEB-T", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(15,2), "LQ", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(6,21), "SEB-1", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(12), "SEB-2", core_sg_scvi$inherited_celltype_lvl_5)
core_sg_scvi$inherited_celltype_lvl_5 <- ifelse(core_sg_scvi$seurat_clusters %in% c(17), "SEB-2L", core_sg_scvi$inherited_celltype_lvl_5)

core_sg_scvi <- JoinLayers(core_sg_scvi)

# ==== Annotate SG subtypes ====
Idents(core_sg_scvi) <- core_sg_scvi$inherited_celltype_lvl_5
markers <- FindMarkers(core_sg_scvi, ident.1 = "SEB-2L", only.pos = TRUE)
head(markers, 20)

# Feature plots for marker genes
FeaturePlot(core_sg_scvi, rownames(markers)[1:12], raster = FALSE)
FeaturePlot(core_sg_scvi, rownames(markers)[11:20], raster = FALSE)
FeaturePlot(core_sg_scvi, rownames(markers)[21:30], raster = FALSE)

# Visualize final annotated clusters
DimPlot(core_sg_scvi, group.by = "inherited_celltype_lvl_5", label = TRUE, raster = FALSE)
DimPlot(core_sg_scvi, group.by = "inherited_celltype_lvl_5", label = TRUE, reduction = "phate")

core_sg_scvi$inherited_celltype_lvl_5 %>% table()

saveRDS(core_sg_scvi,"./Core_subsets/raw/sg_core.rds")

# Transfer lvl_5 cell type annotations from lower hair follicle subset to final PSU dataset
common_cells <- intersect(rownames(core_psu_scvi@meta.data), rownames(core_sg_scvi@meta.data))

core_psu_scvi@meta.data[common_cells, "inherited_celltype_lvl_5"] <- core_sg_scvi@meta.data[common_cells, "inherited_celltype_lvl_5"]

DimPlot(core_psu_scvi, group.by = "inherited_celltype_lvl_4", label = T, raster = F)
DimPlot(core_psu_scvi, group.by = "inherited_celltype_lvl_5", label = T, raster = F)

saveRDS(core_psu_scvi,"./Core_subsets/raw/psu_core.rds")
