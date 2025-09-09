set.seed(42)

library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(SoupX)
library(scDblFinder)
library(BiocParallel)
library(BPCells)
library(harmony)

setwd("~/HSCA_data/core/HRA000395/")

# ===============================================
# Data read-in and initial Seurat object creation
# ===============================================

accession_list <- c(
  "HRS118996",
  "HRS118997",
  "HRS118998",
  "HRS118999",
  "HRS119000",
  "HRS119001",
  "HRS119002",
  "HRS119003",
  "HRS119004"
)

# ==== Read in data as separated Seurat Objects ====
for (i in seq_along(accession_list)) {
  srr <- accession_list[i]
  name <- srr
  print(name)

  # Read matrix, features, and barcodes
  cts <- ReadMtx(
    mtx = paste0("./", srr, "/outs/filtered_feature_bc_matrix/matrix.mtx.gz"),
    features = paste0("./", srr, "/outs/filtered_feature_bc_matrix/features.tsv.gz"),
    cells = paste0("./", srr, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
  )

  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts,  min.features = 0, min.cells = 0))
}


# ==== merge all individual seurat objects into on ====
seurat_obj <- merge(
  HRS118996,
  y = c(
    HRS118997,
    HRS118998,
    HRS118999,
    HRS119000,
    HRS119001,
    HRS119002,
    HRS119003,
    HRS119004
  ),
  add.cell.ids = c(
    "HRS118996",
    "HRS118997",
    "HRS118998",
    "HRS118999",
    "HRS119000",
    "HRS119001",
    "HRS119002",
    "HRS119003",
    "HRS119004"
  ),
  project = "HSCA_HRA000395"
)

seurat_obj <- JoinLayers(seurat_obj)

# ==== Add metadata columns (general info) ====
seurat_obj$sample <- NA
seurat_obj$subject_ID <- NA
seurat_obj$age_years <- NA
seurat_obj$age_range <- NA
seurat_obj$sex <- "Female"
seurat_obj$ethnicity <- NA
seurat_obj$anatomical_region_level1 <- "Head"
seurat_obj$anatomical_region_level2 <- "Face"
seurat_obj$anatomical_region_level3 <- "Upper eyelid"
seurat_obj$sequencing_platform <- "Illumina NovaSeq 6000"
seurat_obj$reference_genome_compact <- "GRCh38"
seurat_obj$cell_ranger_version_compact <- "7"
seurat_obj$single_cell_platform <- "10X 3' v2"
seurat_obj$tissue_sampling_type <- "surgical"
seurat_obj$Accession_source <- "HRA000395"
seurat_obj$Dataset <- "Zou_Liu_2021"
seurat_obj$Condition <- "Healthy"
seurat_obj$Core <- "Yes"

# ==== Map accession IDs to sample metadata ====
source_patterns <- accession_list
target_subject_ID <- paste0("HRA000395_S", 1:9)
target_sample <- source_patterns
target_age_years <- c("18", "22", "23", "44", "47", "48", "70", "73", "76")
target_age_range <- c("13-25",  "13-25", "13-25", "41-60", "41-60", "41-60", "61-80", "61-80", "61-80")

for (i in seq_along(source_patterns)) {
  seurat_obj$sample[grepl(source_patterns[i], colnames(seurat_obj))] <- target_sample[i]
  seurat_obj$subject_ID[grepl(source_patterns[i], colnames(seurat_obj))] <- target_subject_ID[i]
  seurat_obj$age_years[grepl(source_patterns[i], colnames(seurat_obj))] <- target_age_years[i]
  seurat_obj$age_range[grepl(source_patterns[i], colnames(seurat_obj))] <- target_age_range[i]
}

seurat_obj$sample %>% table()

Idents(seurat_obj) <- seurat_obj$sample


# ===============================================
# Preprocessing and quality control
# ===============================================

# Compute mitochondrial read percentage
seurat_obj[["percent_mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Look at quality criteria before QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(seurat_obj, features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

# Split Seurat object by sample for individual QC / SoupX processing
data_split <- SplitObject(seurat_obj, split.by = "sample")

# ==== SoupX: ambient RNA correction ====

source("~/HSCA_code/Build_HSCA/demo_analysis_individual_datasets/core/soupx_helpers.R") # custom helper functions for SoupX

# Run SoupX processing on each sample
data_split <- process_soupx(data_split)

# Inspect counts before and after correction
sum(data_split[[1]]@assays$original.counts@counts) # total counts before correction
# ratio of adjusted counts
sum(GetAssayData(data_split[[1]], assay = "RNA", slot = "counts")) /
  sum(GetAssayData(data_split[[1]], assay = "original.counts", slot = "counts"))

# Merge back into a single Seurat object
seurat_soupx <- Reduce(
  function(x, y) merge(x, y, project = "Skin_filtered"),
  data_split
) %>%
  JoinLayers()

# Basic filtering: remove very low count cells
min_features <- min(seurat_soupx$nFeature_RNA)
min_counts <- min(seurat_soupx$nCount_RNA)

seurat_soupx_filt <- subset(seurat_soupx, subset = nCount_RNA >= 200)

# Quick checks
seurat_soupx
seurat_soupx_filt
min_features
min_counts

# ==== SoupX END ====

# ==== scDblFinder: doublet detection ====
# Vignette adapted from: https://biostatsquid.com/scdblfinder-tutorial/

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_soupx_filt)

# Parallel parameters
bp <- MulticoreParam(3, RNGseed = 1234)

# Run scDblFinder on each sample
sce <- scDblFinder(sce, samples = "sample", BPPARAM = bp)

# Inspect results
table(sce$scDblFinder.class)
sce@colData@listData %>% as.data.frame() %>% head()
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>%
  dplyr::select(starts_with("scDblFinder"))

head(meta_scdblfinder)
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)
table(meta_scdblfinder$scDblFinder.class)

# Add scDblFinder results to Seurat object
seurat_soupx_filt <- AddMetaData(
  object = seurat_soupx_filt,
  metadata = meta_scdblfinder %>% select("scDblFinder.class")
)

table(seurat_soupx_filt$scDblFinder.class)

# QC visualization for doublets vs singlets
VlnPlot(
  seurat_soupx_filt,
  group.by = "sample",
  split.by = "scDblFinder.class",
  features = c("nFeature_RNA", "percent_mito"),
  ncol = 3,
  pt.size = 0
) + theme(legend.position = "right")

# Summarize doublets for reporting
doublets_summary <- seurat_soupx_filt@meta.data %>%
  group_by(sample, scDblFinder.class) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  as.data.frame() %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(countT = sum(total_count)) %>%
  group_by(scDblFinder.class, .add = TRUE) %>%
  mutate(percent = paste0(round(100 * total_count/countT, 2), "%")) %>%
  dplyr::select(-countT)

doublets_summary

write.table(
  doublets_summary,
  file = paste0("./scDblFinder_doublets_summary", seurat_soupx_filt$Accession_source[[1]], "_2.txt"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# Remove doublets and filter by QC metrics
seurat_soupx_scdbl <- subset(seurat_soupx_filt, subset = scDblFinder.class == "singlet")

seurat_soupx_scdbl

# ==== scDblFinder END ====

# stringent filtering based on mitochondrial content and number of features
seurat_obj_clean <- subset(seurat_soupx_scdbl, subset = percent_mito <= 25 & nFeature_RNA >= 200)

seurat_obj_clean
# QC visualization post-filter
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(
  seurat_obj_clean,
  group.by = "sample",
  features = feats,
  pt.size = 0.1,
  ncol = 3
) + NoLegend()

seurat_obj_clean$genes <- rownames(seurat_obj_clean)


# ===============================================
# On-disk matrix storage
# ===============================================

# Save and reload on-disc counts
write_matrix_dir(
  mat = seurat_obj_clean[["RNA"]]$counts,
  dir = "../processed_seurat_objects/on_disc_counts/HRA000395"
)

counts_mat <- open_matrix_dir(
  dir = "../processed_seurat_objects/on_disc_counts/HRA000395"
)

# Replace counts in Seurat object with on-disk matrix
seurat_obj_clean[["RNA"]]$counts <- counts_mat
seurat_obj_clean

# ===============================================
# Post-QC data processing and downstream analysis
# ===============================================

# ==== Standard Seurat workflow: normalization, HVG selection, scaling, PCA ====
set.seed(42)
seurat_obj_clean <- NormalizeData(seurat_obj_clean)
seurat_obj_clean <- FindVariableFeatures(seurat_obj_clean)
seurat_obj_clean <- ScaleData(seurat_obj_clean)
seurat_obj_clean <- RunPCA(seurat_obj_clean)
ElbowPlot(seurat_obj_clean)

# ==== Neighborhood graph, clustering, and UMAP embedding ====
seurat_obj_clean <- FindNeighbors(seurat_obj_clean, dims = 1:30)
seurat_obj_clean <- FindClusters(seurat_obj_clean)
seurat_obj_clean <- RunUMAP(seurat_obj_clean, dims = 1:30)

DimPlot(seurat_obj_clean, group.by = "sample", label = T)
DimPlot(seurat_obj_clean, group.by = "seurat_clusters", label = T)
FeaturePlot(seurat_obj_clean, c("percent_mito", "nFeature_RNA"))
FeaturePlot(seurat_obj_clean, c("KRT1", "KRT15"))

# Batch Correction
set.seed(8)
seurat_obj_har <- RunHarmony(
  seurat_obj_clean,
  group.by.vars = c("sample"),
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 50,
  project.dim = FALSE
)

# UMAP on Harmony embeddings
seurat_obj_har <- RunUMAP(
  seurat_obj_har,
  reduction = "harmony",
  dims = 1:30,
  n.components = 2
)

# Neighborhood graph and clustering on Harmony embeddings
seurat_obj_har <- FindNeighbors(
  seurat_obj_har,
  dims = 1:30,
  reduction = "harmony"
)

seurat_obj_har <- FindClusters(seurat_obj_har, resolution = 1.2)

DimPlot(seurat_obj_har, group.by = "seurat_clusters", label = T)
DimPlot(seurat_obj_har, group.by = "sample")

# === DGEA ====
mark <- FindMarkers(seurat_obj_har, ident.1 = 28, only.pos = T)

mark$gene <- rownames(mark)
head(mark, 10)

# Feature plots for top markers
FeaturePlot(seurat_obj_har, rownames(mark)[1:12])
FeaturePlot(seurat_obj_har, c("percent_mito", "nFeature_RNA"))

# Basal KC
FeaturePlot(seurat_obj_har, c("COL17A1", "KRT15"))
# Prolif. KC
FeaturePlot(seurat_obj_har, c("MKI67", "TOP2A"))
# Spinous KC
FeaturePlot(seurat_obj_har, c("KRT1", "KRT10"))
# Infundibulum
FeaturePlot(seurat_obj_har, c("KRT6A", "S100A8"))
# Isthmus and JZ
FeaturePlot(seurat_obj_har, c("PTN", "CST6", "C1QTNF12", "GATA6")) #FAM132A
# Bulge
FeaturePlot(seurat_obj_har, c("LGR5", "DIO2"))
FeaturePlot(seurat_obj_har, c("CTNND2", "FABP5"))
# Bulb
FeaturePlot(seurat_obj_har_filt, c("KRT23", "KRT71"))
FeaturePlot(seurat_obj_har_filt, c("TCHH", "KRT81"))
FeaturePlot(seurat_obj_har_filt, c("COMP", "ANGPTL7"))
# Basal SG
FeaturePlot(seurat_obj_har, c("NNAT", "IL1R2", "KRT7"))
# SG
FeaturePlot(seurat_obj_har, c("FADS2", "PPARG"))

# ===============================================
# Subclustering: PSU and IFE clusters
# ===============================================

subset_pilo_IFE <- subset(
  seurat_obj_har,
  subset = seurat_clusters %in% c(12,7,9,17,3,1,5,2,10,0,4,8,13,26,6,14,15)
)

DimPlot(subset_pilo_IFE, label = T)

set.seed(42)
subset_pilo_IFE <- FindVariableFeatures(subset_pilo_IFE)
subset_pilo_IFE <- ScaleData(subset_pilo_IFE)
subset_pilo_IFE <- RunPCA(subset_pilo_IFE)
subset_pilo_IFE <- FindNeighbors(subset_pilo_IFE, dims = 1:30)
subset_pilo_IFE <- FindClusters(subset_pilo_IFE)
subset_pilo_IFE <- RunUMAP(subset_pilo_IFE, dims = 1:30)

DimPlot(subset_pilo_IFE, group.by = "seurat_clusters", label = TRUE)
DimPlot(subset_pilo_IFE, group.by = "sample", label = TRUE)

# Batch correction
library(harmony)
set.seed(8)
subset_pilo_IFE_har <- RunHarmony(
  subset_pilo_IFE,
  group.by.vars = c("sample"),
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 50,
  project.dim = FALSE
)

subset_pilo_IFE_har <- RunUMAP(
  subset_pilo_IFE_har,
  reduction = "harmony",
  dims = 1:30,
  n.components = 2
)

subset_pilo_IFE_har <- FindNeighbors(
  subset_pilo_IFE_har,
  dims = 1:30,
  reduction = "harmony"
)
subset_pilo_IFE_har <- FindClusters(subset_pilo_IFE_har, resolution = 1.4)

DimPlot(subset_pilo_IFE_har, group.by = "seurat_clusters", label = TRUE)
DimPlot(subset_pilo_IFE_har, group.by = "sample")

# ==== Marker discovery and QC ====
mark <- FindMarkers(subset_pilo_IFE_har, ident.1 = 17, only.pos = TRUE) %>%
  tibble::rownames_to_column("gene")

head(mark, 10)
FeaturePlot(subset_pilo_IFE_har, mark$gene[1:8])
FeaturePlot(subset_pilo_IFE_har, c("percent_mito", "nFeature_RNA"))

# Identify and remove low-quality (LQ) cluster
FeaturePlot(subset_pilo_IFE_har, c("percent_mito", "nFeature_RNA"))
Idents(subset_pilo_IFE_har) %>% table()
cells_to_remove <- WhichCells(subset_pilo_IFE_har, idents = c(17))
length(cells_to_remove)

subset_pilo_IFE_har_filt <- subset(subset_pilo_IFE_har, cells = setdiff(Cells(subset_pilo_IFE_har), cells_to_remove))
# remove cells also in full seurat object
seurat_obj_har_filt <- subset(seurat_obj_har, cells = setdiff(Cells(seurat_obj_har), cells_to_remove))

DimPlot(subset_pilo_IFE_har_filt, label = T)
DimPlot(seurat_obj_har, label = T)
DimPlot(seurat_obj_har_filt, label = T)

# ==== Cell type annotation (orig_celltype_lvl_3) ====

# Basal KC
FeaturePlot(subset_pilo_IFE_har_filt, c("COL17A1", "KRT15"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(0,1,6),
  "Basal KC",
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters
)

# Prolif. KC
FeaturePlot(subset_pilo_IFE_har_filt, c("MKI67", "TOP2A"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(18,16,12,13),
  "Prolif. KC",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

# Spinous KC
FeaturePlot(subset_pilo_IFE_har_filt, c("KRT1", "KRT10"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(7,21,3,10,9,5,2,11),
  "Spinous KC",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

# Granular KC
FeaturePlot(subset_pilo_IFE_har_filt, c("KRT2", "FLG"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(14),
  "Granular KC",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

# Cornified KC (markers shown but not assigned)
FeaturePlot(subset_pilo_IFE_har_filt, c("LCE1A", "LCE1B"))

# Infundibulum
FeaturePlot(subset_pilo_IFE_har_filt, c("KRT6A", "S100A8"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(8),
  "Infundibulum",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

# Isthmus and JZ
FeaturePlot(subset_pilo_IFE_har_filt, c("PTN", "CST6", "C1QTNF12", "GATA6")) #FAM132A
FeaturePlot(subset_pilo_IFE_har_filt, c("PTN", "KRT79")) #FAM132A
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(15,4),
  "Isthmus",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

# Bulge
FeaturePlot(subset_pilo_IFE_har_filt, c("LGR5", "DIO2"))
FeaturePlot(subset_pilo_IFE_har_filt, c("CTNND2", "FABP5"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(19,20),
  "Bulge",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

# Basal SG
FeaturePlot(subset_pilo_IFE_har_filt, c("NNAT", "IL1R2", "KRT7")) # not present

# SG
FeaturePlot(subset_pilo_IFE_har_filt, c("FADS2", "PPARG"))
subset_pilo_IFE_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_pilo_IFE_har_filt@meta.data$seurat_clusters %in% c(22),
  "SG",
  subset_pilo_IFE_har_filt@meta.data$orig_celltype_lvl_3
)

DimPlot(subset_pilo_IFE_har_filt, group.by = "orig_celltype_lvl_3", label = T)

# Map labels back to original object
labels <- subset_pilo_IFE_har_filt$orig_celltype_lvl_3
seurat_obj_har_filt$orig_celltype_lvl_3 <- NA
seurat_obj_har_filt$orig_celltype_lvl_3[names(labels)] <- labels

DimPlot(seurat_obj_har_filt, group.by = "orig_celltype_lvl_3", label = T)
DimPlot(seurat_obj_har_filt, group.by = "seurat_clusters", label = T)

# ===============================================
# Subclustering: Immune cells
# ===============================================

subset_immune <- subset(
  seurat_obj_har_filt,
  subset = seurat_clusters %in% c(18,19,29)
)

DimPlot(subset_immune, label = T)

set.seed(42)
subset_immune <- FindVariableFeatures(subset_immune)
subset_immune <- ScaleData(subset_immune)
subset_immune <- RunPCA(subset_immune)
ElbowPlot(subset_immune)

# Batch correction
library(harmony)
set.seed(8)
subset_immune_har <- RunHarmony(
  subset_immune,
  group.by.vars = c("sample"),
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 50,
  project.dim = FALSE
)

subset_immune_har<- RunUMAP(
  subset_immune_har,
  reduction = "harmony",
  dims = 1:30,
  n.components = 2
)

subset_immune_har <- FindNeighbors(
  subset_immune_har,
  dims = 1:30,
  reduction = "harmony"
)

subset_immune_har <- FindClusters(subset_immune_har, resolution = 1.1)
DimPlot(subset_immune_har, group.by = "seurat_clusters", label = T)
DimPlot(subset_immune_har, group.by = "sample")

# Remove low-quality cluster
FeaturePlot(subset_immune_har, c("percent_mito", "nFeature_RNA"))
Idents(subset_immune_har) %>% table()
cells_to_remove <- WhichCells(subset_immune_har, idents = c(5))
length(cells_to_remove)

subset_immune_har_filt <- subset(subset_immune_har, cells = setdiff(Cells(subset_immune_har), cells_to_remove))
# remove cells also in full seurat object
seurat_obj_har_filt <- subset(seurat_obj_har_filt, cells = setdiff(Cells(seurat_obj_har_filt), cells_to_remove))

DimPlot(subset_immune_har_filt, label = T)
DimPlot(seurat_obj_har_filt, label = T)

# ==== Marker discovery and QC ====

mark <- FindMarkers(subset_immune_harm_filt, ident.1 = 12, only.pos = TRUE) %>%
  tibble::rownames_to_column("gene")

head(mark, 10)
FeaturePlot(subset_immune_harm_filt, mark$gene[1:4])
FeaturePlot(seurat_obj_har_filt, mark$gene[1:12])


# ==== Cell type annotation (orig_celltype_lvl_3) ====

# Dendritic cells (DC)
FeaturePlot(subset_immune_har_filt,  c("CD207", "FCGBP", "CD1C", "FCER1A", "CLEC10A", "CD1E"))
subset_immune_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_immune_har_filt@meta.data$seurat_clusters %in% c(9,10,3,12,15,17,1),
  "DC",
  subset_immune_har_filt@meta.data$seurat_clusters
)

# T cell
# CD8+ T cells ()
FeaturePlot(subset_immune_har_filt, c("CD8A", "CD3E"))
# cytoxic t cells are also all cd8+
FeaturePlot(subset_immune_har_filt, c("GZMK", "GZMA","CD8A", "FXYD2", "CD8B", "ZNF683", "NKG7", "CCL5", "IFNG", "CST7", "GZMH"))
FeaturePlot(subset_immune_har_filt, c("CCR6", "NTRK2", "IL26", "TOX", "CAMK4", "BATF", "SPOCK2", "DUSP4"))
subset_immune_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_immune_har_filt@meta.data$seurat_clusters %in% c(6,0,2,8,7),
  "T cell",
  subset_immune_har_filt@meta.data$orig_celltype_lvl_3
)

# Neutrophil
FeaturePlot(subset_immune_har_filt, c("S100A8", "S100A9", "FCGR3B", "CMTM2", "AQP9"))
subset_immune_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_immune_har_filt@meta.data$seurat_clusters %in% c(13),
  "Neutrophil",
  subset_immune_har_filt@meta.data$orig_celltype_lvl_3
)

# NK cell
FeaturePlot(subset_immune_har_filt, c("GNLY", "GZMB", "KLRD1", "PRF1", "FGFBP2", "KLRF1", "FCGR3A", "SPON2", "CLICL3"))
subset_immune_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_immune_har_filt@meta.data$seurat_clusters %in% c(16,14),
  "NK cell",
  subset_immune_har_filt@meta.data$orig_celltype_lvl_3
)

# Macrophage (Mph)
FeaturePlot(subset_immune_har_filt, c("CD68", "CTSD", "ACP5",  "SPP1"))
FeaturePlot(subset_immune_har_filt, c("C1QA", "SLENOP", "FOLR2", "RNASE1", "C1QC", "DAB2",
  "F13A1", "C1QB", "CD14", "STAB1", "LGMN"))
subset_immune_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_immune_har_filt@meta.data$seurat_clusters %in% c(4),
  "Mph",
  subset_immune_har_filt@meta.data$orig_celltype_lvl_3
)

# Classical monocyte
FeaturePlot(subset_immune_har_filt, c("S100A12", "FCN1", "RNASE2"))
subset_immune_har_filt$orig_celltype_lvl_3 <- ifelse(
  subset_immune_har_filt@meta.data$seurat_clusters %in% c(11),
  "Monocyte",
  subset_immune_har_filt@meta.data$orig_celltype_lvl_3
)

# ==== Map immune labels back to full object ====
labels <- subset_immune_harm_filt$orig_celltype_lvl_3
seurat_obj_har_filt$orig_celltype_lvl_3[names(labels)] <- labels

# Inspect
DimPlot(subset_immune_harm_filt, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "seurat_clusters", label = TRUE)


# ===============================================
# Subclustering: Endothelial cells
# ===============================================
# Subset endothelial cells
subset_ec <- subset(seurat_obj_har_filt, subset = seurat_clusters %in% c(21))

DimPlot(subset_ec)

set.seed(42)
subset_ec <- FindVariableFeatures(subset_ec)
subset_ec <- ScaleData(subset_ec)
subset_ec <- RunPCA(subset_ec)
ElbowPlot(subset_ec)

library(harmony)
set.seed(8)
subset_ec_har <- RunHarmony(
  subset_ec,
  group.by.vars = c("sample"),
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 50,
  project.dim = FALSE
)

subset_ec_har <- RunUMAP(
  subset_ec_har,
  reduction = "harmony",
  dims = 1:30,
  n.components = 2
)

subset_ec_har <- FindNeighbors(
  subset_ec_har,
  dims = 1:30,
  reduction = "harmony"
)

subset_ec_har <- FindClusters(subset_ec_har, resolution = 1.2)
DimPlot(subset_ec_har, group.by = "seurat_clusters", label = T)
DimPlot(subset_ec_har, group.by = "sample")

FeaturePlot(subset_ec_har, c("percent_mito", "nFeature_RNA"))

# ==== Marker discovery and QC ====
mark <- FindMarkers(subset_ec_harm, ident.1 = 8, only.pos = TRUE) %>%
  tibble::rownames_to_column("gene")
head(mark, 10)
FeaturePlot(subset_ec_harm, mark$gene[1:4])


# ==== Cell type annotation (orig_celltype_lvl_3) ====
# Arterial EC
FeaturePlot(subset_ec_har, c("EFNB2", "DLL4", "NOTCH4", "IGFBP3"))
FeaturePlot(subset_ec_har, c("ARL15", "CXCL12", "HEY1"))
subset_ec_har$orig_celltype_lvl_3 <- ifelse(
  subset_ec_har@meta.data$seurat_clusters %in% c(8,6),
  "Arterial EC",
  subset_ec_har@meta.data$seurat_clusters
)

# capillary EC
FeaturePlot(subset_ec_har, c("APLN", "VWA1", "RGCC", "PLVAP", "MT1M"))
subset_ec_har$orig_celltype_lvl_3 <- ifelse(
  subset_ec_har@meta.data$seurat_clusters %in% c(4,7),
  "Capillary EC",
  subset_ec_har@meta.data$orig_celltype_lvl_3
)

# venule 1
FeaturePlot(subset_ec_har, c("CCL14", "AQP1", "ID1", "ACKR1"))
FeaturePlot(subset_ec_har, c("TPD52L1", "ZNF385D", "OLFM1", "ACKR1"))
subset_ec_har$orig_celltype_lvl_3 <- ifelse(
  subset_ec_har@meta.data$seurat_clusters %in% c(2,0,5,3,10,1,9),
  "Venous EC",
  subset_ec_har@meta.data$orig_celltype_lvl_3
)

# ==== Map endothelial labels back to full object ====
labels <- subset_ec_harm$orig_celltype_lvl_3
seurat_obj_har_filt$orig_celltype_lvl_3[names(labels)] <- labels

# Inspect
DimPlot(subset_ec_harm, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "orig_celltype_lvl_3", label = TRUE)

# ===============================================
# Subclustering: Fibroblast cells
# ===============================================
# Subset fibroblast clusters
subset_fib <- subset(seurat_obj_har_filt, subset = seurat_clusters %in% c(11, 25))
DimPlot(subset_fib, label = TRUE)

set.seed(42)
subset_fib <- FindVariableFeatures(subset_fib)
subset_fib <- ScaleData(subset_fib)
subset_fib <- RunPCA(subset_fib)
ElbowPlot(subset_fib)

library(harmony)
set.seed(8)
subset_fib_har <- RunHarmony(
  subset_fib,
  group.by.vars = c("sample"),
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 50,
  project.dim = FALSE
)

subset_fib_har <- RunUMAP(
  subset_fib_har,
  reduction = "harmony",
  dims = 1:30,
  n.components = 2
)

subset_fib_har <- FindNeighbors(
  subset_fib_har,
  dims = 1:30,
  reduction = "harmony"
)

subset_fib_har <- FindClusters(subset_fib_har, resolution = 0.4)
DimPlot(subset_fib_har, reduction = "umap", group.by = "seurat_clusters", label = T)
DimPlot(subset_fib_har, group.by = "sample")

FeaturePlot(subset_fib_har, c("percent_mito", "nFeature_RNA"))

# ==== Marker discovery and QC ====
mark <- FindMarkers(subset_fib_harm, ident.1 = 8, only.pos = TRUE) %>%
  tibble::rownames_to_column("gene")
head(mark, 10)
FeaturePlot(subset_fib_harm, mark$gene[1:4])


# ==== Cell type annotation (orig_celltype_lvl_3) ====
# A
FeaturePlot(subset_fib_har, c("SFRP2", "ELN", "MMP2", "QPCT"))
subset_fib_har$orig_celltype_lvl_3 <- ifelse(
  subset_fib_har@meta.data$seurat_clusters %in% c(1,0,4),
  "Fibro A",
  subset_fib_har@meta.data$seurat_clusters
)

# B
FeaturePlot(subset_fib_har, c("APOE", "C7", "CYGB", "IGFBP7"))
subset_fib_har$orig_celltype_lvl_3 <- ifelse(
  subset_fib_har@meta.data$seurat_clusters %in% c(2),
  "Fibro B",
  subset_fib_har@meta.data$orig_celltype_lvl_3
)

# C
FeaturePlot(subset_fib_har, c("SFRP1", "TNMD", "DKK3", "TNN"))
subset_fib_har$orig_celltype_lvl_3 <- ifelse(
  subset_fib_har@meta.data$seurat_clusters %in% c(3,5),
  "Fibro C",
  subset_fib_har@meta.data$orig_celltype_lvl_3
)

# C1 (DS) (markers shown but not assigned)
FeaturePlot(subset_fib_har, c("COL11A1", "DPEP1", "TNMD", "WFDC1"))
# C2 (DP) (markers shown but not assigned)
FeaturePlot(subset_fib_har, c("COCH", "CRABP1", "FIBIN", "RSPO4"))
FeaturePlot(subset_fib_har, c("GRIK1"))

# ==== Map fibroblast labels back to full object ====
labels <- subset_fib_harm$orig_celltype_lvl_3
seurat_obj_har_filt$orig_celltype_lvl_3[names(labels)] <- labels

# Inspect
DimPlot(subset_fib_harm, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "seurat_clusters", label = TRUE)


# ===============================================
# Subclustering: Eccrine Sweat Gland (ESG) cells
# ===============================================
# Subset ESG clusters
subset_esg <- subset(seurat_obj_har_filt, subset = seurat_clusters %in% c(20,27))

DimPlot(subset_esg)

set.seed(42)
subset_esg <- FindVariableFeatures(subset_esg)
subset_esg <- ScaleData(subset_esg)
subset_esg <- RunPCA(subset_esg)
ElbowPlot(subset_esg)

# Batch correction
library(harmony)
set.seed(8)
subset_esg_har <- RunHarmony(
  subset_esg,
  group.by.vars = c("sample"),
  reduction = "pca",
  reduction.save = "harmony",
  max_iter = 50,
  project.dim = FALSE
)

subset_esg_har <- RunUMAP(
  subset_esg_har,
  reduction = "harmony",
  dims = 1:30,
  n.components = 2
)

subset_esg_har <- FindNeighbors(
  subset_esg_har,
  dims = 1:30,
  reduction = "harmony"
)

subset_esg_har <- FindClusters(subset_esg_har, resolution = 0.4)
DimPlot(subset_esg_har, group.by = "seurat_clusters", label = T)
DimPlot(subset_esg_har, group.by = "sample")

FeaturePlot(subset_esg_har, c("percent_mito", "nFeature_RNA"))

# ==== Marker discovery and QC ====
mark <- FindMarkers(subset_esg_harm, ident.1 = 8, only.pos = TRUE) %>%
  tibble::rownames_to_column("gene")
head(mark, 10)
FeaturePlot(subset_esg_harm, mark$gene[1:4])


# ==== Cell type annotation (orig_celltype_lvl_3) ====
# ======================
# ESG consists of coil and duct parts
# Coil: dark, clear, and myoepithelial cells
# Duct: luminal and basal cells

# Coil cells
# Dark cells
FeaturePlot(subset_esg_harm, c("DCD", "SCGB2A1", "SCGB2A2", "LIPH"))
# Clear cells
FeaturePlot(subset_esg_harm, c("LCN2"))
FeaturePlot(subset_esg_harm, c("S100A1", "KRT19", "KRT8", "KRT18")) # Secretory coil cells
subset_esg_harm$orig_celltype_lvl_3 <- ifelse(
  subset_esg_harm$seurat_clusters %in% c(1, 6),
  "Coil",
  subset_esg_harm$seurat_clusters
)

# Duct cells
FeaturePlot(subset_esg_harm, c("S100A2", "KRT5", "KRT14", "CCL2", "CCR4", "KRT6A", "KRT16", "KRT77"))
FeaturePlot(subset_esg_harm, c("IVL", "SPINK5", "IFI27"))
subset_esg_harm$orig_celltype_lvl_3 <- ifelse(
  subset_esg_harm$seurat_clusters %in% c(7, 4, 2, 0, 3, 5),
  "Duct",
  subset_esg_harm$orig_celltype_lvl_3
)

# ==== Map ESG labels back to full object ====
labels <- subset_esg_harm$orig_celltype_lvl_3
seurat_obj_har_filt$orig_celltype_lvl_3[names(labels)] <- labels

# Inspect
DimPlot(subset_esg_harm, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "orig_celltype_lvl_3", label = TRUE)
DimPlot(seurat_obj_har_filt, group.by = "seurat_clusters", label = TRUE)


# ==== Other cell type annotations ====
# Lymphatic EC
FeaturePlot(seurat_obj_har_filt,  c("CCL21", "LYVE1"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(23),
  "Lymphatic EC",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Mast cell
FeaturePlot(seurat_obj_har_filt, c("TPSAB1"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(24),
  "Mast cell",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Skeletal muscle
FeaturePlot(seurat_obj_har_filt, c("COX6A2", "KLHL41", "MYL1", "CSRP3"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(30),
  "Skeletal muscle",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# SMC
FeaturePlot(seurat_obj_har_filt, c("MYH9", "TAGLN", "TINAGL1", "RERG"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(22),
  "SMC",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Melanocyte
FeaturePlot(seurat_obj_har_filt, c("MLANA", "PMEL"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(16),
  "Melanocyte",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Other
FeaturePlot(seurat_obj_har_filt, c("MEG3", "TAGLN", "TINAGL1", "RERG"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(28),
  "Other_HRA000395",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Inspect
DimPlot(seurat_obj_har_filt, group.by = "seurat_clusters", label = T)
DimPlot(seurat_obj_har_filt,  group.by = "orig_celltype_lvl_3", label = T)

saveRDS(seurat_obj_har_filt, "../processed_seurat_objects/HRA000395.rds")
