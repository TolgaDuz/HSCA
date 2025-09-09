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

setwd("~/HSCA_data/core/EGAS00001002927/")

# ===============================================
# Data read-in and initial Seurat object creation
# ===============================================

accession_list <- c(
  "EGAF00002191947",
  "EGAF00002191954",
  "EGAF00002191956",
  "EGAF00002191945",
  "EGAF00002191949",
  "EGAF00002191953"
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
  EGAF00002191947,
  y = c(
    EGAF00002191954,
    EGAF00002191956,
    EGAF00002191945,
    EGAF00002191949,
    EGAF00002191953
  ),
  add.cell.ids = c(
    "EGAF00002191947",
    "EGAF00002191954",
    "EGAF00002191956",
    "EGAF00002191945",
    "EGAF00002191949",
    "EGAF00002191953"
  ),
  project = "HSCA_EGAS00001002927"
)

seurat_obj <- JoinLayers(seurat_obj)

# ==== Add metadata columns (general info) ====
seurat_obj$sample <- NA
seurat_obj$subject_ID <- NA
seurat_obj$age_years <- NA
seurat_obj$age_range <- NA
seurat_obj$sex <- NA
seurat_obj$ethnicity <- NA
seurat_obj$anatomical_region_level1 <- NA
seurat_obj$anatomical_region_level2 <- NA
seurat_obj$anatomical_region_level3 <- NA
seurat_obj$sequencing_platform <- NA
seurat_obj$reference_genome_compact <- "GRCh38"
seurat_obj$cell_ranger_version_compact <- "7"
seurat_obj$single_cell_platform <- "10X 3' v2"
seurat_obj$tissue_sampling_type <- "surgical"
seurat_obj$Accession_source <- "EGAS00001002927"
seurat_obj$Dataset <- "Cheng_Cho_2018_1"
seurat_obj$Condition <- "Healthy"
seurat_obj$Core <- "Yes"

# ==== Map accession IDs to sample metadata ====
source_patterns <- accession_list
target_subject_ID <- paste0("EGAS00001002927_S", 1:6)
target_sample <- source_patterns
target_sex <- c("Female", "Female", "Female", "Female", "Male", "Female")
target_region_1 <- c("Trunk", "Trunk", "Trunk", "Head", "Head", "Head")
target_region_2 <- c("Breast", "Abdomen", "Breast", "Scalp", "Scalp", "Scalp")
target_region_3 <- rep(NA, 6)

# Map patterns to sample IDs
for (i in seq_along(source_patterns)) {
  seurat_obj$sample[grepl(source_patterns[i], colnames(seurat_obj))] <- target_sample[i]
  seurat_obj$subject_ID[grepl(source_patterns[i], colnames(seurat_obj))] <- target_subject_ID[i]
  seurat_obj$sex[grepl(source_patterns[i], colnames(seurat_obj))] <- target_sex[i]
  seurat_obj$anatomical_region_level1[grepl(source_patterns[i], colnames(seurat_obj))] <- target_region_1[i]
  seurat_obj$anatomical_region_level2[grepl(source_patterns[i], colnames(seurat_obj))] <- target_region_2[i]
  seurat_obj$anatomical_region_level3[grepl(source_patterns[i], colnames(seurat_obj))] <- target_region_3[i]
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
  dir = "../processed_seurat_objects/on_disc_counts/EGAS00001002927"
)

counts_mat <- open_matrix_dir(
  dir = "../processed_seurat_objects/on_disc_counts/EGAS00001002927"
)

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
seurat_obj_clean <- FindNeighbors(seurat_obj_clean, dims = 1:30)
seurat_obj_clean <- FindClusters(seurat_obj_clean)
seurat_obj_clean <- RunUMAP(seurat_obj_clean, dims = 1:30)

# ==== Neighborhood graph, clustering, and UMAP embedding ====
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
seurat_obj_har <- FindClusters(seurat_obj_har, resolution = 1.8)

DimPlot(seurat_obj_har, reduction = "umap", group.by = "seurat_clusters", label = T)
DimPlot(seurat_obj_har, group.by = "sample")

# === DGEA ====
mark <- FindMarkers(seurat_obj_har, ident.1 = 28, only.pos = T)

mark$gene <- rownames(mark)
head(mark, 10)

# Feature plots for top markers
FeaturePlot(seurat_obj_har, rownames(mark)[1:4])
FeaturePlot(seurat_obj_har, rownames(mark)[5:8])

FeaturePlot(seurat_obj_har, c("percent_mito", "nFeature_RNA"))

seurat_obj_har_filt <- subset(
  seurat_obj_har,
  subset = seurat_clusters %in% c(17, 8, 23),
  invert = T
)

DimPlot(seurat_obj_har_filt, group.by = "seurat_clusters", label = T)

# Initialize cell type column
seurat_obj_har_filt$orig_celltype_lvl_3 <- NA

# ===============================================
# Assign cell types and visualize markers
# ===============================================

# Basal KC
FeaturePlot(seurat_obj_har_filt, c("COL17A1", "KRT15"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(3,0,2,15,22),
  "Basal KC",
  seurat_obj_har_filt@meta.data$seurat_clusters
)

# Prolif. KC
FeaturePlot(seurat_obj_har_filt, c("MKI67", "TOP2A"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(11,20,21),
  "Prolif. KC",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Spinous KC
FeaturePlot(seurat_obj_har_filt, c("KRT1", "KRT10"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(5,4,1,6,19,10),
  "Spinous KC",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Granular KC
FeaturePlot(seurat_obj_har_filt, c("KRT2", "FLG"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(27),
  "Granular KC",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Cornified KC
FeaturePlot(seurat_obj_har_filt, c("LCE1A", "LCE1B"))

# Infundibulum
FeaturePlot(seurat_obj_har_filt, c("KRT6A", "S100A8", "KRT79", "SERPINB3"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(12,28),
  "Infundibulum",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Isthmus and JZ
FeaturePlot(seurat_obj_har_filt, c("PTN", "CST6", "C1QTNF12", "GATA6")) #FAM132A
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(14,7),
  "Isthmus",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Bulge
FeaturePlot(seurat_obj_har_filt, c("PDK4", "DIO2"))
FeaturePlot(seurat_obj_har_filt, c("CTNND2", "FABP5"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(9),
  "Bulge",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Bulb
FeaturePlot(seurat_obj_har_filt, c("KRT23", "KRT71"))
FeaturePlot(seurat_obj_har_filt, c("TCHH", "KRT81"))
FeaturePlot(seurat_obj_har_filt, c("COMP", "ANGPTL7"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(16,25),
  "Bulb",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Basal SG
FeaturePlot(seurat_obj_har_filt, c("NNAT", "IL1R2", "KRT7"))

# SG
FeaturePlot(seurat_obj_har_filt, c("FADS2", "PPARG"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(24),
  "SG",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# DC
FeaturePlot(seurat_obj_har_filt,  c("CD207", "FCGBP"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(30,26),
  "DC",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# T cell
FeaturePlot(seurat_obj_har_filt, c("PTPRC", "IL32", "SRGN", "ARHGAP15"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(29),
  "T cell",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Melanocyte
FeaturePlot(seurat_obj_har_filt, c("MLANA", "PMEL"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(18,13),
  "Melanocyte",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Melanocyte
FeaturePlot(seurat_obj_har_filt, c("KRT6A", "KRT16", "KRT77"))
seurat_obj_har_filt$orig_celltype_lvl_3 <- ifelse(
  seurat_obj_har_filt@meta.data$seurat_clusters %in% c(16),
  "Duct",
  seurat_obj_har_filt@meta.data$orig_celltype_lvl_3
)

# Inspect
DimPlot(seurat_obj_har_filt, group.by = "orig_celltype_lvl_3", label = T)

saveRDS(seurat_obj_har_filt, "../processed_seurat_objects/EGAS00001002927.rds")
