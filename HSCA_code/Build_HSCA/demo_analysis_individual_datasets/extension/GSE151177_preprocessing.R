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

setwd("~/HSCA_data/extended/GSE151177/")

# ===============================================
# Data read-in and initial Seurat object creation
# ===============================================

accession_list <- c(
  "GSM4567877",
  "GSM4567878",
  "GSM4567879",
  "GSM4567880",
  "GSM4567881",
  "GSM4567882"
)

sample_list <- c(
  "Control01",
  "Control02",
  "Control03",
  "Control04",
  "Control05",
  "Control05F"
)

# Read in data as seperated Seurat Objects# ==== Read in data as separated Seurat Objects ====
for (i in seq_along(accession_list)){
  srr <- accession_list[i]
  samplenr <- sample_list[i]
  name <- paste0(srr, "_", samplenr)
  print(name)

  # Read matrix, features, and barcodes
  cts <- ReadMtx(mtx = paste0(name, "_matrix.mtx.gz"),
                 features = paste0(name, "_features.tsv.gz"),
                 cells = paste0(name, "_barcodes.tsv.gz"))

  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts,  min.features = 0, min.cells = 0))
}

# ==== merge all individual seurat objects into on ====
seurat_obj <- merge(
  GSM4567877_Control01,
  y = c(
    GSM4567878_Control02,
    GSM4567879_Control03,
    GSM4567880_Control04,
    GSM4567881_Control05,
    GSM4567882_Control05F
  ),
  add.cell.ids = c(
    "GSM4567877_Control01",
    "GSM4567878_Control02",
    "GSM4567879_Control03",
    "GSM4567880_Control04",
    "GSM4567881_Control05",
    "GSM4567882_Control05F"
  ),
  project = "Skin_sc"
)

# free storage from unused objects
rm(
  GSM4567877_Control01,
  GSM4567878_Control02,
  GSM4567879_Control03,
  GSM4567880_Control04,
  GSM4567881_Control05,
  GSM4567882_Control05F
)
gc()

seurat_obj <- JoinLayers(seurat_obj)

# ==== Add metadata columns (general info) ====
seurat_obj$sample <- NA
seurat_obj$subject_ID <- NA
seurat_obj$age_years <- NA
seurat_obj$age_range <- NA
seurat_obj$sex <- "Male"
seurat_obj$ethnicity <- NA
seurat_obj$anatomical_region_level1 <- NA
seurat_obj$anatomical_region_level2 <- NA
seurat_obj$anatomical_region_level3 <- NA
seurat_obj$sequencing_platform <- "Illumina NextSeq 500"
seurat_obj$reference_genome_compact <- "GRCh38"
seurat_obj$cell_ranger_version_compact <- "3"
seurat_obj$single_cell_platform <- "10X 3' v2"
seurat_obj$tissue_sampling_type <- "punch biopsies"
seurat_obj$Accession_source <- "GSE151177"
seurat_obj$Dataset <- "Kim_Krüger_2022"
seurat_obj$Condition <- "Healthy"
seurat_obj$Core <- "No"

# ==== Map accession IDs to sample metadata ====
source_patterns <- c(
  "GSM4567877_Control01",
  "GSM4567878_Control02",
  "GSM4567879_Control03",
  "GSM4567880_Control04",
  "GSM4567881_Control05",
  "GSM4567882_Control05F"
)

target_subject_ID <- c(
  "GSE151177_S1",
  "GSE151177_S2",
  "GSE151177_S3",
  "GSE151177_S4",
  "GSE151177_S5",
  "GSE151177_S5"
)
target_sample <- accession_list
target_age <- c("39", "38", "67", "41", "51", "51")
target_age_range <- c("26-40", "26-40", "61-80", "41-60", "41-60", "41-60")

for (i in seq_along(source_patterns)) {
  seurat_obj$sample[grepl(source_patterns[i], colnames(seurat_obj))] <- target_sample[i]
  seurat_obj$subject_ID[grepl(source_patterns[i], colnames(seurat_obj))] <- target_subject_ID[i]
  seurat_obj$age_years[grepl(source_patterns[i], colnames(seurat_obj))] <- target_age[i]
  seurat_obj$age_range[grepl(source_patterns[i], colnames(seurat_obj))] <- target_age_range[i]
}

seurat_obj$sample %>% table()
Idents(seurat_obj) <- seurat_obj$sample

# ===============================================
# Preprocessing and quality control
# ===============================================

# Compute mitochondrial read percentage
seurat_obj[["percent_mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Look at quality criterias before QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(seurat_obj, group.by = "sample", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

seurat_obj_lenient_qc <- subset(seurat_obj, subset = nCount_RNA >= 200)

# ==== scDblFinder: doublet detection ====
# Vignette adapted from: https://biostatsquid.com/scdblfinder-tutorial/

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj_lenient_qc)

# Parallel parameters
bp <- MulticoreParam(3, RNGseed = 1234)

# Run scDblFinder
sce <- scDblFinder(sce, samples = "sample", BPPARAM = bp)

# Inspect results
table(sce$scDblFinder.class)
sce@colData@listData %>% as.data.frame() %>% head()

# Explore results and add to seurat object
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>%
  dplyr::select(starts_with("scDblFinder"))

head(meta_scdblfinder)
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)
table(meta_scdblfinder$scDblFinder.class)

# Add scDblFinder results to Seurat object
seurat_obj_lenient_qc <- AddMetaData(
  object = seurat_obj_lenient_qc,
  metadata = meta_scdblfinder %>% select("scDblFinder.class")
)

table(seurat_obj_lenient_qc$scDblFinder.class)

# QC visualization for doublets vs singlets
VlnPlot(
  seurat_obj_lenient_qc,
  group.by = "sample",
  split.by = "scDblFinder.class",
  features = c("nFeature_RNA", "percent_mito"),
  ncol = 3,
  pt.size = 0
) + theme(legend.position = "right")


doublets_summary <- seurat_obj_lenient_qc@meta.data %>%
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
  file = paste0("./scDblFinder_doublets_summary", seurat_obj_lenient_qc$Accession_source[[1]], "_2.txt"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# Remove doublets and filter by QC metrics
seurat_scdbl <- subset(seurat_obj_lenient_qc, subset = scDblFinder.class == "singlet")

seurat_scdbl

# ==== scDblFinder END ====

# stringent filtering based on mitochondrial content and number of features
seurat_obj_clean <- subset(seurat_scdbl, subset = percent_mito <= 25 & nFeature_RNA >= 200)

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
write_matrix_dir(
  mat = seurat_obj_clean[["RNA"]]$counts,
  dir = "../processed_seurat_objects/on_disc_counts/GSE151177"
)

counts_mat <- open_matrix_dir(
  dir = "../processed_seurat_objects/on_disc_counts/GSE151177"
)

seurat_obj_clean[["RNA"]]$counts <- counts_mat
seurat_obj_clean

saveRDS(seurat_obj_clean, "../processed_seurat_objects/GSE151177.rds")

rm(list = ls(all.names = TRUE))
gc()
