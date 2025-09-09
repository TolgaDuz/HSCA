# setting working dir
setwd("~/HSCA_data/core/processed_seurat_objects/")

set.seed(42)

# load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(BPCells)

# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)

# ======================
# Define input files
# ======================

rds_files <- c(
  HRA000395 = "./HRA000395.rds",
  GSE138669_2021 = "./GSE138669_2021.rds",
  GSE138669_2020 = "./GSE138669_2020.rds",
  GSE132802_three_prime = "./GSE132802_three_prime.rds",
  GSE275491 = "./GSE275491.rds",
  GSE186476 = "./GSE186476.rds",
  GSE147424 = "./GSE147424.rds",
  GSE129611 = "./GSE129611.rds",
  GSE159929 = "./GSE159929.rds",
  EGAS00001002927 = "./EGAS00001002927.rds",
  GSE191067 = "./GSE191067.rds",
  GSE274955 = "./GSE274955.rds",
  EMTAB13084_2 = "./EMTAB13084_2.rds"
)

# Initialize lists for loaded objects and non-integer datasets
seurat_objects <- list()
non_integer_studies <- list()

check_counts <- FALSE # Set TRUE to check for non-integer counts

# ======================
# Load Seurat objects
# ======================

for (name in names(rds_files)) {

  seurat_obj <- readRDS(rds_files[[name]])

  print(paste0("./", name))

  if (check_counts) {
    counts_mat <- seurat_obj[["RNA"]]$counts
    counts_mat <- as.matrix(counts_mat)

    test <- rowSums(counts_mat %% 1 == 0)
    nonint_genes <- rownames(counts_mat)[test != ncol(counts_mat)]

    if (length(nonint_genes) > 0) {
      print("WARNING: THIS DATASET HAS NON-INTEGER VALUES!!!")
      non_integer_studies[[name]] <- nonint_genes
      next
    } else {
      print("Counts are integers.")
    }
  }

  seurat_objects[[name]] <- seurat_obj
}

seurat_objects

# Print the names of the loaded Seurat objects
print(names(seurat_objects))

# ======================
# Merge datasets
# ======================

merged_seurat <- merge(
  seurat_objects$HRA000395,
  y = c(
    seurat_objects$GSE138669_2021,
    seurat_objects$GSE138669_2020,
    seurat_objects$GSE132802_three_prime,
    seurat_objects$GSE275491,
    seurat_objects$GSE186476,
    seurat_objects$GSE147424,
    seurat_objects$GSE129611,
    seurat_objects$GSE159929,
    seurat_objects$EGAS00001002927,
    seurat_objects$GSE191067,
    seurat_objects$GSE274955,
    seurat_objects$EMTAB13084_2
  )
)

merged_seurat

# Free memory
rm(seurat_objects, seurat_obj)
gc()

merged_seurat <- JoinLayers(merged_seurat)

merged_seurat


# ======================
# Inspect metadata
# ======================

merged_seurat$orig_celltype_lvl_3 %>% table()
merged_seurat$sample %>% table()
merged_seurat$sample %>% unique() %>% length()
merged_seurat$subject_ID %>% table()
merged_seurat$subject_ID %>% unique() %>% length()
merged_seurat$age_years %>% table()
merged_seurat$ethnicity %>% table()
merged_seurat$anatomical_region_level1 %>% table()
merged_seurat$anatomical_region_level2 %>% table()
merged_seurat$anatomical_region_level3 %>% table()
merged_seurat$sequencing_platform %>% table()
merged_seurat$reference_genome_compact %>% table()
merged_seurat$cell_ranger_version_compact %>% table()
merged_seurat$single_cell_platform %>% table()
merged_seurat$tissue_sampling_type %>% table()
merged_seurat$Accession_source %>% table()
merged_seurat$Dataset %>% table()
merged_seurat$Condition %>% table()
merged_seurat$Core %>% table()

merged_seurat

# Ensure counts are sparse matrix
merged_seurat[["RNA"]]$counts <- as(merged_seurat[["RNA"]]$counts, Class = "dgCMatrix")

# Remove genes with zero counts across all cells
merged_seurat <- subset(
  merged_seurat,
  features = rownames(merged_seurat)[Matrix::rowSums(merged_seurat[["RNA"]]@layers$counts) > 0]
)

merged_seurat
gc()

saveRDS(
  object = merged_seurat,
  file = "./HSCA_core_raw.rds"
)
