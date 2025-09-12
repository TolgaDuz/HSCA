# setting working dir
setwd("~/HSCA_data/extended/processed_seurat_objects/")

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
  PRJNA754272 = "../Extension/PRJNA754272.rds",
  GSE130973 = "../Extension/GSE130973.rds",
  GSE132802_five_prime = "../Extension/GSE132802_five_prime.rds",
  GSE179633 = "../Extension/GSE179633.rds",
  GSE176415 = "../Extension/GSE176415.rds",
  GSE182208 = "../Extension/GSE182208.rds",
  GSE173205 = "../Extension/GSE173205.rds",
  GSE153760 = "../Extension/GSE153760.rds",
  PRJNA797897 = "../Extension/PRJNA797897.rds",
  GSE151177 = "../Extension/GSE151177.rds",
  GSE175990 = "../Extension/GSE175990.rds",
  GSE181316 = "../Extension/GSE181316.rds",
  GSE125422 = "../Extension/GSE125422.rds",
  GSE147482 = "../Extension/GSE147482.rds",
  GSE202352 = "../Extension/GSE202352.rds",
  GSE162183 = "../Extension/GSE162183.rds",
  GSE165816 = "../Extension/GSE165816.rds",
  GSE241132 = "../Extension/GSE241132.rds",
  GSE265972 = "../Extension/GSE265972.rds",
  EGAS00001002927_2 = "../Extension/EGAS00001002927_2.rds",
  EMTAB13084_1 = "../Extension/EMTAB13084_1.rds"
)

# Initialize lists for loaded objects and non-integer datasets
seurat_objects <- list()
non_integer_studies <- list()

check_counts <- FALSE # Set TRUE to check for non-integer counts

# ======================
# Load Seurat objects
# ======================

# Function to make cell names unique
make_unique_cellnames <- function(seurat_obj, prefix) {
  colnames(seurat_obj) <- paste0(prefix, "_", colnames(seurat_obj))
  seurat_obj
}

for (name in names(rds_files)) {

  seurat_obj <- readRDS(rds_files[[name]])

  print(paste0("../Extension/", name))

  counts_mat <- seurat_obj[["RNA"]]$counts

  if (check_counts) {
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

  seurat_objects[[name]] <- make_unique_cellnames(seurat_obj, name)
}

seurat_objects

# Print the names of the loaded Seurat objects
print(names(seurat_objects))

# ======================
# Merge datasets
# ======================

merged_seurat <- merge(
  seurat_objects$PRJNA754272,
  y = c(
    seurat_objects$GSE130973,
    seurat_objects$GSE132802_five_prime,
    seurat_objects$GSE179633,
    seurat_objects$GSE176415,
    seurat_objects$GSE182208,
    seurat_objects$GSE173205,
    seurat_objects$GSE153760,
    seurat_objects$PRJNA797897,
    seurat_objects$GSE151177,
    seurat_objects$GSE175990,
    seurat_objects$GSE181316,
    seurat_objects$GSE125422,
    seurat_objects$GSE147482,
    seurat_objects$GSE202352,
    seurat_objects$GSE162183,
    seurat_objects$GSE165816,
    seurat_objects$GSE241132,
    seurat_objects$GSE265972,
    seurat_objects$EGAS00001002927_2,
    seurat_objects$EMTAB13084_1
  )
)

merged_seurat

merged_seurat <- JoinLayers(merged_seurat)

merged_seurat

# ======================
# Inspect metadata
# ======================

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

cols_to_keep <- !grepl("^scDblFinder", colnames(merged_seurat@meta.data))

merged_seurat@meta.data <- merged_seurat@meta.data[, cols_to_keep]

merged_seurat

# ======================
# Remove genes with incorrect Excel-style formatting
# ======================

# Some gene names in the original count matrices were converted to dates by Excel (e.g., "1-Mar" instead of "MARCH1").
# This block identifies such genes using a regular expression and removes them from the Seurat object.

merged_seurat$genes <- rownames(merged_seurat[["RNA"]]$counts)

genes <- data.frame(genes = rownames(merged_seurat))

"1-Mar" %in% rownames(merged_seurat[["RNA"]]$counts)
"MARCH1" %in% rownames(merged_seurat[["RNA"]]$counts)

# 1. Get the gene names
all_genes <- rownames(merged_seurat)

# 2. Create a regular expression for Excel date errors like "7-Mar", "3-Sep", etc.
month_abbr <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# 3. Find genes that match this pattern
date_pattern <- paste0("^\\d{1,2}-(", paste(month_abbr, collapse = "|"), ")")
date_pattern
is_date <- grepl(date_pattern, all_genes, ignore.case = FALSE)

# 4. Show affected genes (optional)
cat("Removed genes:\n")
print(all_genes[is_date])

# 5. Filter the object
merged_seurat <- subset(merged_seurat, features = all_genes[!is_date])
merged_seurat

saveRDS(
  object = merged_seurat,
  file = "./HSCA_only_extended_datasets_raw.rds"
)
