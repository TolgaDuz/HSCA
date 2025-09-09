###############################################################################
# Objective:
# This script reconstructs the atlas object in Seurat using the saved data
# structures generated during the transfer learning workflow
# (see vignette: "label_transfer_HSCA.ipynb").
# The scVI embedding obtained from transfer learning is employed for
# dimensionality reduction, clustering, and visualization.
###############################################################################

set.seed(42)

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SCpubr)
library(SeuratWrappers)
library(BPCells)
library(reticulate)

# Configure reticulate to use scvi conda environment
reticulate::use_condaenv(
  "scvi",
  conda = "/opt/conda/condabin/conda",
  required = TRUE
)


# Working directory containing the converted intermediate files.
# These were exported from the original Python atlas object
# (since transfer learning with scVI was performed in Python).
# Key structures (e.g., matrix, metadata, scVI embedding)
# were converted into Seurat-compatible formats to reconstruct
# the atlas object in R.
# Note: due to their large size, these intermediate files are
# not included in the GitHub repository.
setwd("~/HSCA_data/extended/processed_seurat_objects/")

# Load raw 10X Genomics data and metadata
Raw_data <- Read10X(data.dir = "matrix_files")
meta_data <- read.csv("metadata.csv")

# Create Seurat object with raw counts and metadata
atlas <- CreateSeuratObject(counts = Raw_data, meta.data = meta_data)

atlas

DefaultAssay(atlas)

# ==== Load precomputed scVI embedding (from external Python pipeline) ====
scvi_df <- read.csv("scvi_embedding.csv", row.names = 1)

# Ensure ordering of cells in embedding matches Seurat object
scvi_df <- scvi_df[colnames(atlas), ]
scvi_df

# Confirm row/column consistency
head(rownames(scvi_df))
head(colnames(atlas))
stopifnot(all(rownames(scvi_df) == colnames(atlas)))

# Convert embedding DataFrame to matrix
scvi_mat <- as.matrix(scvi_df)
scvi_mat

# Rename columns to 1, 2, …, ncol (Seurat requires proper numeric dimension keys)
colnames(scvi_mat) <- seq_len(ncol(scvi_mat))
colnames(scvi_mat)
rownames(scvi_mat) <- rownames(scvi_df)

# Create Seurat DimReduc object from scVI embedding
scvi_reduction <- CreateDimReducObject(
  embeddings = scvi_mat,
  key = "SCVI_",
  assay = DefaultAssay(atlas)
)

# Add scVI embedding to Seurat atlas object
atlas[["X_scvi_emb"]] <- scvi_reduction


# ================================================
# Downstream analysis: UMAP, neighbors, clustering
# ================================================

atlas
set.seed(8)
atlas <- RunUMAP(
  atlas,
  reduction = "X_scvi_emb",
  reduction.name = "umap",
  return.model = T,
  dims = 1:10
)
atlas <- FindNeighbors(atlas, dims = 1:10, reduction = "X_scvi_emb")
set.seed(8)
atlas <- FindClusters(atlas, resolution = 3.5)
atlas <- NormalizeData(atlas)

# Visualization
FeaturePlot(atlas, "CCER2", raster = F)
FeaturePlot(atlas, "inherited_celltype_lvl_3_transfer_uncert", raster = F)
DimPlot(atlas, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)
DimPlot(atlas, reduction = "umap", group.by = "celltype_lvl_3", label = T, raster = F)
DimPlot(atlas, reduction = "umap", group.by = "Dataset", label = T, raster = F)

table(is.na(atlas$inherited_celltype_lvl_3_transferred_label))

atlas$inherited_celltype_lvl_3_transferred_label %>% table()


# =========================================================================
# Function: inherited_celltypes()
# Iterates over levels 1–5 and cleans annotation columns by replacing empty
# values with appropriate defaults ("Extended" or "Unknown").
# =========================================================================

inherited_celltypes <- function(seurat_obj) {
  for (lvl in 1:5) {
    for (suffix in c("", "_transferred_label")) {
      label_col <- paste0("inherited_celltype_lvl_", lvl, suffix)

      # Skip if column does not exist in metadata
      if (!label_col %in% colnames(seurat_obj@meta.data)) {
        message(paste("Column", label_col, "not found. Skipping..."))
        next
      }

      # Define replacement value depending on column type
      replacement <- ifelse(suffix == "_transferred_label", "Unknown", "Extended")

      # Replace empty values with replacement
      current_values <- seurat_obj@meta.data[[label_col]]
      current_values[current_values == ""] <- replacement
      seurat_obj@meta.data[[label_col]] <- current_values
    }
  }
  return(seurat_obj)
}

atlas <- inherited_celltypes(atlas)

# Visualization of level 3 transferred labels
DimPlot(
  atlas,
  reduction = "umap",
  group.by = "inherited_celltype_lvl_3_transferred_label",
  label = TRUE,
  raster = FALSE
)

# Save raw HSCA extended
saveRDS(atlas, "./HSCA_extended_raw.rds")
