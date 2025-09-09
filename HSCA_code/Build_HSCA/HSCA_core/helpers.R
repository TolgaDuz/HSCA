# ===============================================
# helpers.R
# ===============================================
# Helper functions for HSCA workflow
# ===============================================

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
reticulate::use_condaenv("scvi", conda = "/opt/conda/condabin/conda",required = TRUE)

# ===============================================
# 1. Process a Seurat object with scVI
# ===============================================

# Function to process a Seurat object through feature selection, scaling, PCA,
# and scVI-based layer integration followed by UMAP and neighbor graph calculation
process_with_scvi <- function(
    seurat_obj,
    conda_env_path = "/opt/conda/envs/scvi",
    seed_integrate = 42,
    seed_umap = 8
) {

  # Feature selection, scaling, PCA
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)

  # Integrate with scVI
  set.seed(seed_integrate)
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    method = scVIIntegration,
    assay = "RNA",
    orig.reduction = "pca",
    new.reduction = "pca.scvi",
    conda_env = conda_env_path,
    verbose = TRUE,
    seed = seed_integrate
  )

  # UMAP and neighbors on scVI latent space
  set.seed(seed_umap)
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca.scvi",
    reduction.name = "umap",
    return.model = TRUE,
    dims = 1:30
  )

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "pca.scvi")

  return(seurat_obj)
}

# ===============================================
# 2. Remove low-quality (LQ) cells based on inherited cell type annotations
# ===============================================

remove_LQ_cells <- function(seurat_obj) {
  has_lvl4 <- "inherited_celltype_lvl_4" %in% colnames(seurat_obj[[]])
  has_lvl5 <- "inherited_celltype_lvl_5" %in% colnames(seurat_obj[[]])

  total_cells <- ncol(seurat_obj)
  keep <- rep(TRUE, total_cells)

  if (has_lvl4) {
    keep <- keep & !grepl("LQ", seurat_obj$inherited_celltype_lvl_4)
  }
  if (has_lvl5) {
    keep <- keep & !grepl("LQ", seurat_obj$inherited_celltype_lvl_5)
  }

  removed_cells <- total_cells - sum(keep)
  message(removed_cells, " cells removed.")

  return(seurat_obj[, keep])
}

# ===============================================
# 3. Clean numeric-only labels in inherited_celltype_lvl_5
# ===============================================

clean_lvl5_numeric <- function(object) {
  if ("inherited_celltype_lvl_5" %in% colnames(object@meta.data)) {
    object$inherited_celltype_lvl_5[
      grepl("^\\d+$", object$inherited_celltype_lvl_5)
    ] <- NA
  }
  return(object)
}

# ===============================================
# 4. Inherit lvl4 annotations to lvl5, appending suffix "_4" if not already present
# ===============================================

inherit_lvl4_to_lvl5 <- function(seurat_obj) {
  if (!"inherited_celltype_lvl_5" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$inherited_celltype_lvl_5 <- NA
  }

  # Only for cells where inherited_celltype_lvl_5 is still empty
  empty_lvl5 <- is.na(seurat_obj$inherited_celltype_lvl_5) | seurat_obj$inherited_celltype_lvl_5 == ""

  lvl4 <- seurat_obj$inherited_celltype_lvl_4[empty_lvl5]

  # Check if lvl4 ends with "_<number/letter>"
  ends_with_underscore <- grepl("_[^_]+$", lvl4)

  # If true: use lvl4 directly
  # If false: append "_4"
  new_lvl5 <- ifelse(ends_with_underscore, lvl4, paste0(lvl4, "_4"))

  # Write the result into the corresponding positions in lvl5
  seurat_obj$inherited_celltype_lvl_5[empty_lvl5] <- new_lvl5

  return(seurat_obj)
}
