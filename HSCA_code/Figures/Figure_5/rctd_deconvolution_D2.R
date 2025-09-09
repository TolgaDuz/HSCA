set.seed(42)
library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(hdf5r)
library(ggplot2)
library(gridExtra)
library(pals)
library(data.table)
library(tidyr)
library(spacexr)
library(pheatmap)
library(viridis)

setwd("~/HSCA_code/Figures/Figure_5/")

# read spatial dataset as a query
visiumHD_lower_hf <- readRDS("/path/to/figshare/visiumHD_lower_hairfollicle.rds")

# read scRNA-seq dataset as a reference
ref <- readRDS("./bulge_core.rds")

# Remove LQ cells
ref <- subset(ref, subset = inherited_celltype_lvl_5 == "LQ", invert = TRUE)
# Deconvolution on tissue section D2
visiumHD_lower_hf_D2 <- subset(visiumHD.lower.hf, subset = Sample == "D2")

DimPlot(ref, group.by = "inherited_celltype_lvl_5")

### ==== RCTD deconvolution (adapted from Seurat vignette) ====
### This analysis follows and adapts the Seurat Visium HD vignette:
### https://satijalab.org/seurat/articles/visiumhd_analysis_vignette
### We use RCTD to map single-cell references onto Visium HD spots.

Idents(ref) <- "inherited_celltype_lvl_5"
counts <- ref[["RNA"]]$counts
celltype <- as.factor(ref$inherited_celltype_lvl_5)
nUMI <- ref$nCount_RNA

# create the RCTD reference object
reference <- Reference(counts, celltype, nUMI)

counts_hd <- visiumHD_lower_hf_D2[["Spatial.008um"]]$counts
visiumHD_lower_hf_D2_cells_hd <- colnames(visiumHD_lower_hf_D2[["Spatial.008um"]])
coords <- GetTissueCoordinates(visiumHD_lower_hf_D2)[visiumHD_lower_hf_D2_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
visiumHD_lower_hf_D2 <- AddMetaData(visiumHD_lower_hf_D2, metadata = RCTD@results$results_df)
visiumHD_lower_hf_D2 <- AddMetaData(visiumHD_lower_hf_D2, metadata = RCTD@results$weights)

# cells that are part of the reject class
table(visiumHD_lower_hf_D2@meta.data$first_type)

# predicted cells that are part of the reject class
table(visiumHD_lower_hf_D2@meta.data$first_type[visiumHD_lower_hf_D2@meta.data$spot_class == "reject"])

# predicted cells that are part of the doublet uncertain class
table(visiumHD_lower_hf_D2@meta.data$first_type[visiumHD_lower_hf_D2@meta.data$spot_class == "doublet_uncertain"])

# predicted cells that are part of the singlet class
table(visiumHD_lower_hf_D2@meta.data$first_type[visiumHD_lower_hf_D2@meta.data$spot_class == "singlet"])

# predicted cells that are part of the doublet class
table(visiumHD_lower_hf_D2@meta.data$first_type[visiumHD_lower_hf_D2@meta.data$spot_class == "doublet_certain"])

rownames(visiumHD_lower_hf_D2)

# 1. extract celltype
celltypes <- unique(visiumHD_lower_hf_D2@meta.data$first_type)

# 2. name of prediction columns
score_cols <- paste0(celltypes)

# 3. Check if the cols exist in metadata
score_cols <- score_cols[score_cols %in% colnames(visiumHD_lower_hf_D2@meta.data)]
score_cols

# 4. extract pred values and format as a matrix
pred_matrix <- as.matrix(visiumHD_lower_hf_D2@meta.data[, score_cols])
pred_matrix
rownames(pred_matrix) <- colnames(visiumHD_lower_hf_D2)  # ensure correct order

# 5. create new assay and save
pred_assay <- CreateAssayObject(data = t(pred_matrix))
visiumHD_lower_hf_D2[["predictions"]] <- pred_assay

pred_matrix <- visiumHD_lower_hf_D2@assays[["predictions"]]@data

trans_mat <- as.data.frame(t(pred_matrix))
trans_mat

head(rownames(trans_mat))
head(rownames(visiumHD_lower_hf_D2@meta.data))

merged_df <- cbind(visiumHD_lower_hf_D2@meta.data, trans_mat)

rownames(visiumHD_lower_hf_D2@meta.data)
rownames(merged_df)

dim(visiumHD_lower_hf_D2@meta.data)
dim(merged_df)

visiumHD_lower_hf_D2@meta.data <- merged_df

DefaultAssay(visiumHD_lower_hf_D2)
DefaultAssay(visiumHD_lower_hf_D2) <- "predictions"

rownames(visiumHD_lower_hf_D2@assays[["predictions"]]@data)

### ==== RCTD vignette end ==== ###

# === START Create HEATMAP ====
data_matrix <- visiumHD_lower_hf_D2@assays[["predictions"]]@data
cell_types <- visiumHD_lower_hf_D2$celltype

rownames(data_matrix)

cell_ids <- colnames(data_matrix)
cell_ids

cell_ids_by_type <- split(cell_ids, cell_types[cell_ids])
cell_ids_by_type
length(cell_ids_by_type)

sum_values <- matrix(NA, nrow = nrow(data_matrix), ncol = length(unique(cell_types)))
colnames(sum_values) <- unique(cell_types)

for (cell_type in unique(cell_types)) {
  selected_columns <- cell_ids_by_type[[cell_type]]

  sum_values[, cell_type] <- rowSums(data_matrix[, selected_columns, drop = FALSE])
}

rownames(sum_values) <- rownames(data_matrix)

pheatmap(
  sum_values,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row"

)

rownames(sum_values)
colnames(sum_values)

row_order <- c(
  "Outer bulge",
  "MUCL1+ Inner bulge",
  "CST6+ Inner bulge",
  "Basal ORS",
  "Supr. ORS",
  "Companion layer",
  "Henle's layer (upper)",
  "Matrix (Late) & Henle's layer (lower)",
  "Huxley's layer",
  "Cuticle IRS",
  "Cuticle hair shaft",
  "LEF1+ Cortex",
  "KRT83+ Cortex",
  "Matrix (Early)"
)

setdiff(row_order, rownames(sum_values))
visiumHD_lower_hf_D2$celltype %>% table()

col_order <- c(
  "Outer bulge",
  "MUCL1+ Inner bulge",
  "CST6+ Inner bulge",
  "Basal ORS",
  "Supr. ORS",
  "Companion layer",
  "Henle's layer (upper)",
  "Henle's layer (lower)",
  "Huxley's layer",
  "Cuticle",
  "LEF1+ Cortex",
  "KRT83+ Cortex",
  "Matrix (Late)",
  "Matrix (Early)"
)

sum_values_sorted <- sum_values[row_order, col_order, drop = FALSE]

p1 <-  pheatmap::pheatmap(
  sum_values_sorted,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row"

)

p1

# save as PDF
ggsave(
  filename = "./rctd_deconv_D2_zscore.pdf",
  plot = p1,
  device = "pdf",
  width = 8,
  height = 6
)

#### Min max scale
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

data <- sum_values_sorted

scaled_data <- t(apply(data, 1, min_max_scale))

scaled_data_df <- as.data.frame(scaled_data)
colnames(scaled_data_df) <- colnames(data)
rownames(scaled_data_df) <- rownames(data)

# Create Figure 5b
p <- pheatmap::pheatmap(
  scaled_data_df,
  color = viridis::viridis(256),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 14,
  angle_col = "90"
)

# save Figure 5b
ggsave(
  filename = "./rctd_deconv_D2_minmax.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 6
)
