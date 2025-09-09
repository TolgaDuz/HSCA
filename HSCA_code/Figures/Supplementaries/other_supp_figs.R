set.seed(42)

# load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SCpubr)
library(SeuratWrappers)
library(BPCells)
library(viridis)

setwd("~/HSCA_code/Figures/Figure_7")

# ===============================
# Load and inspect HSCA extended other lineages
# ===============================

extended_other_cells <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/other_cells_extended.rds")

extended_other_cells

DimPlot(extended_other_cells, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_other_cells, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_other_cells, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_other_cells, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_other_cells, group.by = "Core")
DimPlot(extended_other_cells, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_other_cells) <- "RNA"
extended_other_cells <- NormalizeData(extended_other_cells)

Idents(extended_other_cells) <- extended_other_cells$inherited_celltype_lvl_2_extended

Idents(extended_other_cells) %>% unique()

new_order_ident <- c(
  "Melanocyte",
  "Neuron",
  "Neural progenitor",
  "Schwann cell",
  "Schwann cell (Ex vivo)",
  "Chondrocyte",
  "Erythrocyte"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_other_cells))
setdiff(Idents(extended_other_cells), new_order_ident)

levels(Idents(extended_other_cells))
levels(extended_other_cells) <- new_order_ident
levels(Idents(extended_other_cells))

genes <- unique(c(
  "MLANA", "QPCT", "TYRP1", "PMEL",
  "NRXN1", "STARD13", "XKR4", "NTM",
  "RBP1", "STMN2", "HOXB5", "DCX", "BASP1", "PCSK1N",
  "PRX", "GLDN", "MPZ", "MLIP",
  "CRLF1", "NGFR", "GAP43", "HLA-DRB1", "HLA-DRA", "S100B", "NES", "JUN", "POU3F1", "SOX10", "MPZ", "GAP43", "CDH19",
  "ACAN", "COL2A1",
  "ALAS2", "HBA1", "HBB", "HBA2"
))

p <- SCpubr::do_DotPlot(
  sample = extended_other_cells,
  features = genes,
  cluster = F,
  dot.scale = 9,
  legend.position = "right",
  font.size = 18,
  legend.length = 5,
  min.cutoff = 0,
  zscore.data = T,
  flip = F
)


# Keep only data points where percentage expression is at least 5%
p$data <- p$data %>% filter(P.Exp >= 5)

# Set circle border thickness in the plot to 0.2
p[["layers"]][[1]][["geom"]][["default_aes"]][["stroke"]] <- 0.2

# Everything that is not upregulated (Z-score ≤ 0) is hidden
# Set negative Z-scores to NA
p$data$Avg.Exp[p$data$Avg.Exp <= 0] <- NA

# Also set P.Exp (dot size) to NA if the corresponding expression is NA
p$data$P.Exp[is.na(p$data$Avg.Exp)] <- NA

p

# Determine the upper limit of Avg.Exp for color scaling
upper_limit <- max(p$data$Avg.Exp, na.rm = TRUE)
upper_limit

# Create Supp. Fig. X
p2 <- p + ggplot2::scale_fill_gradientn(
  colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  limits = c(0, upper_limit),
  oob = scales::squish
)

p2

# Save Supp. Fig. X
ggsave(
  filename = "./dotplot_other.pdf",
  plot = p2,
  device = "pdf",
  width = 15,
  height = 4
)

extended_other_cells$inherited_celltype_lvl_2_extended %>% table()

# Define a custom color scheme for each level 5 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_other_cells,
  group.by = "inherited_celltype_lvl_2_extended",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  legend.position = "right",
  colors.use = c(
    "Melanocyte"              = "#d6ccb8",
    "Neuron"                  = "#F18A85",
    "Neural progenitor"       = "#B2DF8A",
    "Schwann cell"            = "#ffff7f",
    "Schwann cell (Ex vivo)"  = "#A6CEE3",
    "Chondrocyte"             = "#8B4513",
    "Erythrocyte"             = "#E31A1C"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./other_umap_extended.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)

# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_other_cells) <- extended_other_cells$inherited_celltype_lvl_2

# Assign cells not originating from the HSCA core to "Extended"
extended_other_cells$inherited_celltype_lvl_2[extended_other_cells$inherited_celltype_lvl_2 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_2' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_2 <- extended_other_cells$inherited_celltype_lvl_2 %>% unique()
inherited_celltype_lvl_2

extended_other_cells$inherited_celltype_lvl_2 %>% table()

extended_celltype_2 <- extended_other_cells$inherited_celltype_lvl_2_extended %>% unique()

setdiff(inherited_celltype_lvl_2, extended_celltype_2)
setdiff(extended_celltype_2, inherited_celltype_lvl_2)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_2[!inherited_celltype_lvl_2 %in% c(
  "Eccrine sweat Gland",
  "Muscle",
  "Fibroblast",
  "Interfollicular Epidermis",
  "Lymphoid",
  "Pilosebaceous Unit",
  "Extended"
)]

celltypes_clean

table(Idents(extended_other_cells), useNA = "always")

Idents(extended_other_cells) %>% table()
Idents(extended_other_cells) %>% unique()

Idents(extended_other_cells) <- extended_other_cells$inherited_celltype_lvl_2

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_other_cells,
  idents.keep = celltypes_clean,
  label = F,
  legend.position = "right",
  border.size = 3,
  pt.size = 0.2,
  na.value = "#363636"
)

p

# Highlight core immune cells on top of the extended dataset using distinct colors
# Create Supplementary Figure X
p2 <- p + scale_color_manual(values = c(
  "Melanocyte"              = "#d6ccb8",
  "Neuron"                  = "#F18A85",
  "Neural progenitor"       = "#B2DF8A",
  "Schwann cell"            = "#ffff7f",
  "Schwann cell (Ex vivo)"  = "#A6CEE3",
  "Chondrocyte"             = "#8B4513",
  "Erythrocyte"             = "#E31A1C"
))

p2

# Save figure X
ggsave(
  filename = "./other_core_highlighted.png",
  plot = p2,
  width = 15,       # in inches
  height = 9,
  dpi = 600
)
