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
# Load and inspect HSCA extended endothelial cells
# ===============================

extended_ec <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/ec_extended.rds")

extended_ec

DimPlot(extended_ec, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_ec, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_ec, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_ec, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_ec, group.by = "Core")
DimPlot(extended_ec, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_ec) <- "RNA"
extended_ec <- NormalizeData(extended_ec)
Idents(extended_ec) <- extended_ec$inherited_celltype_lvl_5_extended

Idents(extended_ec) %>% unique()

new_order_ident <- c(
  "Arterial EC_3",
  "Capillary EC_3",
  "Venous 1 EC_4",
  "Venous 2 EC_4",
  "RSG5+ EC_4",
  "LYVE1+ LEC_4",
  "SCG3+ LEC_4",
  "NEO1+ LEC_4"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_ec))
setdiff(Idents(extended_ec), new_order_ident)

levels(Idents(extended_ec))
levels(extended_ec) <- new_order_ident
levels(Idents(extended_ec))

genes <- unique(c(
  "IGFBP3", "CXCL12", "HEY1",
  "PCSK5", "NEBL", "SEMA3G", "FBLN5",
  "APLN", "RGCC", "H19",
  "FABP4", "CD36", "BTNL9", "RBP7",
  "G0S2", "SELE", "ACKR1", "CSF3", "VCAM1", "IL6",
  "CCL14", "AQP1", "ID1",
  "LDB2", "TLL1",
  "RGS5", "ANPEP",
  "CCL21", "LYVE1", "NRP2", "PTX3",
  "SCG3", "NRXN3", "LYPD6", "RADIL", "HGF", "ADM",
  "NEO1", "SCN3A", "DPP4", "SCG3"
))

p <- SCpubr::do_DotPlot(
  sample = extended_ec,
  features = genes,
  cluster = F,
  dot.scale = 12,
  legend.position = "right",
  font.size = 18,
  legend.length = 10,
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
  filename = "./dotplot_ec.pdf",
  plot = p2,
  device = "pdf",
  width = 18,
  height = 4
)


extended_ec$inherited_celltype_lvl_5_extended %>% table()

# Define a custom color scheme for each level 5 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_ec,
  group.by = "inherited_celltype_lvl_5_extended",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  legend.position = "right",
  colors.use = c(
    "Arterial EC_3"  = "#08306B",
    "Capillary EC_3" = "#1F78B4",
    "Venous 1 EC_4"  = "#B2DF8A",
    "Venous 2 EC_4"  = "#006D2C",
    "RSG5+ EC_4"     = "#B22222",
    "LYVE1+ LEC_4"   = "#8B4513",
    "SCG3+ LEC_4"    = "#FFFF99",
    "NEO1+ LEC_4"    = "#6A3D9A"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./ec_umap_extended.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)

# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_ec) <- extended_ec$inherited_celltype_lvl_5

# Assign cells not originating from the HSCA core to "Extended"
extended_ec$inherited_celltype_lvl_5[extended_ec$inherited_celltype_lvl_5 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_5' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_5 <- extended_ec$inherited_celltype_lvl_5 %>% unique()
inherited_celltype_lvl_5

extended_ec$inherited_celltype_lvl_5 %>% table()

extended_celltype_5 <- extended_ec$inherited_celltype_lvl_5_extended %>% unique()

setdiff(inherited_celltype_lvl_5, extended_celltype_5)
setdiff(extended_celltype_5, inherited_celltype_lvl_5)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_5[!inherited_celltype_lvl_5 %in% c(
  "Extended", "SEB-T", "Supr. Infund. KC_4", "Spinous KC_3", "FRMD3+ SMC", "Anagen DP (C5)_4"
)]

celltypes_clean

table(Idents(extended_ec), useNA = "always")

Idents(extended_ec) %>% table()
Idents(extended_ec) %>% unique()

Idents(extended_ec) <- extended_ec$inherited_celltype_lvl_5

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_ec,
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
p2 <- p +
  scale_color_manual(values = c(
    "Arterial EC_3"   = "#08306B",
    "Capillary EC_3"  = "#1F78B4",
    "Venous 1 EC_4"   = "#B2DF8A",
    "Venous 2 EC_4"   = "#006D2C",
    "RSG5+ EC_4"      = "#B22222",
    "LYVE1+ LEC_4"    = "#8B4513",
    "SCG3+ LEC_4"     = "#FFFF99"
  ))

p2

# Save figure X
ggsave(
  filename = "./ec_core_highlighted.png",
  plot = p2,
  width = 15,
  height = 9,
  dpi = 600
)
