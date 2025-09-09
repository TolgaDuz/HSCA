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
# Load and inspect HSCA extended eccrine sweat glands
# ===============================

extended_esg <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/esg_extended.rds")

extended_esg

DimPlot(extended_esg, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_esg, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_esg, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_esg, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_esg, group.by = "Core")
DimPlot(extended_esg, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_esg) <- "RNA"
extended_esg <- NormalizeData(extended_esg)

Idents(extended_esg) <- extended_esg$inherited_celltype_lvl_4_extended

new_order_ident <- c(
  "Dark cell",
  "Clear cell 1",
  "Clear cell 2",
  "Basal luminal",
  "S100P+ luminal",
  "S100P- luminal",
  "Myoepithelial",
  "CNTNAP2+ Neuron"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_esg))
setdiff(Idents(extended_esg), new_order_ident)

levels(Idents(extended_esg))
levels(extended_esg) <- new_order_ident
levels(Idents(extended_esg))

genes <- unique(c(
  "KRT19", "KRT8", "KRT18",
  "DCD", "SCGB2A1", "SCGB2A2", "LIPH",
  "LCN2", "WFDC2", "ANKRD36C",
  "CHRM3", "NRG3", "CLDN10",
  "S100A2", "KRT5", "KRT14", "CCL2", "CCR4",
  "KRT6A", "KRT16", "KRT77", "S100P",
  "IVL", "SPINK5", "IFI27",
  "MYH11", "ACTG2", "TAGLN", "ACTA2",
  "CNTNAP2", "CNTN5", "CSMD1", "PTPRD"
))

p <- SCpubr::do_DotPlot(
  sample = extended_esg,
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
  filename = "./dotplot_esg.pdf",
  plot = p2,
  device = "pdf",
  width = 14,
  height = 4
)


extended_esg$inherited_celltype_lvl_4_extended %>% table()

# Define a custom color scheme for each level 4 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_esg,
  group.by = "inherited_celltype_lvl_4_extended",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  legend.position = "right",
  colors.use = c(
    "Dark cell"      = "#732191",
    "Clear cell 1"   = "#1F78B4",
    "Clear cell 2"   = "#B2DF8A",
    "Basal luminal"  = "#FDB462",
    "S100P+ luminal" = "#B22222",
    "S100P- luminal" = "#8B4513",
    "Myoepithelial"  = "#d6ccb8",
    "CNTNAP2+ Neuron" = "#FFFF99"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./esg_umap_extended.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)


# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_esg) <- extended_esg$inherited_celltype_lvl_4

# Assign cells not originating from the HSCA core to "Extended"
extended_esg$inherited_celltype_lvl_4[extended_esg$inherited_celltype_lvl_4 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_5' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_4 <- extended_esg$inherited_celltype_lvl_4 %>% unique()
inherited_celltype_lvl_4

extended_esg$inherited_celltype_lvl_4 %>% table()

extended_celltype_lvl_4 <- extended_esg$inherited_celltype_lvl_4_extended %>% unique()
extended_celltype_lvl_4

setdiff(inherited_celltype_lvl_4, extended_celltype_lvl_4)
setdiff(extended_celltype_lvl_4, inherited_celltype_lvl_4)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_4[!inherited_celltype_lvl_4 %in% c(
  "Spinous KC_3",
  "Infund. KC",
  "Basal JZ",
  "Outer Bulge",
  "DES+ SMC",
  "Prolif. KC_3",
  "Granular KC_3",
  "Basal SG",
  "Supr. SG",
  "Schwann cell_2",
  "LC",
  "Melanocyte_2",
  "Basal KC_3",
  "Supr. Infund. KC",
  "Anti-inflammatory Mph",
  "A1",
  "A2",
  "Merkel cell_3",
  "Extended"
)]

celltypes_clean

table(Idents(extended_esg), useNA = "always")

Idents(extended_esg) <- extended_esg$inherited_celltype_lvl_4
Idents(extended_esg) %>% table()
Idents(extended_esg) %>% unique()

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_esg,
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
    "Dark cell"       = "#732191",
    "Clear cell 1"    = "#1F78B4",
    "Clear cell 2"    = "#B2DF8A",
    "Basal luminal"   = "#FDB462",
    "S100P+ luminal"  = "#B22222",
    "S100P- luminal"  = "#8B4513",
    "Myoepithelial"   = "#d6ccb8",
    "CNTNAP2+ Neuron" = "#FFFF99"
  ))


p2

# Save figure X
ggsave(
  filename = "./esg_core_highlighted.png",
  plot = p2,
  width = 15,
  height = 9,
  dpi = 600
)
