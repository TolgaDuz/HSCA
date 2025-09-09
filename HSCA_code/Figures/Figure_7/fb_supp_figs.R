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
# Load and inspect HSCA extended fibroblasts
# ===============================

extended_fb <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/fb_extended.rds")

DimPlot(extended_fb, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_fb, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_fb, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_fb, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_fb, group.by = "Core")
DimPlot(extended_fb, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_fb) <- "RNA"
extended_fb <- NormalizeData(extended_fb)

Idents(extended_fb) <- extended_fb$inherited_celltype_lvl_4_extended


new_order_ident <- c(
  "A1",
  "A2",
  "A3",
  "A4",
  "B1",
  "B2",
  "B3",
  "B4",
  "Dermal sheath (C1)",
  "Outer bulge DP (C2)",
  "C3",
  "Anagen DP (C5)",
  "D1",
  "D2",
  "RAMP1+ Fibro (E1)",
  "CAF1",
  "CAF2",
  "Merkel cell_3"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_fb))
setdiff(Idents(extended_fb), new_order_ident)

levels(Idents(extended_fb))
levels(extended_fb) <- new_order_ident
levels(Idents(extended_fb))

genes <- unique(c(
  "IGFBP6", "PI116", "SLPI",
  "WISP2", "SEMA3B", "LGR5",
  "APCDD1", "COL18A1", "COMP", "NKD2",
  "HSPB3", "COL6A5",
  "RGCC", "SGCA", "WIF1",
  "SOSTDC1", "CORIN",
  "FBN1", "PCOLCE2", "PRG4", "SFRP4",
  "SCARA5", "TRAC",
  "CCL2",  "SPSB1", "TNFAIP6",
  "CCDC146", "CCL19", "CD74", "TNFSF13B",
  "IL33", "SCN4B", "CTSH", "RBP5", "ACHE",
  "EFEMP1", "ITM2A", "MYOC", "GDF10",
  "COL11A1", "DPEP1", "WFDC1", "MEF2C",
  "COCH", "CRABP1", "FIBIN", "RSPO4", "NDNF", "SLITRK6",
  "POSTN", "LTBP2", "LRRC15",
  "IGFBP3", "SLC5A3", "WNT5A", "LUZP2", "DKK2", "RSPO4", "EDN3", "DIO2",
  "ANGPTL7", "APOD", "ECRG4",
  "ITGA6", "ITGB4", "TNNC1",
  "IGFBP2", "COL26A1",  "RAMP1", "OLFML2A",
  "CCER2", "KRT20"
))

p <- SCpubr::do_DotPlot(
  sample = extended_fb,
  features = genes,
  cluster = F,
  dot.scale = 9,
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
  filename = "./dotplot_fibro.pdf",
  plot = p2,
  device = "pdf",
  width = 20,
  height = 8
)


# Define a custom color scheme for each level 4 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_fb,
  group.by = "inherited_celltype_lvl_4_extended",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  legend.position = "right",
  colors.use = c(
    "A1" = "#66c2a5",
    "A2" = "#fc8d62",
    "A3" = "#8da0cb",
    "A4" = "#e78ac3",
    "Anagen DP (C5)" = "#bfffad",
    "B1" = "#ffd92f",
    "B2" = "#e5c494",
    "B3" = "#b3b3b3",
    "B4" = "#a1d99b",
    "C3" = "#9ecae1",
    "CAF1" = "#fdd0a2",
    "CAF2" = "#d4b9da",
    "D1" = "#3CB371",
    "D2" = "#fcbba1",
    "Dermal sheath (C1)" = "#c994c7",
    "Merkel cell_3" = "#cc00c6",
    "Outer bulge DP (C2)" = "#cab2d6",
    "RAMP1+ Fibro (E1)" = "#e07a5f"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./fibro.png",
  plot = p1,
  width = 13,       # in inches
  height = 9,
  dpi = 600
)

# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_fb) <- extended_fb$inherited_celltype_lvl_4

# The labels in 'inherited_celltype_lvl_4' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
celltypes_lvl4 <- extended_fb$inherited_celltype_lvl_4 %>% unique()
celltypes_lvl4
celltypes_clean <- celltypes_lvl4[!celltypes_lvl4 %in% c(
  "Extended", "Muscle progenitor", "STEAP4+ SMC", "Melanocyte_2",
  "DES+ SMC", "Supr. SG", "RERGL+ SMC", "Matrix", "RSG5+ EC"
)]

celltypes_clean

table(Idents(extended_fb), useNA = "always")

Idents(extended_fb) %>% table()
Idents(extended_fb) %>% unique()

Idents(extended_fb) <- extended_fb$inherited_celltype_lvl_4

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_fb,
  idents.keep = celltypes_clean,
  label = F,
  legend.position = "none",
  border.size = 3,
  pt.size = 0.2,
  na.value = "#363636"

)

p

# Highlight core immune cells on top of the extended dataset using distinct colors
# Create Supplementary Figure X
p2 <- p +
  scale_color_manual(values = c(
    "A1" = "#66c2a5",
    "A2" = "#fc8d62",
    "A3" = "#8da0cb",
    "A4" = "#e78ac3",
    "Anagen DP (C5)" = "#bfffad",
    "B1" = "#ffd92f",
    "B2" = "#e5c494",
    "B3" = "#b3b3b3",
    "B4" = "#a1d99b",
    "C3" = "#9ecae1",
    "CAF1" = "#fdd0a2",
    "CAF2" = "#d4b9da",
    "D1" = "#3CB371",
    "D2" = "#fcbba1",
    "Dermal sheath (C1)" = "#c994c7",
    "Merkel cell_3" = "#cc00c6",
    "Outer bulge DP (C2)" = "#cab2d6",
    "RAMP1+ Fibro (E1)" = "#e07a5f"
  ))

p2

# Save figure X
ggsave(
  filename = "./fibro_core_highlighted.png",
  plot = p2,
  width = 12,       # in inches
  height = 9,
  dpi = 600
)
