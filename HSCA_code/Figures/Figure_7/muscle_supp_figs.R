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
# Load and inspect HSCA extended muscle cells
# ===============================

extended_muscle <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/muscle_extended.rds")


extended_muscle

DimPlot(extended_muscle, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_muscle, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_muscle, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_muscle, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_muscle, group.by = "Core")
DimPlot(extended_muscle, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_muscle) <- "RNA"
extended_muscle <- NormalizeData(extended_muscle)

Idents(extended_muscle) <- extended_muscle$inherited_celltype_lvl_5_extended

new_order_ident <-  c(
  "RERGL+ SMC_4",
  "FRMD3+ SMC",
  "TM4SF1+ SMC",
  "CYP26B1+ SMC",
  "DES+ SMC_4",
  "Muscle progenitor_4",
  "Skeletal Muscle_3"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_muscle))
setdiff(Idents(extended_muscle), new_order_ident)

levels(Idents(extended_muscle))
levels(extended_muscle) <- new_order_ident
levels(Idents(extended_muscle))

genes <- unique(c(
  "ACTA2", "TAGLN",
  "RERGL", "PLN", "SBCG", "MYH11", "BCAM",  "NRGN",
  "STEAP4", "COL6A3",
  "FRMD3", "TSHZ2", "CD36",
  "TM4SF1", "GGT5", "CFD",
  "CYP26B1", "CFH", "P2RY14", "BASP1", "ADAMTS5", "APCDD1", "COL23A1", "GPM6B",
  "DES", "PCP4", "ACTG2", "PDE4D", "SLC8A1", "PRUNE2", "PCDH7", "ROBO2",
  "MYF5", "PAX7", "RBP1", "CHRNA1", "CADM2",  # muscle progenitor
  "COX6A2", "MYL1", "CSRP3", "MB", "TCAP", "MYLPF" # skeletal muscle
))

p <- SCpubr::do_DotPlot(
  sample = extended_muscle,
  features = genes,
  cluster = F,
  dot.scale = 10,
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
  filename = "./dotplot_muscle.pdf",
  plot = p2,
  device = "pdf",
  width = 17,
  height = 4
)



extended_muscle$inherited_celltype_lvl_5_extended %>% table()

# Define a custom color scheme for each level 5 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_muscle,
  group.by = "inherited_celltype_lvl_5_extended",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  legend.position = "right",
  colors.use = c(
    "RERGL+ SMC_4"        = "#CAB2D6",
    "FRMD3+ SMC"          = "#1F78B4",
    "TM4SF1+ SMC"         = "#B2DF8A",
    "CYP26B1+ SMC"        = "#FDB462",
    "DES+ SMC_4"          = "#B22222",
    "Muscle progenitor_4" = "#8B4513",
    "Skeletal Muscle_3"   = "#FFF999"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./muscle_umap_extended.png",
  plot = p1,
  width = 14,
  height = 9,
  dpi = 600
)

# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_muscle) <- extended_muscle$inherited_celltype_lvl_5

# Assign cells not originating from the HSCA core to "Extended"
extended_muscle$inherited_celltype_lvl_5[extended_muscle$inherited_celltype_lvl_5 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_5' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_5 <- extended_muscle$inherited_celltype_lvl_5 %>% unique()
inherited_celltype_lvl_5

extended_muscle$inherited_celltype_lvl_5 %>% table()

extended_celltype_5 <- extended_muscle$inherited_celltype_lvl_5_extended %>% unique()

setdiff(inherited_celltype_lvl_5, extended_celltype_5)
setdiff(extended_celltype_5, inherited_celltype_lvl_5)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_5[!inherited_celltype_lvl_5 %in% c(
  "Myoepithelial_4", "A2_4", "RSG5+ EC_4", "CNTNAP2+ Neuron_4",
  "Clear cell 2_4", "D2_4", "Basal luminal_4", "Venous 2 EC_4", "B3_4",
  "B2_4", "Dermal sheath (C1)_4", "Anagen DP (C5)_4", "SEB-T", "Capillary EC_3",
  "B1_4", "Extended"
)]

celltypes_clean

table(Idents(extended_muscle), useNA = "always")

Idents(extended_muscle) %>% table()
Idents(extended_muscle) %>% unique()

Idents(extended_muscle) <- extended_muscle$inherited_celltype_lvl_5

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_muscle,
  idents.keep = celltypes_clean,
  label = F,
  legend.position = "right",
  border.size = 3,
  pt.size = 0.2,
  na.value ="#363636"
)

p

# Highlight core immune cells on top of the extended dataset using distinct colors
# Create Supplementary Figure X
p2 <- p + scale_color_manual(values = c(
    "RERGL+ SMC_4"            = "#CAB2D6",
    "FRMD3+ SMC"              = "#1F78B4",
    "TM4SF1+ SMC"             = "#B2DF8A",
    "CYP26B1+ SMC"            = "#FDB462",
    "DES+ SMC_4"              = "#B22222",
    "Muscle progenitor_4"     = "#8B4513",
    "Skeletal muscle_3"       = "#FFF999"
  ))


p2

# Save figure X
ggsave(
  filename = "./muscle_core_highlighted.png",
  plot = p2,
  width = 15,
  height = 9,
  dpi = 600
)
