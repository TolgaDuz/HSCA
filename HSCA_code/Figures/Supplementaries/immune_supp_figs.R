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
# Load and inspect HSCA extended immune cells
# ===============================

extended_immune <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/immune_extended.rds")

extended_immune

DimPlot(extended_immune, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_immune, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_immune, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_immune, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_immune, group.by = "Core")
DimPlot(extended_immune, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_immune) <- "RNA"
extended_immune <- NormalizeData(extended_immune)

Idents(extended_immune) <- extended_immune$inherited_celltype_lvl_5_extended

new_order_ident <-  c(
  "FGFBP2+ NK_4",
  "SPINK2+ NK",
  "GZMK+ NK",
  "cMono_4",
  "ncMono_4",
  "EREG+ Mono",
  "Neutrophil_3",
  "GD-T cell_4",
  "CD4+ T cell_4",
  "Prolif. T cell_4",
  "CD8+ T cell_4",
  "naive CD4+ T cell",
  "naive CD8+ T cell",
  "Reg. T cell_4",
  "LC_4",
  "Prolif. DC_4",
  "Mature DC_4",
  "cDC1",
  "cDC2",
  "pDC",
  "Plasma cell_4",
  "Naive B cell_4",
  "Anti-inflammatory Mph_4",
  "Inflammatory Mph_4",
  "C3+ Mph",
  "LPL+ Mph",
  "Mast cell_3"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_immune))
setdiff(Idents(extended_immune), new_order_ident)

levels(Idents(extended_immune))
levels(extended_immune) <- new_order_ident
levels(Idents(extended_immune))

# Marker genes for skin cell types
genes <- unique(c(
  "PRF1", "SPON2", "FGFBP2", "MYOM2", "XCL1", "SPINK2", "GZMK", "KLRC1",
  "FCN1", "S100A12", "RNASE2",
  "MTSS1",
  "RHOC", "CDKN1C", "EREG",
  "AQP9", "FCGR3B", "CMTM2",
  "SPOCK2",
  "FXYD2",  "TRGC2", "ZNF683",
  "CCR6", "NTRK2", "IL26", "TOP2A",
  "CD8A",  "CD8B", "ZNF683",
  "SELL", "LEF1", "CD8A", "NELL2",
  "IKZF2", "CTLA4", "FOXP3",
  "CD207", "FCGBP", "CD1A",
  "MKI67",
  "LAMP3",  "WNT5B", "MARCKSL1", "CERS6", "SLCO5A1", "NAV1",
  "DNASE1L3", "CADM1",  "CLEC9A",
  "CD1C", "CLEC10A",
  "LILRA4", "CLEC4C", "LRRC26", "LAMP5", "SCT", "SMIM5", "JCHAIN",
  "IGKC", "CD79A", "MS4A1", "IGHM", "PAX5",
  "SELENOP", "F13A1", "C1QA",
  "CXCL2", "CXCL3", "SELENOP", "CTSL", "C1QA",
  "TREM2", "SPP1", "C3", "LPL",
  "TPSB2", "TPSAB1"
))

p <- SCpubr::do_DotPlot(
  sample = extended_immune,
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
  filename = "./dotplot_immun.pdf",
  plot = p2,
  device = "pdf",
  width = 22,
  height = 8
)


extended_immune$inherited_celltype_lvl_5_extended %>% table()

# Define a custom color scheme for each level 5 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_immune,
  group.by = "inherited_celltype_lvl_5_extended",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  legend.position = "right",
  colors.use = c(
    "FGFBP2+ NK_4"            = "#A6CEE3",
    "SPINK2+ NK"              = "#1F78B4",
    "GZMK+ NK"                = "#B2DF8A",
    "cMono_4"                 = "#FB9A99",
    "ncMono_4"                = "#E31A1C",
    "EREG+ Mono"              = "#8B4513",
    "Neutrophil_3"            = "#CAB2D6",
    "GD-T cell_4"             = "#6A3D9A",
    "CD4+ T cell_4"           = "#FFFF99",
    "Prolif. T cell_4"        = "#B15928",
    "CD8+ T cell_4"           = "#FF7F00",
    "naive CD4+ T cell"       = "#CCEBC5",
    "naive CD8+ T cell"       = "#B22222",
    "Reg. T cell_4"           = "#FDB462",
    "LC_4"                    = "#80B1D3",
    "Prolif. DC_4"            = "#d6ccb8",
    "Mature DC_4"             = "#8DD3C7",
    "cDC1"                    = "#006D2C",
    "cDC2"                    = "#999999",
    "pDC"                     = "#D9D9D9",
    "Plasma cell_4"           = "white",
    "Naive B cell_4"          = "#ECBF83",
    "Anti-inflammatory Mph_4" = "#FFED6F",
    "Inflammatory Mph_4"      = "#D95F02",
    "C3+ Mph"                 = "#A6761D",
    "LPL+ Mph"                = "#F781BF",
    "Mast cell_3"             = "#08306B"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./immune_umap_extended.png",
  plot = p1,
  width = 15,       # in inches
  height = 9,
  dpi = 600
)

# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_immune) <- extended_immune$inherited_celltype_lvl_5

# Assign cells not originating from the HSCA core to "Extended"
extended_immune$inherited_celltype_lvl_5[extended_immune$inherited_celltype_lvl_5 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_5' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_5 <- extended_immune$inherited_celltype_lvl_5 %>% unique()
inherited_celltype_lvl_5

extended_immune$inherited_celltype_lvl_5 %>% table()

extended_celltype_5 <- extended_immune$inherited_celltype_lvl_5_extended %>% unique()

setdiff(inherited_celltype_lvl_5, extended_celltype_5)
setdiff(extended_celltype_5, inherited_celltype_lvl_5)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_5[!inherited_celltype_lvl_5 %in% c(
  "Basal JZ_4",
  "Prolif. KC_3",
  "Venous 1 EC_4",
  "SEB-T",
  "Infund. KC_4",
  "Spinous KC_3",
  "CYP26B1+ SMC",
  "Matrix (Early)",
  "Venous 2 EC_4",
  "TM4SF1+ SMC",
  "Prolif. SEB-B",
  "A2_4",
  "Clear cell 2_4",
  "Basal SD_4",
  "Supr. Infund. KC_4",
  "Huxley's layer",
  "Dark cell_4",
  "Outer bulge",
  "Extended"
)]

setdiff(celltypes_clean, extended_celltype_5)
setdiff(extended_celltype_5, celltypes_clean)

celltypes_clean

table(Idents(extended_immune), useNA = "always")

Idents(extended_immune) %>% table()
Idents(extended_immune) %>% unique()

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_immune,
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
  "FGFBP2+ NK_4"            = "#A6CEE3",
  "SPINK2+ NK"              = "#1F78B4",
  "GZMK+ NK"                = "#B2DF8A",
  "Monocyte_3"              = "#FB9A99",
  "Neutrophil_3"            = "#CAB2D6",
  "GD-T cell_4"             = "#6A3D9A",
  "CD4+ T cell_4"           = "#FFFF99",
  "CD8+ T cell_4"           = "#FF7F00",
  "Reg. T cell_4"           = "#FDB462",
  "LC_4"                    = "#80B1D3",
  "Prolif. DC_4"            = "#d6ccb8",
  "Mature DC_4"             = "#8DD3C7",
  "cDC1"                    = "#006D2C",
  "cDC2"                    = "#999999",
  "Plasma cell_4"           = "white",
  "Naive B cell & pDC_4"    = "#ECBF83",
  "Anti-inflammatory Mph_4" = "#FFED6F",
  "Inflammatory Mph_4"      = "#D95F02",
  "C3+ Mph"                 = "#A6761D",
  "LPL+ Mph"                = "#F781BF",
  "Mast cell_3"             = "#08306B"
))


p2

# Save figure X
ggsave(
  filename = "./immune_core_highlighted.png",
  plot = p2,
  width = 15,
  height = 9,
  dpi = 600
)
