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
library(scCustomize)

setwd("~/HSCA_code/Figures/Figure_7")

# ===============================
# Load and inspect HSCA extended psu
# ===============================

extended_psu <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/psu_extended.rds")

extended_psu

DimPlot(extended_psu, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_psu, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_psu, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_psu, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_psu, group.by = "Core")
DimPlot(extended_psu, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_psu) <- "RNA"
extended_psu <- NormalizeData(extended_psu)

Idents(extended_psu) <- extended_psu$inherited_celltype_lvl_5_extended

new_order_ident <-  c(
  "Supr. Infund. KC_4",
  "Infund. KC_4",
  "Supr. JZ & SD_4",
  "Basal SD_4",
  "Basal JZ_4",
  "SG progenitors_4",
  "Basal KC_3",
  "SEB-B",
  "Prolif. SEB-B",
  "SEB-T",
  "SEB-1",
  "SEB-2",
  "SEB-2L",
  "Outer bulge",
  "MUCL1+ Inner bulge",
  "CST6+ Inner bulge",
  "WFDC3+ Inner bulge",
  "Basal ORS",
  "Supr. ORS",
  "Companion layer",
  "Henle's layer (upper)",
  "Matrix (Late) & Henle's layer (lower)",
  "Matrix (Early)",
  "Cuticle IRS",
  "Cuticle hair shaft",
  "LEF1+ Cortex",
  "KRT83+ Cortex"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_psu))
setdiff(Idents(extended_psu), new_order_ident)

levels(Idents(extended_psu))
levels(extended_psu) <- new_order_ident
levels(Idents(extended_psu))

genes <- unique(c(
  "FLG", "SERPINB4", "PI3", "KRT1", "IVL",
  "S100A8", "S100A7", "S100A9", "KRT6A",
  "CALB2",  "KRT79", "GATA6",
  "KLK6", "LCN2", "RHCG", "CDA",
  "PTN", "C1QTNF12", "KRT15",  "POSTN", "KRT79", "LGR6",
  "NNAT", "IL1R2", "MKI67", "TINAGL1",
  "WFDC2", "KRT7", "TBC1D4",
  "HAO2", "ACO1",
  "PNPLA3", "PLIN5",
  "CRAT", "SEC14L6", "DNASE1L2",
  "DIO2", "WIF1",
  "MUCL1", "CTSV", "CST6", "WFDC3",
  "ANGPTL7", "COMP", "FGF18", "CRLF1",
  "SERPINA3",
  "KRT75", "NES", "GABRP",
  "VEGFA", "CA9", "MEST",
  "CHI3L1", "COMP",
  "CAPN8",
  "LPAR6", "MYCN",
  "TOP2A",
  "KRT71", "KRT72", "KRT73",
  "KRT82", "PROCR", "ACTBL2",
  "SELENBP1", "KRT85", "GPNMB", "LEF1", "KRT35",
  "VSIG8", "KRT81", "KRT83", "KRT86"
))

p <- SCpubr::do_DotPlot(
  sample = extended_psu,
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
  filename = "./dotplot_psu.pdf",
  plot = p2,
  device = "pdf",
  width = 22,
  height = 8
)

extended_psu$inherited_celltype_lvl_5_extended %>% table()

# Define a custom color scheme for each level 5 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_psu,
  group.by = "inherited_celltype_lvl_5_extended",
  label = FALSE,
  border.size = 1.3,
  pt.size = 2.2,
  legend.position = "right",
  colors.use = c(
    "Supr. Infund. KC_4"                    = "#ae2012",
    "Infund. KC_4"                          = "#800080",
    "Supr. JZ & SD_4"                       = "#FF0000",
    "Basal SD_4"                            = "#FC6A03",
    "Basal JZ_4"                            = "#F18A85",
    "SG progenitors_4"                      = "#d6ccb8",
    "Basal KC_3"                            = "#005f73",
    "SEB-B"                                 = "#00B4D8",
    "Prolif. SEB-B"                         = "#3CB043",
    "SEB-T"                                 = "#8B80C4",
    "SEB-1"                                 = "#d7961d",
    "SEB-2"                                 = "#B22222",
    "SEB-2L"                                = "#FFE130",
    "Outer bulge"                           = "#8B4513",
    "MUCL1+ Inner bulge"                    = "#C6C6C6",
    "CST6+ Inner bulge"                     = "#8B80C4",
    "WFDC3+ Inner bulge"                    = "#361509",
    "Basal ORS"                             = "#e9d8a6",
    "Supr. ORS"                             = "#de9bfd",
    "Companion layer"                       = "#3CA0EC",
    "Henle's layer (upper)"                 = "#a30000",
    "Matrix (Late) & Henle's layer (lower)" = "#fb8804",
    "Matrix (Early)"                        = "#084943",
    "Cuticle IRS"                           = "#556b2f",
    "Cuticle hair shaft"                    = "#D6cf00",
    "LEF1+ Cortex"                          = "#Fff70f",
    "KRT83+ Cortex"                         = "#02007a",
    "SHG"                                   = "#00f323"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./psu_umap_extended.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)

# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_psu) <- extended_psu$inherited_celltype_lvl_5

# Assign cells not originating from the HSCA core to "Extended"
extended_psu$inherited_celltype_lvl_5[extended_psu$inherited_celltype_lvl_5 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_5' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_5 <- extended_psu$inherited_celltype_lvl_5 %>% unique()
inherited_celltype_lvl_5

extended_psu$inherited_celltype_lvl_5 %>% table()

extended_celltype_5 <- extended_psu$inherited_celltype_lvl_5_extended %>% unique()

setdiff(inherited_celltype_lvl_5, extended_celltype_5)
setdiff(extended_celltype_5, inherited_celltype_lvl_5)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_5[!inherited_celltype_lvl_5 %in%  c(
  "Basal_JZ_4",
  "Spinous KC_3",
  "Basal luminal_4",
  "Prolif. KC_3",
  "Huxley's layer",
  "TM4SF1+ SMC",
  "S100P+ luminal_4",
  "Clear cell 2_4",
  "Cornified KC_3",
  "Granular KC_3",
  "A2_4",
  "LC_4",
  "Myoepithelial_4",
  "Clear cell 1_4",
  "Extended"
)]

celltypes_clean

table(Idents(extended_psu), useNA = "always")

Idents(extended_psu) %>% table()
Idents(extended_psu) %>% unique()

Idents(extended_psu) <- extended_psu$inherited_celltype_lvl_5

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_psu,
  idents.keep = celltypes_clean,
  label = F,
  legend.position = "right",
  border.size = 1.3,
  pt.size = 2.2,
  na.value = "#363636"
)

p

# Highlight core immune cells on top of the extended dataset using distinct colors
# Create Supplementary Figure X
p2 <- p + scale_color_manual(values = c(
  "Supr. Infund. KC_4"                    = "#ae2012",
  "Infund. KC_4"                          = "#800080",
  "Supr. JZ & SD_4"                       = "#FF0000",
  "Basal SD_4"                            = "#FC6A03",
  "Basal JZ_4"                            = "#F18A85",
  "SG progenitors_4"                      = "#d6ccb8",
  "Basal KC_3"                            = "#005f73",
  "SEB-B"                                 = "#00B4D8",
  "Prolif. SEB-B"                         = "#3CB043",
  "SEB-T"                                 = "#8B80C4",
  "SEB-1"                                 = "#d7961d",
  "SEB-2"                                 = "#B22222",
  "SEB-2L"                                = "#FFE130",
  "Outer bulge"                           = "#8B4513",
  "MUCL1+ Inner bulge"                    = "#C6C6C6",
  "CST6+ Inner bulge"                     = "#8B80C4",
  "WFDC3+ Inner bulge"                    = "#361509",
  "Basal ORS"                             = "#e9d8a6",
  "Supr. ORS"                             = "#de9bfd",
  "Companion layer"                       = "#3CA0EC",
  "Henle's layer (upper)"                 = "#a30000",
  "Matrix (Late) & Henle's layer (lower)" = "#fb8804",
  "Matrix (Early)"                        = "#06402B",
  "Cuticle IRS"                           = "#556b2f",
  "Cuticle hair shaft"                    = "#D6cf00",
  "LEF1+ Cortex"                         = "#Fff70f",
  "KRT83+ Cortex"                         = "#02007a",
  "SHG"                                   = "#00f323"
))

p2

# Save figure X
ggsave(
  filename = "./psu_core_highlighted.png",
  plot = p2,
  width = 20,
  height = 9,
  dpi = 600
)

# ==============================================================
# Label transfer uncertainty grouped by cell type (level 2) and aggregated from level 5 annotations
# ==============================================================

extended_atlas <- readRDS("./atlas/final_filt/hsca_extended.rds")

extended_atlas_only <- subset(extended_atlas, subset = Core == "No")

extended_atlas_only$inherited_celltype_lvl_1_extended %>% table()

colors_list <- c(
  "#b30030",  # FB
  "#FFFFFF", # Mela
  "gray", # LEC
  "#003166", # VEC
  "#0062cc", # Lymphoid
  "#363636", # Erythrocyte
  "#007313", # Muscle
  "#556b2f", # IFE
  "#00B4D8", # Myeloid
  "#F18A85", # Neuron
  "#ffff7f", # Schwann cell
  "#363636", # Schwann (ex vivo)
  "#363636", # Chondrocyte
  "#732191", # ESG
  "#ECBF83", # PSU
  "#363636", # Neural progenitor
  "#FF0000"  # adipocyte
)

Idents(extended_atlas_only) <- extended_atlas_only$inherited_celltype_lvl_2_extended

Idents(extended_atlas_only)

# Create Figure 7x
p1 <- Stacked_VlnPlot(
  seurat_object = extended_atlas_only,
  features = "inherited_celltype_lvl_5_transfer_uncert",
  x_lab_rotate = TRUE,
  colors_use = colors_list
)

p1_with_lines <- p1 & geom_hline(yintercept = 0.2, linetype = "dashed", color = "red")

p1_with_lines

ggsave(
  filename = "./violin_labeltransfer.pdf",
  plot = p1_with_lines,
  width = 15,
  height = 7
)
