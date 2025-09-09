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
# Load and inspect HSCA extended IFE
# ===============================

extended_ife <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/ife_extended.rds")

extended_ife

DimPlot(extended_ife, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_ife, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_ife, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_ife, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_ife, group.by = "Core")
DimPlot(extended_ife, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_ife) <- "RNA"
extended_ife <- NormalizeData(extended_ife)

Idents(extended_ife) <- extended_ife$inherited_celltype_lvl_3_extended

Idents(extended_ife) %>% unique()

new_order_ident <-  c(
  "Basal KC",
  "Prolif. KC",
  "Spinous KC",
  "Granular KC",
  "Cornified KC",
  "Adipocyte_2"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_ife))
setdiff(Idents(extended_ife), new_order_ident)

levels(Idents(extended_ife))
levels(extended_ife) <- new_order_ident
levels(Idents(extended_ife))

genes <- unique(c(
  "KRT15", "KRT5", "COL17A1", "POSTN",
  "MKI67", "TOP2A", "CENPF", "PTTG1",
  "APOE", "KRT1", "KRT10", "DMKN",
  "KRT2", "IVL",
  "LOR", "FLG", "FLG2",
  "LCE1A", "LCE1B", "LCE1C", "SPRR1B", "SPRR1A",
  "ADIPOQ", "PLIN1", "FABP4"
))

p <- SCpubr::do_DotPlot(
  sample = extended_ife,
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
  filename = "./dotplot_ife.pdf",
  plot = p2,
  device = "pdf",
  width = 13,
  height = 8
)

extended_ife$inherited_celltype_lvl_3_extended %>% table()

# Define a custom color scheme for each level 5 cell type
# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_ife,
  group.by = "inherited_celltype_lvl_3_extended",
  label = FALSE,
  border.size = 1.7,
  pt.size = 2,
  legend.position = "right",
  colors.use = c(
    "Basal KC"     = "#005f73",
    "Prolif. KC"   = "#3CA0EC",
    "Spinous KC"   = "#556b2f",
    "Granular KC"  = "#ae2012",
    "Cornified KC" = "#C6C6C6",
    "Adipocyte_2"  = "#FF0000"
  )
)

p1

# Save Fig. X
ggsave(
  filename = "./ife_umap_extended.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)


# Assign identities based on the HSCA core cell type annotations (inherited_celltype_lvl_X)
Idents(extended_ife) <- extended_ife$inherited_celltype_lvl_3

# Assign cells not originating from the HSCA core to "Extended"
extended_ife$inherited_celltype_lvl_3[extended_ife$inherited_celltype_lvl_3 == ""] <- "Extended"

# The labels in 'inherited_celltype_lvl_3' were derived during the HSCA core atlas construction.
# After de novo integration and annotation of the HSCA extended, some HSCA core labels may not align perfectly
# with the final HSCA extended cell type clusters.These cases are rare, and such mismatched core cells
# will not be highlighted to ensure consistency.
inherited_celltype_lvl_3 <- extended_ife$inherited_celltype_lvl_3 %>% unique()
inherited_celltype_lvl_3

extended_ife$inherited_celltype_lvl_3 %>% table()

extended_celltype_3 <- extended_ife$inherited_celltype_lvl_3_extended %>% unique()

setdiff(inherited_celltype_lvl_3, extended_celltype_3)
setdiff(extended_celltype_3, inherited_celltype_lvl_3)

# Remove mismatched or extended-only labels to create a clean set of core cell types,
# that will be highlighted on top of the HSCA extended
celltypes_clean <- inherited_celltype_lvl_3[inherited_celltype_lvl_3 %in%  c(
  "Basal KC",
  "Prolif. KC",
  "Spinous KC",
  "Granular KC",
  "Cornified KC",
  "Adipocyte_2",
  "Infundibulum",
  "Isthmus"
)]

celltypes_clean

table(Idents(extended_ife), useNA = "always")

Idents(extended_ife) %>% table()
Idents(extended_ife) %>% unique()

Idents(extended_ife) <- extended_ife$inherited_celltype_lvl_3

# Plot cells with NA values (extended-only cells) highlighted in black
p <- SCpubr::do_DimPlot(
  sample = extended_ife,
  idents.keep = celltypes_clean,
  label = F,
  legend.position = "right",
  border.size = 1.7,
  pt.size = 2,
  na.value = "#363636"

)

p

# Highlight core immune cells on top of the extended dataset using distinct colors
# Create Supplementary Figure X
p2 <- p +
  scale_color_manual(values = c(
    "Basal KC"        = "#005f73",
    "Prolif. KC"      = "#3CA0EC",
    "Spinous KC"      = "#556b2f",
    "Granular KC"     = "#ae2012",
    "Cornified KC"    = "#C6C6C6",
    "Infundibulum"    = "#363636",
    "Isthmus"         = "#363636",
    "Adipocyte_2"     = "#FF0000"
  ))

p2

# Save figure X
ggsave(
  filename = "./ife_core_highlighted.png",
  plot = p2,
  width = 15,
  height = 9,
  dpi = 600
)

# ===============================
# Highlight anatomical region
# ===============================

extended_ife$anatomical_region_level2 %>% table()
extended_ife$anatomical_region_level2 %>% unique()

Idents(extended_ife) <- extended_ife$anatomical_region_level2

# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_ife,
  group.by = "anatomical_region_level2",
  label = FALSE,
  border.size = 1.7,
  pt.size = 2,
  legend.position = "right",
  colors.use = c(
    "Face"            = "#005f73",
    "Forearm"         = "cyan",
    "Flank"           = "#3CB043",
    "Lower Back"      = "#8B80C4",
    "Arm"             = "yellow",
    "diverse"         = "#B22222",
    "nan"             = "#363636",
    "Scalp"           = "#00f323",
    "Breast"          = "#999999",
    "Abdomen"         = "#CCEBC5",
    "Pelvic Region"   = "#8B4513",
    "Trunk"           = "#e9d8a6",
    "Elbow"           = "#de9bfd",
    "Neck"            = "#3CA0EC",
    "Genitalia"       = "red",
    "Hip"             = "#fb8804",
    "Foot"            = "#02007a",
    "Hand"            = "#556b2f"
  )
)

p1

# Save Figure X
ggsave(
  filename = "./ife_umap_anatomical_regions.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)

# ===============================
# Highlight Accesion source
# ===============================

extended_ife$Accession_source %>% unique()

# Create Figure X
p1 <- SCpubr::do_DimPlot(
  extended_ife,
  group.by = "Accession_source",
  label = FALSE,
  border.size = 1.7,
  pt.size = 2,
  legend.position = "right",
  colors.use = c(
    "HRA000395"       = "#A6CEE3",
    "GSE138669"       = "#1F78B4",
    "GSE132802"       = "#B2DF8A",
    "GSE275491"       = "#FB9A99",
    "GSE186476"       = "#E31A1C",
    "GSE147424"       = "#8B4513",
    "GSE129611"       = "#CAB2D6",
    "GSE159929"       = "#6A3D9A",
    "EGAS00001002927" = "#FFFF99",
    "GSE191067"       = "#B15928",
    "GSE274955"       = "#FF7F00",
    "EMTAB13084"      = "#CCEBC5",
    "PRJNA754272"     = "#B22222",
    "GSE130973"       = "#FDB462",
    "GSE179633"       = "#80B1D3",
    "GSE176415"       = "#d6ccb8",
    "GSE182208"       = "#8DD3C7",
    "GSE173205"       = "#006D2C",
    "GSE153760"       = "#FFED6F",
    "PRJNA797897"     = "#D9D9D9",
    "GSE151177"       = "white",
    "GSE181316"       = "#ECBF83",
    "GSE147482"       = "#999999",
    "GSE202352"       = "#D95F02",
    "GSE162183"       = "#A6761D",
    "GSE165816"       = "#F781BF",
    "GSE241132"       = "#08306B",
    "GSE265972"       = "black"
  )
)

p1

ggsave(
  filename = "./ife_umap_Accession_source.png",
  plot = p1,
  width = 15,
  height = 9,
  dpi = 600
)
