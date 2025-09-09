set.seed(42)

library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SCpubr)
library(SeuratWrappers)
library(viridis)

setwd("~/HSCA_code/Figures/Figure_2/")

# Input the PSU core object
psu <- readRDS("./Core_psu.rds")

# Inspect object
DimPlot(
  psu,
  reduction = "umap",
  group.by = "inherited_celltype_lvl_4",
  label = TRUE,
  raster = FALSE,
  repel = TRUE
)

DimPlot(
  psu,
  reduction = "umap",
  group.by = "Dataset",
  label = TRUE,
  raster = FALSE
)

table(Idents(psu))
DefaultAssay(psu) <- "RNA"
psu <- NormalizeData(psu)

Idents(psu) <- psu$inherited_celltype_lvl_4
Idents(psu) %>% unique()

new_order_ident <- c(
  "Supr. Infund. KC",
  "Infund. KC",
  "Supr. JZ & SD",
  "Basal SD",
  "Basal JZ",
  "SG progenitors",
  "Basal KC_3",
  "Basal SG",
  "Supr. SG",
  "ORS",
  "Outer bulge",
  "Inner bulge",
  "IRS",
  "Hair shaft",
  "Matrix"
)

# Check for marker genes
Idents(psu) <- psu$celltype_lvl_4
markers <- FindMarkers(psu, ident.1 = "IRS", only.pos = TRUE)

head(markers, 20)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(psu))
setdiff(Idents(psu), new_order_ident)

levels(Idents(psu))
levels(psu) <- new_order_ident
levels(Idents(psu))

# Marker genes for PSU cell types
genes <- unique(
  c(
    "FLG", "SERPINB4", "PI3", "KRT1", "IVL",
    "S100A8", "S100A7", "S100A9", "KRT6A",
    "CALB2",  "KRT79", "GATA6",
    "KLK6", "LCN2", "RHCG", "CDA",
    "PTN", "C1QTNF12", "KRT15",  "POSTN", "KRT79", "LGR6",
    "NNAT", "IL1R2", "TINAGL1",
    "WFDC2", "KRT7",
    "PNPLA3", "PLIN5",
    "ANGPTL7", "COMP", "FGF18", "CRLF1", "SERPINA3",
    "DIO2", "WIF1",
    "KRT75", "NES", "GABRP",
    "TCHH", "KRT25", "KRT71", "KRT72", "KRT73",
    "KRT82", "ACTBL2",
    "SELENBP1", "KRT85", "GPNMB", "LEF1", "KRT35",
    "KRT83",
    "LPAR6", "MYCN",
    "TOP2A", "MKI67"
  )
)

p <- SCpubr::do_DotPlot(
  sample = psu,
  features = genes,
  cluster = FALSE,
  dot.scale = 10,
  legend.position = "right",
  font.size = 18,
  legend.length = 10,
  min.cutoff = 0,
  zscore.data = TRUE,
  flip = FALSE
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

p2 <- p + ggplot2::scale_fill_gradientn(
  colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  limits = c(0, upper_limit),
  oob = scales::squish
)

p2

# Create Supp. fig. 4
ggsave(
  filename = "./psu_dotplot.pdf",
  plot = p2,
  device = "pdf",
  width = 18,
  height = 8
)
