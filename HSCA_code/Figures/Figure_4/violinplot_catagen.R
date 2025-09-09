library(Seurat)
library(scCustomize)
library(tidyverse)

setwd("~/HSCA_code/Figures/Figure_4/")

visiumHD_lower_hf <- readRDS("/path/to/figshare/visiumHD_lower_hairfollicle.rds")

# Inspect
DimPlot(visiumHD_lower_hf, group.by = "seurat_clusters", label = T)
DimPlot(visiumHD_lower_hf, group.by = "celltype", label = T)

Idents(visiumHD_lower_hf) <- visiumHD_lower_hf$celltype

table(Idents(visiumHD_lower_hf))
DefaultAssay(visiumHD_lower_hf) <- "Spatial.008um"
visiumHD_lower_hf <- NormalizeData(visiumHD_lower_hf)

# Set order for violinplot visualization
new_order_ident <- c(
  "Outer bulge",
  "MUCL1+ Inner bulge",
  "CST6+ Inner bulge",
  "Basal ORS",
  "Supr. ORS",
  "Companion layer",
  "Henle's layer (upper)",
  "Henle's layer (middle)",
  "Henle's layer (lower)",
  "Huxley's layer",
  "Cuticle",
  "LEF1+ Cortex",
  "KRT83+ Cortex",
  "Cortex/Medulla",
  "Melanocyte",
  "Matrix (Late)",
  "Matrix (Early)",
  "Dermal papilla",
  "Dermal sheath",
  "lower catagen bulge"
)

levels(visiumHD_lower_hf) <- new_order_ident
Idents(visiumHD_lower_hf)

# Define color palette for human hair follicle cell types
human_colors_list <- c(
  "#8B4513", # Outer bulge
  "#8B80C4", # MUCL1+ inner bulge
  "#C6C6C6", # CST6+
  "#e9d8a6", # Basal ORS
  "#556b2f", # Supr. ORS
  "#3CA0EC", # Companion layer
  "#a30000", # Henle's layer (upper)
  "#ff2a1f", # Henle's layer (middle)
  "#fb8804", # Henle's layer (lower)
  "#ffb8b8", # Huxley's layer
  "#D6cf00", # Cuticle
  "#Fff70f", # LEF1+ cortex
  "#02007a", # KRT83+ cortex
  "#E0115F", # Cortex/Medulla
  "#FFFFFF", # Melanocyte
  "#196fb0", # Matrix (Late)
  "#000000", # Matrix (Early)
  "#bfffad", # Dermal papilla
  "#94d2db", # Dermal sheath
  "#00fbff" # Lower catagen bulge
)


# ==== DGEA ====
DefaultAssay(visiumHD_lower_hf) <- "SCT"

# Run differential gene expression
visiumHD_lower_hf <- PrepSCTFindMarkers(visiumHD_lower_hf)
mark_spatial <- FindMarkers(
  visiumHD_lower_hf,
  ident.1 = c("lower catagen bulge"),
  only.pos = T,
  assay = "SCT"
)

head(mark_spatial, 15)
FeaturePlot(visiumHD_lower_hf, rownames(mark_spatial)[1:12])
FeaturePlot(visiumHD_lower_hf, rownames(mark_spatial)[13:24])

# Catagen markers
gene_list_plot <- c(
  "SULT1E1",
  "THBS1",
  "PAPLN",
  "CLU"
)

# Create Figure 4f
p1 <- Stacked_VlnPlot(
  seurat_object = visiumHD_lower_hf,
  features = gene_list_plot,
  x_lab_rotate = TRUE,
  colors_use = human_colors_list
)

p1

# Save Figure 4f
ggsave(
  filename = "./violin_lower_catagen.pdf",
  plot = p1,
  device = "pdf",
  width = 8,
  height = 6
)
