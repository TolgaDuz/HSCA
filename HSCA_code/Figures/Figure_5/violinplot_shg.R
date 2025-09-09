library(Seurat)
library(scCustomize)
library(tidyverse)

setwd("~/HSCA_code/Figures/Figure_5/")

# ==== Check SHG for HSCA core dataset ====
lower_hf_core <- readRDS("/path/to/figshare/visiumHD_lower_hairfollicle.rds")

lower_hf_core$inherited_celltype_lvl_5 %>% table()

# Remove low quality cells
lower_hf_core <- subset(lower_hf_core, subset = inherited_celltype_lvl_5 != "LQ")
DimPlot(lower_hf_core, group.by = "inherited_celltype_lvl_5", label = T)

Idents(lower_hf_core) <- lower_hf_core$inherited_celltype_lvl_5

DefaultAssay(lower_hf_core) <- "RNA"
lower_hf_core <- NormalizeData(lower_hf_core)

# Set order for violinplot visualization
new_order_ident <- c(
  "Basal JZ_4",
  "Outer bulge",
  "MUCL1+ Inner bulge",
  "CST6+ Inner bulge",
  "WFDC3+ Inner bulge",
  "Basal ORS",
  "Supr. ORS",
  "Companion layer",
  "Henle's layer (upper)",
  "Matrix (Late) & Henle's layer (lower)",
  "Huxley's layer",
  "Cuticle IRS",
  "Cuticle hair shaft",
  "LEF1+ Cortex",
  "KRT83+ Cortex",
  "Matrix (Early)",
  "SHG"
)

levels(lower_hf_core) <- new_order_ident
DimPlot(lower_hf_core)

human_colors_list <- c(
  "#F18A85", # Basal JZ
  "#8B4513", # Outer bulge
  "#C6C6C6", # MUCL1 inner bulge
  "#8B80C4", # CST6+
  "#361509", # WFDC3+
  "#e9d8a6", # Basal ORS
  "#de9bfd", # Supr. ORS
  "#3CA0EC", # Companion layer
  "#a30000", # Henle's layer (upper)
  "#fb8804", # Henles layer (lower)
  "#ffb8b8", # Huxley
  "#556b2f", # Cuticle IRS
  "#D6cf00", # Cuticle hair shaft
  "#Fff70f", # LEF1+ cortex
  "#02007a", # KRT83+ cortex
  "#084943", # Matrix early
  "#00f323"  # SHG
)


# ==== DGEA ====
DefaultAssay(lower_hf_core)
mark <- FindMarkers(lower_hf_core, ident.1 = "SHG", only.pos = T)

mark$gene <- rownames(mark)
head(mark, 15)

FeaturePlot(lower_hf_core, rownames(mark)[1:12])
FeaturePlot(lower_hf_core, rownames(mark)[5:8])

# Marker genes
gene_list_plot <- c(
  "CUX2",
  "KREMEN2",
  "TOX",
  "TLL1",
  "FAT3",
  "RARB",
  "NAV2",
  "PTCH1",
  "LAMB1",
  "TMTC1",
  "FRAS1",
  "PTCH2"
)

# Create Figure 5c
p1 <- Stacked_VlnPlot(
  seurat_object = lower_hf_core,
  features = gene_list_plot,
  x_lab_rotate = TRUE,
  colors_use = human_colors_list
)

p1

# Save Figure 5c
ggsave(
  filename = "./violin_shg_core.pdf",
  plot = p1,
  device = "pdf",
  width = 7,
  height = 8
)


# ==== Check SHG for Visium HD dataset ====
visiumHD_lower_hf <- readRDS("../Figure_4/visiumHD_lower_hairfollicle.rds")

DimPlot(visiumHD_lower_hf, group.by = "celltype", label = T)

Idents(visiumHD_lower_hf) <- visiumHD_lower_hf$celltype

table(Idents(visiumHD_lower_hf))
DefaultAssay(visiumHD_lower_hf) <- "Spatial.008um"
visiumHD_lower_hf <- NormalizeData(visiumHD_lower_hf)

Idents(visiumHD_lower_hf) <- visiumHD_lower_hf$celltype

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
DimPlot(visiumHD_lower_hf)

human_colors_list <- c(
  "#8B4513", # Outer bulge
  "#C6C6C6", # MUCL1+ inner bulge
  "#8B80C4", # CST6+
  "#e9d8a6", # Basal ORS
  "#de9bfd", # Supr. ORS
  "#3CA0EC", # Companion layer
  "#a30000", # Henle's layer (upper)
  "#ff2a1f", # Henle's layer (middle)
  "#fb8804", # Henle's layer (lower)
  "#ffb8b8", # Huxleys layer
  "#D6cf00", # Cuticle
  "#Fff70f", # LEF1+ cortex
  "#02007a", # KRT83+ cortex
  "#E0115F", # cortex medulla
  "#FFFFFF", # Melanocyte
  "#196fb0", # Matrix late
  "#084943", # Matrix early
  "#bfffad", # Dermal papilla
  "#94d2db", # Dermal sheath
  "#00fbff"  # Lower catagen bulge
)


gene_list_plot <- c(
  "CUX2",
  "KREMEN2",
  "TOX",
  "TLL1",
  "FAT3",
  "RARB",
  "NAV2",
  "PTCH1",
  "LAMB1",
  "TMTC1",
  "FRAS1",
  "PTCH2"
)

# Create Figure 5d
p1 <- Stacked_VlnPlot(
  seurat_object = visiumHD_lower_hf,
  features = gene_list_plot,
  x_lab_rotate = TRUE,
  colors_use = human_colors_list
)

p1

# Save Figure 5d
ggsave(
  filename = "./violin_shg_visiumhd.pdf",
  plot = p1,
  device = "pdf",
  width = 7,
  height = 8
)
