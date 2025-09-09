library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCpubr)

setwd("~/HSCA_code/Figures/Figure_2/")

# Input the HSCA core object (available on Figshare)
core <- readRDS("/path/to/figshare/HSCA_core.rds")

# Inspect
core$Dataset %>% table() %>% sort(decreasing = TRUE)
core$inherited_celltype_lvl_2 %>% unique()

p1 <- SCpubr::do_DimPlot(core,
  group.by = "inherited_celltype_lvl_2",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  colors.use = c(
    "Pilosebaceous Unit" = c("#ECBF83"),
    "Interfollicular Epidermis" = c("#556b2f"),
    "Vascular EC" = c("#003166"),
    "Melanocyte" = c("#FFFFFF"),
    "Fibroblast" = c("#b30030"),
    "Eccrine sweat Gland" = c("#732191"),
    "Lymphoid" = c("#0062cc"),
    "Myeloid" = c("#00B4D8"),
    "Lymphatic EC" = c("gray"),
    "Neuron" = c("#F18A85"),
    "Schwann cell" = c("#ffff7f"),
    "Adipocyte" = c("#FF0000"),
    "Muscle" = c("#007313")
  )
)

p1

ggsave(
  filename = "./HSCA_core_lvl2.pdf",
  plot = p1,
  device = "pdf",
  width = 8,
  height = 6,
)

p2 <- SCpubr::do_BarPlot(sample = core,
  group.by = "celltype_lvl_2",
  split.by = "Accession_source",
  plot.title = "Number of cells per cluster",
  flip = TRUE,
  order = TRUE,
  colors.use = c(
    "Pilosebaceous Unit" = c("#ECBF83"),
    "Interfollicular Epidermis" = c("#556b2f"),
    "Vascular EC" = c("#003166"),
    "Melanocyte" = c("#FFFFFF"),
    "Fibroblast" = c("#b30030"),
    "Eccrine sweat Gland" = c("#732191"),
    "Lymphoid" = c("#0062cc"),
    "Myeloid" = c("#00B4D8"),
    "Lymphatic EC" = c("gray"),
    "Neuron" = c("#F18A85"),
    "Schwann cell" = c("#ffff7f"),
    "Adipocyte" = c("#FF0000"),
    "Muscle" = c("#007313")
  )
)

p2

ggsave(
  filename = "./HSCA_core_lvl2_barplot.pdf",
  plot = p2,
  device = "pdf",
  width = 8,
  height = 6,
)


# Highlight PSU
p1 <- SCpubr::do_DimPlot(core,
  group.by = "inherited_celltype_lvl_2",
  label = FALSE,
  border.size = 3,
  pt.size = 0.2,
  colors.use = c(
    "Pilosebaceous Unit" = c("#ECBF83"),
    "Interfollicular Epidermis" = c("#D3D3D3"),
    "Vascular EC" = c("#D3D3D3"),
    "Melanocyte" = c("#D3D3D3"),
    "Fibroblast" = c("#D3D3D3"),
    "Eccrine sweat Gland" = c("#D3D3D3"),
    "Lymphoid" = c("#D3D3D3"),
    "Myeloid" = c("#D3D3D3"),
    "Lymphatic EC" = c("gray"),
    "Neuron" = c("#D3D3D3"),
    "Schwann cell" = c("#D3D3D3"),
    "Adipocyte" = c("#D3D3D3"),
    "Muscle" = c("#D3D3D3")
  )
)

p1

ggsave(
  filename = "HSCA_core_highlight_psu_gray.png",
  plot = p1,
  width = 11,
  height = 9,
  dpi = 600
)
