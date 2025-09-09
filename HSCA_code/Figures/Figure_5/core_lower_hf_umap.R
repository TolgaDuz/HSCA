library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCpubr)

setwd("~/HSCA_code/Figures/Figure_5/")

lower_hf_core <- readRDS("/path/to/figshare/visiumHD_lower_hairfollicle.rds")

# Inspect
lower_hf_core$inherited_celltype_lvl_5 %>% table()
DimPlot(lower_hf_core, group.by = "inherited_celltype_lvl_5", label = TRUE)

# Remove low quality cluster
lower_hf_core_sub <- subset(
  lower_hf_core,
  subset = inherited_celltype_lvl_5 %in% c("LQ"),
  invert = TRUE
)

# Create Figure 5a
p1 <- SCpubr::do_DimPlot(
  lower_hf_core_sub,
  group.by = "inherited_celltype_lvl_5",
  label = FALSE,
  border.size = 1.7,
  pt.size = 2,
  label.box = TRUE,
  legend.position = "none",
  repel = TRUE,
  colors.use = c(
    "Basal ORS" = c("#e9d8a6"),
    "CST6+ Inner bulge" = c("#8B80C4"),
    "Cuticle hair shaft" = c("#D6cf00"),
    "Henle's layer (upper)" = c("#a30000"),
    "Huxley's layer" = c("#ffb8b8"),
    "KRT83+ Cortex" = c("#02007a"),
    "LEF1+ Cortex" = c("#Fff70f"),
    "Matrix (Early)" = c("#084943"),
    "MUCL1+ Inner bulge" = c("#C6C6C6"),
    "Outer bulge" = c("#8B4513"),
    "Basal JZ_4" = c("#F18A85"),
    "Cuticle IRS" = c("#556b2f"),
    "Matrix (Late) & Henle's layer (lower)" = c("#fb8804"),
    "Supr. ORS" = c("#de9bfd"),
    "SHG" = c("#00f323"), #ändern!!!
    "Companion layer" = c("#3CA0EC"),
    "WFDC3+ Inner bulge" = c("#361509")
  )
)

p1

# Save Figure 5a
ggsave(
  filename = "./lower_hf_core_umap.pdf",
  plot = p1,
  width = 13,
  height = 9
)
