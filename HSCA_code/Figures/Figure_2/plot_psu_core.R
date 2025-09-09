library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCpubr)

setwd("~/HSCA_code/Figures/Figure_2/")

# Input the PSU core object (available on Figshare)
psu <- readRDS("/path/to/figshare/PSU_core.rds")

psu$inherited_celltype_lvl_4 %>% unique()

p1 <- SCpubr::do_DimPlot(
  psu,
  group.by = "inherited_celltype_lvl_4",
  label = FALSE,
  legend.position = "none",
  border.size = 1.3,
  pt.size = 2.2,
  colors.use = c(
    "Basal JZ" = c("#F18A85"),
    "Supr. JZ & SD" = c("#FF0000"),
    "Infund. KC" = c("#800080"),
    "Supr. SG" = c("#d7961d"),
    "Basal SD" = c("#FC6A03"),
    "SG progenitors" = c("#d6ccb8"),
    "Basal KC_3" = c("#005f73"),
    "Basal SG" = c("#00B4D8"),
    "Supr. Infund. KC" = c("#ae2012"),
    "Prolif. basal SG" = c("#3CB043"),
    "Matrix" = c("#084943"),
    "Hair shaft" = c("#000000"),
    "IRS" = c("#98BF64"),
    "ORS" = c("#e9d8a6"),
    "Outer bulge" = c("#8B4513"),
    "Inner bulge" = c("#8B80C4")
  )
)

p1

ggsave(
  filename = "./PSU_celltypes.png",
  plot = p1,
  width = 13,
  height = 9,
  dpi = 600
)
