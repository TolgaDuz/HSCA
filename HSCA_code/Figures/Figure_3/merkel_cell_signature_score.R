set.seed(42)

library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(SCpubr)
options(warn = 1)

setwd("~/HSCA_code/Figures/Figure_3/")

visiumHD_full <- readRDS("/path/to/figshare/VisiumHD_full.rds")

DimPlot(visiumHD_full, group.by = "celltype")

signature_full <- c(
  "CCER2", "KRT20", "NELL1", "CCK", "KCNMB2", "DACH2", "NYAP2", "ISL1", "ANO5",
  "GNG4", "HEPACAM2", "SOX2", "PCP4", "SCGN", "VIP"
)

signature_full <- c(
  "CCER2", "KRT20", "NELL1", "CCK", "KCNMB2", "DACH2", "NYAP2", "ISL1", "ANO5",
  "HEPACAM2", "SOX2", "SCGN", "VIP"
)

DefaultAssay(visiumHD_full) <- "SCT"

visiumHD_full <- AddModuleScore(
  object = visiumHD_full,
  features = list(signature_full_pcp4_gng4),
  name = "signature_full_pcp4_gng4"
)

DefaultAssay(visiumHD_full) <- "Spatial.008um"
visiumHD_full <- NormalizeData(visiumHD_full)

visiumHD_full <- AddModuleScore(
  object = visiumHD_full,
  features = list(signature_full_pcp4_gng4),
  name = "signature_full_pcp4_gng4"
)

colnames(visiumHD_full@meta.data)

FeaturePlot(visiumHD_full, "signature_full1")
FeaturePlot(visiumHD_full, "CCER2")

# For Supp. Fig 7
saveRDS(visiumHD_full, "visiumHD_full_merkel_sig.rds")
