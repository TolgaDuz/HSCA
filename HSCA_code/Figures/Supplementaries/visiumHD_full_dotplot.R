set.seed(42)

library(Seurat)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(SCpubr)
library(viridis)

setwd("~/HSCA_code/Figures/Figure_2/")

visiumHD_full <- readRDS("/path/to/figshare/VisiumHD_full.rds")

# Inspect
visiumHD_full$Sample %>% table()
visiumHD_full$celltype %>% table()

# Mean and maximum number of features in Sample D1
mean(visiumHD_full$nFeature_Spatial.008um[visiumHD_full$Sample == "D1"])
max(visiumHD_full$nFeature_Spatial.008um[visiumHD_full$Sample == "D1"])

# Mean and maximum number of features in Sample D2
mean(visiumHD_full$nFeature_Spatial.008um[visiumHD_full$Sample == "D2"])
max(visiumHD_full$nFeature_Spatial.008um[visiumHD_full$Sample == "D2"])

DimPlot(visiumHD_full, group.by = "seurat_clusters", label = TRUE)
DimPlot(visiumHD_full, group.by = "celltype", label = TRUE)

# Aggregate cell types
visiumHD_full$celltype <- ifelse(visiumHD_full@meta.data$seurat_clusters %in% c(56), "Outer bulge", visiumHD_full@meta.data$celltype)
visiumHD_full$celltype <- ifelse(visiumHD_full@meta.data$celltype == "Henle's layer (upper)", "IRS", visiumHD_full@meta.data$celltype)
visiumHD_full$celltype <- ifelse(visiumHD_full@meta.data$celltype == "Companion layer", "ORS", visiumHD_full@meta.data$celltype)

Idents(visiumHD_full) <- visiumHD_full$celltype

DefaultAssay(visiumHD_full) <- "Spatial.008um"
visiumHD_full <- NormalizeData(visiumHD_full)

Idents(visiumHD_full) %>% table()

new_order_ident <- c(
  "Infundibular KC",
  "Isthmus",
  "Sebaceous duct & JZ",
  "Outer bulge",
  "CST6+ Inner bulge",
  "MUCL1+ Inner bulge",
  "ORS",
  "KRT6A+ Inner bulge",
  "IRS",
  "Hair shaft",
  "Matrix",
  "SEB-B",
  "SEB-1",
  "Late SEB-2.",
  "Bursted seb.",
  "Basal KC",
  "Early spinous KC",
  "Late spinous KC",
  "Granular KC",
  "Inflamed IFE",
  "Dermal sheath",
  "COL1A1+ FB",
  "activated inflamed FB",
  "FB & immune cell (periglandular)",
  "Papillary FB",
  "SMC",
  "APM",
  "Adipocyte",
  "Endothelial",
  "ESG coil",
  "ESG duct",
  "Melanocyte",
  "LMT cell"
)

levels(Idents(visiumHD_full))
levels(visiumHD_full) <- new_order_ident
levels(Idents(visiumHD_full))

# Check consistency between new order and current identities
setdiff(new_order_ident, levels(visiumHD_full))
setdiff(levels(visiumHD_full), new_order_ident)

DimPlot(visiumHD_full, group.by = "celltype", label = TRUE)


# ==== Conduct DGEA to find marker genes ====
DefaultAssay(visiumHD_full) <- "SCT"

visiumHD_full$celltype %>% table()

visiumHD_full <- PrepSCTFindMarkers(visiumHD_full)

mark_spatial <- FindMarkers(
  visiumHD_full,
  ident.1 = c("Adipocyte"),
  only.pos = TRUE,
  assay = "SCT"
)

head(mark_spatial, 45)

FeaturePlot(visiumHD_full, rownames(mark_spatial)[1:12])
FeaturePlot(visiumHD_full, rownames(mark_spatial)[13:24])
FeaturePlot(visiumHD_full, c("percent_mito", "nFeature_Spatial.008um"))

# Marker genes for skin cell types
genes <- c(
  "S100A8", "S100A9", "PTN",  "C1QTNF12", "GATA6",
  "DIO2", "CXCL14", "KRT6B",
  "CST6", "MUCL1", "CTSV", "COMP", "SELENOP",
  "KRT6A", "KRT75",
  "TCHH", "KRT71", "KRT25",
  "KRT35", "DSG4", "HIST1H1B",  "TOP2A",
  "FASN",
  "PLIN5",
  "SEC14L6",
  "KRT15", "KRT1",
  "KRTDAP", "FLG", "AQP3",
  "SFN", "FOSL1", "COL11A1", "COMP", "APOD", "DCN",
  "FOSB", "IRF1",
  "SOCS3",  "PTGDS",
  "C3",
  "MYL9", "DES", "CNN1",
  "PLIN4", "ADIPOQ",
  "CLDN5", "AQP1",
  "SCGB2A2", "DCD", "MMP7",
  "TYRP1", "DCT",
  "CD74", "CCL19"
)

genes <- unique(genes)

DefaultAssay(visiumHD_full)

p <- SCpubr::do_DotPlot(
  sample = visiumHD_full,
  features = genes,
  cluster = FALSE,
  dot.scale = 9,
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

ggsave(
  filename = "./visiumHD_full_dotplot.pdf",
  plot = p2,
  device = "pdf",
  width = 22,
  height = 10
)
