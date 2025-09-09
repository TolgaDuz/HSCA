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
library(viridis)

setwd("~/HSCA_code/Figures/Figure_4/")

visiumHD_lower_hf <- readRDS("/path/to/figshare/visiumHD_lower_hairfollicle.rds")

# Inspect
visiumHD_lower_hf$celltype %>% table()
DimPlot(visiumHD_lower_hf, group.by = "celltype", label = T)
DimPlot(visiumHD_lower_hf, group.by = "seurat_clusters", label = T)

Idents(visiumHD_lower_hf) <- visiumHD_lower_hf$celltype

# Without catagen bulge
visiumHD_lower_hf <- subset(
  visiumHD_lower_hf,
  subset = celltype %in% c("lower catagen bulge"),
  invert = TRUE
)

DimPlot(visiumHD_lower_hf, group.by = "celltype", label = T)

DefaultAssay(visiumHD_lower_hf) <- "Spatial.008um"
visiumHD_lower_hf <- NormalizeData(visiumHD_lower_hf)

# Set order for dotplot visualization
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
  "Dermal sheath"
)

levels(Idents(visiumHD_lower_hf))
levels(visiumHD_lower_hf) <- new_order_ident
levels(Idents(visiumHD_lower_hf))

genes <- c(
  "DIO2", "PDK4",  # Outer bulge
  "MUCL1", "CTSV", # MUCL1 bulge
  "CST6", # CST6+
  "ANGPTL7", "COMP", # Basal ORS
  "SERPINA3",  "FABP5", # Supr. ORS
  "KRT6A", "KRT75", "VEGFA", "CA9", "MEST", "CHI3L1", # Companion layer
  "NES", "GAL", "GABRP", #  Henle's upper
  "CRCT1", "CPA4", "CNFN",  # Henle's middle
  "KRT75", "SLPI", "RHCG", "CAPN8",  # Henle's lower
  "TGFA", "SPINK5", "CDSN", "CLDN4",  # Huxleys layer
  "KRT73", "KRT82", "PROCR",   "CYP26B1", # Cuticle
  "KRT32", "KRT35", "KRT82", "KRT85", # Cuticle
  "SELENBP1", "KRT85",  "KRT85", "LEF1", "DSG4", # LEF1+ cortex
  "KRT83", "NEU2", "ENGASE", "VSIG8", "KRT86", # KRT83+ cortex
  "KRTAP10-2", "SCYGR4", "SCYGR8", # Medulla
  "MLANA", # Melanocyte
  "LPAR6", "SCEL", # Matrix (Late)
  "TOP2A",  # Matrix (Early)
  "COMP", "COL1A1", # Dermal sheath
  "IGFBP3", "DKK2", "RSPO4", "WIF1", "EDN3", "DIO3" # Dermal papilla
)

# Create dotplot for marker genes
p <- SCpubr::do_DotPlot(
  sample = visiumHD_lower_hf,
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

# Create Figure 4d
p2 <- p + ggplot2::scale_fill_gradientn(
  colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  limits = c(0, upper_limit),
  oob = scales::squish
)

p2

# Save Figure 4d
ggsave(
  filename = "./dotplot_lower_hairfollicle.pdf",
  plot = p2,
  device = "pdf",
  width = 20,
  height = 8
)
