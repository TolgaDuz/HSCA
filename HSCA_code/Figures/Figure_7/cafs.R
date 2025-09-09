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
# Load and inspect HSCA extended fibroblasts
# ===============================

extended_fibro <- readRDS("~/HSCA_data/extended/processed_seurat_objects/Extended_subsets/final/fb_extended.rds")

extended_fibro

DimPlot(extended_fibro, group.by = "inherited_celltype_lvl_5_extended",  label = T, raster = F, repel = T)
DimPlot(extended_fibro, group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_fibro, group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_fibro, group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)
DimPlot(extended_fibro, group.by = "Core")
DimPlot(extended_fibro, group.by = "sample", label = T, repel = F)

DefaultAssay(extended_fibro) <- "RNA"
extended_fibro <- NormalizeData(extended_fibro)

extended_fibro <- subset(extended_fibro, subset = inherited_celltype_lvl_4_extended != "Merkel cell_3")
DimPlot(extended_fibro, group.by = "inherited_celltype_lvl_4_extended",  label = T, raster = T, repel = T)

Idents(extended_fibro) <- extended_fibro$inherited_celltype_lvl_4_extended

new_order_ident <- c(
  "A1",
  "A2",
  "A3",
  "A4",
  "B1",
  "B2",
  "B3",
  "B4",
  "Dermal sheath (C1)",
  "Outer bulge DP (C2)",
  "C3",
  "Anagen DP (C5)",
  "D1",
  "D2",
  "RAMP1+ Fibro (E1)",
  "CAF1",
  "CAF2"
)

# Check consistency between new order and current identities
setdiff(new_order_ident, Idents(extended_fibro))
setdiff(Idents(extended_fibro), new_order_ident)

levels(Idents(extended_fibro))
levels(extended_fibro) <- new_order_ident
levels(Idents(extended_fibro))

genes <- unique(c(
  "RGS5", "NDUFA4L2", "GMFG", "MCAM", "PDGFA", "TAGLN", "COL4A1", "TPM2", "LOXL2",
  "RRAD", "CHI3L1", "TMEM158", "CA12", "WISP1", "BCAT1", "CD82", "TGFBI", "SLC16A3"
))

p <- SCpubr::do_DotPlot(
  sample = extended_fibro,
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

# Create Supp. Fig. X
p2 <- p + ggplot2::scale_fill_gradientn(
  colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  limits = c(0, upper_limit),
  oob = scales::squish
)

p2


# Save Supp. Fig. X
ggsave(
  filename = "./dotplot_cafs.pdf",
  plot = p2,
  device = "pdf",
  width = 10,
  height = 6,
)


# =========================
# Explore CAF clusters by sample origin
# =========================

DimPlot(extended_fibro, group.by = "sample", label = T)
DimPlot(extended_fibro, group.by = "seurat_clusters", label = T)


# The metadata is stored in extended_fibro[[ ]].
# We filter for clusters of interest (74 and 75), which belong to CAF1 and CAF2
df_sub <- df %>% filter(seurat_clusters %in% c(74, 75))

# Count cells per cluster and sample
cluster_counts <- df_sub %>%
  group_by(seurat_clusters, sample) %>%
  summarise(Count = n(), .groups = "drop")

# Add a highlight column to distinguish one sample of interest
# (here: GSM5050540). All other samples are grouped as "Other".
cluster_counts$highlight <- ifelse(cluster_counts$sample == "GSM5050540", "GSM5050540", "Other")

# ==== Bar plot of cluster composition ====
# Create Figure X
p3 <- ggplot(cluster_counts, aes(x = factor(seurat_clusters),
                                 y = Count,
                                 fill = highlight)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("GSM5050540" = "darkblue",
                               "Other" = "grey80")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.title = element_blank()
  ) +
  labs(
    x = "Seurat clusters",
    y = "Number of cells",
    title = "Highlighting sample GSM5050540 in clusters 74 and 75"
  )

# Save figure X
ggsave(
  filename = "./dotplot_cafs_sample_origin.pdf",
  plot = p3,
  device = "pdf",
  width = 8,
  height = 6,
)
