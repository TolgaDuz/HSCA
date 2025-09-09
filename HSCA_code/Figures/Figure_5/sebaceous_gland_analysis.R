library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCpubr)
library(scCustomize)

setwd("~/HSCA_code/Figures/Figure_5/")

sg_core <- readRDS("/path/to/figshare/sg_core.rds")

# Remove low quality cluster
sg_core <- subset(sg_core, subset = inherited_celltype_lvl_5 != "LQ")

# Inspect object
DimPlot(sg_core, group.by = "inherited_celltype_lvl_5", label = T)
DimPlot(sg_core, group.by = "seurat_clusters", label = T)
DimPlot(sg_core, group.by = "inherited_celltype_lvl_5",
        label = T, reduction = "phate")

sg_core$inherited_celltype_lvl_5 %>% table()

# Create Figure 5e or Supp. Fig. 5
p1 <- SCpubr::do_DimPlot(
  sg_core,
  group.by = "inherited_celltype_lvl_5",
  reduction = "phate", # Or change to umap
  label = F,
  border.size = 1.7,
  pt.size = 2,
  label.box = T,
  legend.position = "none",
  repel = T,
  colors.use = c(
    "Basal JZ_4" = c("#F18A85"),
    "SG progenitors_4" = c("#d6ccb8"),
    "Prolif. SEB-B" = c("#3CB043"),
    "SEB-B" = c("#00B4D8"),
    "SEB-T" = c("#8B80C4"),
    "SEB-1" = c("#d7961d"),
    "SEB-2" = c("#B22222"),
    "SEB-2L" = c("#D6cf00")
  )
)

p1

# Save Figure 5e or Supp. Fig. 5
ggsave(
  filename = "./Figs/SG_phate.pdf",
  plot = p1,
  width = 13,
  height = 9
)


# ==== Accession source pie chart ====
# Subset: exclude Basal JZ and SG progenitors
sg_core_sub <- subset(
  sg_core,
  subset = inherited_celltype_lvl_5 %in% c("Basal JZ_4", "SG progenitors_4"),
  invert = T
)

sg_core_sub

# 1. Count cells per accession source and sort by size
df <- sg_core_sub@meta.data %>%
  count(Accession_source) %>%
  arrange(desc(n)) %>%
  mutate(
    prop = n / sum(n),
    label = paste0(Accession_source, "\n(", round(prop * 100, 1), "%)")
  )

# 2. Create pie chart (Fugure 5g)
p2 <- ggplot(
  df,
  aes(x = "", y = prop, fill = factor(Accession_source, levels = Accession_source))) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Cell distribution by accession sourc", fill = "Dataset") +
  scale_fill_brewer(palette = "Set3", labels = df$label) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

p2

# Save Figure 5g
ggsave(
  filename = "./sg_pie_chart_accession_source.pdf",
  plot = p2,
  width = 13,
  height = 9
)

# ==== celltype pie chart ====

celltype_column <- "inherited_celltype_lvl_5"

# remove Basal JZ cells since they are not of interest here
sg_core_wo_jz <- subset(
  sg_core,
  subset = inherited_celltype_lvl_5 != "Basal JZ_4"
)

# Create dataframe with counts per cell type
df <- sg_core_wo_jz@meta.data %>%
  count(!!sym(celltype_column)) %>%      # count cells per type
  arrange(desc(n)) %>%                   # order by cell number
  mutate(
    celltype = !!sym(celltype_column),   # store cell type name
    label = paste0(celltype, "\n(n = ", n, ")")  # add label with counts
  )

celltype_colors <- c(
  "SG progenitors_4"  = "#d6ccb8",
  "Prolif. SEB-B"      = "#3CB043",
  "SEB-B"              = "#00B4D8",
  "SEB-T"              = "#8B80C4",
  "SEB-1"              = "#d7961d",
  "SEB-2"              = "#B22222",
  "SEB-2L"             = "#FFE130"
)

# Create pie chart (Figure 5f)
p <- ggplot(df, aes(x = "", y = n, fill = factor(celltype, levels = celltype))) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Zellverteilung nach Zelltyp", fill = "Zelltyp") +
  scale_fill_manual(values = celltype_colors, labels = df$label) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

p

# Save Figure 5f
ggsave(
  filename = "./sg_pie_chart_cellcount.pdf",
  plot = p,
  width = 13,
  height = 9
)


# ==== Inspect SG differentiation markers ====
# Basal SD
FeaturePlot(sg_core, c("KRT15", "PTN", "COL17A1", "C1QTNF12"))
FeaturePlot(sg_core, c("KRT5", "SERPINB4", "COL17A1", "POSTN"))

# Supr. SD
FeaturePlot(sg_core, c("KRT6A", "IVL", "KRT79", "KRT1"))
FeaturePlot(sg_core, c("KRT6A", "IVL", "KRT79", "KRT1"), reduction = "phate")

# SEB-B
FeaturePlot(sg_core, c("NNAT", "IL1R2", "MKI67", "TINAGL1"))
FeaturePlot(sg_core, c("NNAT", "IL1R2", "KRT7", "TINAGL1"), reduction = "phate")

# SEB-T
FeaturePlot(sg_core, c("WFDC2", "KRT7", "TBC1D4", "MKI67"))
FeaturePlot(sg_core, c("WFDC2", "KRT7", "TBC1D4", "MKI67"), reduction = "phate")

# SEB-1
FeaturePlot(sg_core, c("HAO2", "ACO1", "HSD11B1", "FASN"))
FeaturePlot(sg_core, c("HAO2", "ACO1", "HSD11B1", "FASN"), reduction = "phate")

# SEB-2
FeaturePlot(sg_core, c("PNPLA3", "PLIN5", "PPARG", "FADS2"))
FeaturePlot(sg_core, c("PNPLA3", "PLIN5", "PPARG", "FADS2"), reduction = "phate")

# SEB-2L
FeaturePlot(sg_core, c("CRAT", "SEC14L6", "DNASE1L2")) #
FeaturePlot(sg_core, c("CRAT", "SEC14L6", "DNASE1L2"), reduction = "phate")

DimPlot(sg_core, group.by = "seurat_clusters", label = T, reduction = "phate")
DimPlot(sg_core, group.by = "seurat_clusters", label = T, reduction = "umap")
DimPlot(sg_core, group.by = "Dataset", label = T)


# Find new marker genes for refined signature (Supp. Fig 6)
sg_core <- JoinLayers(sg_core)
Idents(sg_core) <- sg_core$inherited_celltype_lvl_5
Idents(sg_core) %>% table()
markers <- FindMarkers(object = sg_core, ident.1 = "SEB-1",  only.pos = TRUE)
head(markers, 40)

FeaturePlot(sg_core, c(rownames(markers)[1:12]), raster = F)
FeaturePlot(sg_core, c(rownames(markers)[13:24]), raster = F)
FeaturePlot(sg_core, c(rownames(markers)[25:36]), raster = F)

sg_core$inherited_celltype_lvl_5 %>% table()

# === Violinplots SG differentiation markers====

DefaultAssay(sg_core) <- "RNA"
sg_core <- NormalizeData(sg_core)

Idents(sg_core) <- sg_core$inherited_celltype_lvl_5
DimPlot(sg_core, group.by = "inherited_celltype_lvl_5", label = T)

new_order_ident <- c(
  "Basal JZ_4",
  "SG progenitors_4",
  "Prolif. SEB-B",
  "SEB-B",
  "SEB-T",
  "SEB-1",
  "SEB-2",
  "SEB-2L"
)

levels(sg_core) <- new_order_ident
DimPlot(sg_core)

human_colors_list <- c(
  "#F18A85",
  "#d6ccb8",
  "#3CB043",
  "#00B4D8",
  "#8B80C4",
  "#d7961d",
  "#B22222",
  "#FFE130"
)

DefaultAssay(sg_core)

gene_list_plot <- c(
  "NNAT",
  "IL1R2",
  "TINAGL1",
  "WFDC2",
  "KRT7",
  "TBC1D4",
  "HAO2",
  "ACO1",
  "PNPLA3",
  "PLIN5",
  "SEC14L6",
  "CRAT",
  "DNASE1L2"
)

# Progenitors
gene_list_plot <- c(
  "KRT15",
  "PTN",
  "C1QTNF12"
)

# Create Figure 5h
p1 <- Stacked_VlnPlot(
  seurat_object = sg_core,
  features = gene_list_plot,
  x_lab_rotate = TRUE,
  colors_use = human_colors_list
)

p1

# Save Figure 5h
ggsave(
  filename = "./sg_progenitors_vln.pdf",
  plot = p1,
  device = "pdf",
  width = 7,
  height = 8
)


# ==== refined SG differentiation signature ====
new_order_ident <- c(
  "Basal JZ_4",
  "SG progenitors_4",
  "Prolif. SEB-B",
  "SEB-B",
  "SEB-T",
  "SEB-1",
  "SEB-2",
  "SEB-2L"
)

levels(Idents(sg_core))
levels(sg_core) <- new_order_ident
levels(Idents(sg_core))

genes <- c(
  "NNAT",
  "IL1R2",
  "TINAGL1",
  "WFDC2",
  "KRT7",
  "TBC1D4",
  "HAO2",
  "ACO1",
  "PNPLA3",
  "PLIN5",
  "SEC14L6",
  "CRAT",
  "DNASE1L2",
  "ACACB",
  "ACACA",
  "PPARG",
  "FASN",
  "ELOVL5",
  "ACSBG1",
  "IDH1",
  "FAR2",
  "TF",
  "PLEKHH1",
  "MOGAT1",
  "PMFBP1",
  "LY86-AS1",
  "PDE6A",
  "ENTPD1",
  "ROS1",
  "SCHLAP1",
  "AZGP1",
  "FADS1",
  "HGD",
  "PM20D1",
  "GPD1",
  "FA2H",
  "GAL",
  "SOD2-OT1",
  "AZGP1",
  "AP002387.2",
  "C3orf84",
  "ACAD8",
  "TERB2",
  "GPX3",
  "CLPSL2",
  "GJC3"
)

p <- SCpubr::do_DotPlot(
  sample = sg_core,
  features = genes,
  cluster = F,
  dot.scale = 12,
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
p[["layers"]][[1]][["geom"]][["default_aes"]][["stroke"]] <- 0.1

# Everything that is not upregulated (Z-score ≤ 0) is hidden
# Set negative Z-scores to NA
p$data$Avg.Exp[p$data$Avg.Exp <= 0] <- NA

# Also set P.Exp (dot size) to NA if the corresponding expression is NA
p$data$P.Exp[is.na(p$data$Avg.Exp)] <- NA

p

# Determine the upper limit of Avg.Exp for color scaling
upper_limit <- max(p$data$Avg.Exp, na.rm = TRUE)
upper_limit

# Create Supp. Fig 7
p2 <- p + ggplot2::scale_fill_gradientn(
  colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  limits = c(0, upper_limit),
  oob = scales::squish
)

p2

# Save Supp. Fig 7
ggsave(
  filename = "./dotplot_SG_new.pdf",
  plot = p2,
  device = "pdf",
  width = 20,
  height = 5
)
