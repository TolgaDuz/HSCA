library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCpubr)

setwd("~/HSCA_code/Figures/Figure_6/")

# Input the HSCA core object
core <- readRDS("/path/to/figshare/Core_atlas.rds")

core$Dataset %>% table() %>% sort(decreasing = TRUE)

# Create Figure 6a
p <- SCpubr::do_FeaturePlot(
  sample = core,
  features = "CCER2"
)

p

# Save Figure 6a
ggsave(
  filename = "~/core_ccer2.png",
  plot = p,
  width = 10,
  height = 9,
  dpi = 600
)

DimPlot(core, group.by = "celltype_lvl_3", label = T)

Idents(core) <- core$inherited_celltype_lvl_3
core$inherited_celltype_lvl_3 %>% table()
core$inherited_celltype_lvl_3 %>% unique()

DefaultAssay(core) <- "RNA"
core <- NormalizeData(core)

DimPlot(core, group.by = "inherited_celltype_lvl_3", label = T)

setdiff(core$inherited_celltype_lvl_3 %>% unique(), new_order_ident)

# ==== streamlined violin plot of Merkel cell gene signature ====

DimPlot(core, group.by = "inherited_celltype_lvl_2", label = T)
Idents(core) <- core$inherited_celltype_lvl_2
core$inherited_celltype_lvl_2 %>% table()

# Initialize new column with values from lvl_2
core$inherited_celltype_lvl_2_modified <- core$inherited_celltype_lvl_2

core$inherited_celltype_lvl_2_modified[
  core$inherited_celltype_lvl_3 == "Merkel cell"
] <- "Merkel cell"

Idents(core) <- core$inherited_celltype_lvl_2_modified
DimPlot(core, group.by = "inherited_celltype_lvl_2_modified", label = T)

table(Idents(core))
DefaultAssay(core) <- "RNA"
core <- NormalizeData(core)

core$inherited_celltype_lvl_2_modified %>% unique()

new_order_ident <- c(
  "Merkel cell",
  "Schwann cell",
  "Neuron",
  "Pilosebaceous Unit",
  "Interfollicular Epidermis",
  "Muscle",
  "Vascular EC",
  "Melanocyte",
  "Fibroblast",
  "Eccrine sweat Gland",
  "Lymphoid",
  "Myeloid",
  "Lymphatic EC",
  "Adipocyte"
)

setdiff(new_order_ident, core$inherited_celltype_lvl_2 %>% unique())

levels(core) <- new_order_ident
Idents(core)

DimPlot(core)

human_colors_list <- c(
  "#cc00c6", # Merkel cell
  "#F18A85", # Schwann cell
  "#ffff7f", # Neuron
  "#000000", # Other
  "#000000",
  "#000000",
  "#000000",
  "#000000",
  "#000000",
  "#000000",
  "#000000",
  "#000000",
  "#000000",
  "#000000"

)

DefaultAssay(core)
core <- JoinLayers(core)

# DGEA
mark <- FindMarkers(core, ident.1 = "Merkel cell", only.pos = T)

mark$gene <- rownames(mark)
head(mark, 15)

FeaturePlot(core, rownames(mark)[1:12])
FeaturePlot(core, rownames(mark)[5:8])

gene_list_plot <- c(
  "CCER2",
  "KRT20",
  "NELL1",
  "CCK",
  "KCNMB2",
  "DACH2",
  "NYAP2",
  "ISL1",
  "ANO5",
  "GNG4",
  "HEPACAM2",
  "SOX2",
  "PCP4",
  "SCGN",
  "VIP"
)

# Create Figure 6b
p1 <- Stacked_VlnPlot(
  seurat_object = core,
  features = gene_list_plot,
  x_lab_rotate = TRUE,
  colors_use = human_colors_list
)

p1


# Save Figure 6b
ggsave(
  filename = "./Figs/violin_merkelcells_main_all_lvl3.pdf",
  plot = p1,
  device = "pdf",
  width = 7,
  height = 8
)


# ==== Enrichment analysis - merkel cells ====
set.seed(42)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(viridis)
library(ggplot2)
library(scales)

DefaultAssay(core)
core <- JoinLayers(core)

Idents(core) <- core$inherited_celltype_lvl_3

# DGEA
mark <- FindMarkers(core, ident.1 = "Merkel cell", only.pos = T)

mark$gene <- rownames(mark)
head(mark, 15)

FeaturePlot(core, rownames(mark)[1:12])
FeaturePlot(core, rownames(mark)[5:8])

# Filter gene list
core_mark <- subset(mark, p_val < 0.05)
genes_to_test <-  core_mark$gene

# Conduct GO Analysis
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db",
                       keyType = "SYMBOL", ont = "BP")

# Create dataframe from GO results
GO_df <- as.data.frame(GO_results)

# Save the dataframe
write.csv(
  GO_df,
  file = paste0("Merkel_cells_GO_results.csv"),
  row.names = FALSE
)

dotplot(GO_results, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment – Merkel cells")

df <- read.csv("./Merkel_cells_GO_results.csv")

names_to_filter <- GO_df$Description[1:10]

# filter dataframe based on names_to_filter
filtered_df <- df %>%
  filter(Description %in% names_to_filter) %>%
  # convert generatio from string to numeric value
  mutate(
    GeneRatio = sapply(
      strsplit(as.character(GeneRatio), "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2])
    )
  ) %>%  # Sort in ascending order
  arrange(GeneRatio)

# calculate minimal adjusted p value
min_p_adjust <- min(filtered_df$p.adjust, na.rm = TRUE)

# scientific formatting of p values
min_p_adjust_label <- format(min_p_adjust, scientific = TRUE, digits = 2)

# Create Figure 6c
plot <- ggplot(
  filtered_df,
  aes(
    x = GeneRatio,
    y = reorder(Description, GeneRatio),
    color = p.adjust,
    size = Count
  )
) +
  geom_point() +
  scale_color_distiller(
    name = paste("p.adjust\n(min: ", min_p_adjust_label, ")", sep = ""),
    palette = "YlGnBu",
    direction = 1,
    trans = "reverse",
    labels = scientific_format(digits = 2)
  ) +
  labs(x = "Gene Ratio", y = "GO Term") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("GO Enrichment Analysis")

plot

plot

# Save Figure 6c
ggsave(
  "./Figs/GO_Enrichment_Analysis_10.pdf",
  plot = plot,
  width = 6.3,
  height = 8
)
