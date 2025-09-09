# setting working dir
setwd("~/HSCA_code/Figures/Figure_7")

set.seed(42)

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SCpubr)
library(SeuratWrappers)
library(BPCells)
library(lintr)

# Configure reticulate to use scvi conda environment
reticulate::use_condaenv(
  "scvi",
  conda = "/opt/conda/condabin/conda",
  required = TRUE
)

# ===============================
# Load and inspect HSCA extended atlas
# ===============================

extended_atlas <- readRDS("~/HSCA_data/extended/processed_seurat_objects/final_atlas/HSCA_extended.rds")

DimPlot(extended_atlas, group.by = "inherited_celltype_lvl_1_extended", label = TRUE)

extended_atlas

# ==== Summary statistics of the atlas ====

table(extended_atlas$Core)

length(unique(extended_atlas$sample))        # Number of distinct samples
length(unique(extended_atlas$subject_ID))    # Number of distinct donors

length(unique(extended_atlas$inherited_celltype_lvl_5_extended))
length(unique(extended_atlas$inherited_celltype_lvl_5[extended_atlas$Core == "Yes"]))

# Count unique samples and datasets split by Core vs Extended-only
length(unique(extended_atlas$sample[extended_atlas$Core == "Yes"]))
length(unique(extended_atlas$sample[extended_atlas$Core == "No"]))
length(unique(extended_atlas$Dataset[extended_atlas$Core == "Yes"]))
length(unique(extended_atlas$Dataset[extended_atlas$Core == "No"]))
length(unique(extended_atlas$subject_ID[extended_atlas$Core == "Yes"]))
length(unique(extended_atlas$subject_ID[extended_atlas$Core == "No"]))


# =====================================================
# Visualize HSCA extended atlas by cell type, laber transfer uncert. and dataset
# =====================================================

extended_atlas$inherited_celltype_lvl_2 %>% table()

# Extract unique cell types at level 2 from the HSCA core
celltypes_lvl2 <- extended_atlas$inherited_celltype_lvl_2 %>% unique()
celltypes_lvl2
celltypes_clean <- celltypes_lvl2[celltypes_lvl2 != "Extended"]
celltypes_clean <- celltypes_lvl2

celltypes_clean

Idents(extended_atlas) <- extended_atlas$inherited_celltype_lvl_2

table(Idents(extended_atlas), useNA = "always")

sum(is.na(extended_atlas$inherited_celltype_lvl_2))

Idents(extended_atlas) %>% table()
Idents(extended_atlas) %>% unique()


# ==== Highlighting the HSCA core cells on top of the extended datasets ====

p <- SCpubr::do_DimPlot(
  sample = extended_atlas,
  idents.keep = celltypes_clean,
  label = FALSE,
  legend.position = "none",
  border.size = 3,
  pt.size = 0.2,
  na.value = "#363636" # Color for extended datasets
)

p

# Define a custom color scheme for each level 2 cell type
# Create Figure 7b
p2 <- p +
  scale_color_manual(values =  c(
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
  ))


p2

# Save the UMAP highlighting core vs extended cell types (Figure 7b)
ggsave(
  filename = "./hsca_extended_core_highlight.png",
  plot = p2,
  width = 12,
  height = 9,
  dpi = 600
)

# ==== FeaturePlot showing label transfer uncertainty ====

# Create Figure 7c
p <- SCpubr::do_FeaturePlot(
  sample = extended_atlas,
  features = "inherited_celltype_lvl_3_transfer_uncert"
)

p

# Save Figure 7c
ggsave(
  filename = "./HSCA_extended_umap_uncerts.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 600
)


# ==== HSCA extended colored by dataset ====

p <- SCpubr::do_DimPlot(
  sample = extended_atlas,
  group.by = "Dataset",
  label = FALSE,
  legend.position = "none",
  border.size = 3,
  pt.size = 0.2
)

p

ggsave(
  filename = "./HSCA_extended_umap_datasets.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 600
)
