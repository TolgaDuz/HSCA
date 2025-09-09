set.seed(42)

setwd("~/HSCA_data/extended/processed_seurat_objects/")

# Load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(mlr)
library(RColorBrewer)
library(devtools)
library(SeuratWrappers)
library(BPCells)

# Import helper functions
source("helpers.R")

# ======================
# Load raw HSCA extended object
# ======================

extended_scvi <- readRDS("./HSCA_extended_raw.rds")

# UMAP visualization by different metadata

DimPlot(extended_scvi, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)
DimPlot(extended_scvi, reduction = "umap", group.by = "Dataset",  label = T, raster = F)
DimPlot(extended_scvi, reduction = "umap", group.by = "anatomical_region_level2",  label = T, raster = F, repel = T)

# Low quality cluster removal
extended_scvi <- subset(extended_scvi, subset = seurat_clusters %in% c(26,45,107,56,52,81), invert = T)
# Remove samples with strong wound healing signatures
extended_scvi <- subset(extended_scvi, subset = sample %in% c("GSM8238438", "GSM8238437"), invert = T)
extended_scvi

DimPlot(extended_scvi, reduction = "umap", group.by = "seurat_clusters", label = T, raster = F)

# ===============================================
# Load Seurat objects for different extension subsets
# ===============================================

extended_immu_scvi <- readRDS("./Extended_subsets/raw/immune_extended.rds")
extended_ec_scvi <- readRDS("./Extended_subsets/raw/ec_extended.rds")
extended_muscle_scvi <- readRDS("./Extended_subsets/raw/muscle_extended.rds")
extended_fibro_scvi <- readRDS("./Extended_subsets/raw/fb_extended.rds")
extended_esg_scvi <- readRDS("./Extended_subsets/raw/esg_extended.rds")
extended_IFE_scvi <- readRDS("./Extended_subsets/raw/ife_extended.rds")
extended_psu_scvi <- readRDS("./Extended_subsets/raw/psu_extended.rds")
extended_other_scvi <- readRDS("./Extended_subsets/raw/other_cells_extended.rds")

# ==== Initial visualization of UMAPs by cell type ====
DimPlot(extended_immu_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_ec_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_muscle_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_fibro_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(extended_esg_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_IFE_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_other_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)


# ===============================================
# Clean numeric-only level 5 cell types
# ===============================================

extended_immu_scvi <- clean_lvl5_numeric(extended_immu_scvi)
extended_ec_scvi <- clean_lvl5_numeric(extended_ec_scvi)
extended_muscle_scvi <- clean_lvl5_numeric(extended_muscle_scvi)
extended_fibro_scvi <- clean_lvl5_numeric(extended_fibro_scvi)
extended_psu_scvi <- clean_lvl5_numeric(extended_psu_scvi)
extended_esg_scvi <- clean_lvl5_numeric(extended_esg_scvi)
extended_IFE_scvi <- clean_lvl5_numeric(extended_IFE_scvi)
extended_other_scvi <- clean_lvl5_numeric(extended_other_scvi)


# Check cleaned metadata
extended_immu_scvi$inherited_celltype_lvl_4_extended %>% table()
extended_immu_scvi$inherited_celltype_lvl_5_extended %>% table()
extended_muscle_scvi$inherited_celltype_lvl_5_extended %>% table()
extended_psu_scvi$inherited_celltype_lvl_5_extended %>% table()

# ===============================================
# Remove low-quality (LQ) cells
# ===============================================

extended_immu_scvi <- remove_LQ_cells(extended_immu_scvi)
extended_ec_scvi <- remove_LQ_cells(extended_ec_scvi)
extended_muscle_scvi <- remove_LQ_cells(extended_muscle_scvi)
extended_fibro_scvi <- remove_LQ_cells(extended_fibro_scvi)
extended_psu_scvi <- remove_LQ_cells(extended_psu_scvi)
extended_esg_scvi <- remove_LQ_cells(extended_esg_scvi)
extended_IFE_scvi <- remove_LQ_cells(extended_IFE_scvi)
extended_other_scvi <- remove_LQ_cells(extended_other_scvi)

# ===============================================
# Inherit level 4 annotations to level 5 where missing
# ===============================================

extended_immu_scvi <- inherit_lvl4_to_lvl5(extended_immu_scvi)
extended_ec_scvi <- inherit_lvl4_to_lvl5(extended_ec_scvi)
extended_muscle_scvi <- inherit_lvl4_to_lvl5(extended_muscle_scvi)
extended_fibro_scvi <- inherit_lvl4_to_lvl5(extended_fibro_scvi)
extended_psu_scvi <- inherit_lvl4_to_lvl5(extended_psu_scvi)
extended_esg_scvi <- inherit_lvl4_to_lvl5(extended_esg_scvi)
extended_IFE_scvi <- inherit_lvl4_to_lvl5(extended_IFE_scvi)
extended_other_scvi <- inherit_lvl4_to_lvl5(extended_other_scvi)

extended_psu_scvi$inherited_celltype_lvl_5_extended <- ifelse(
  extended_psu_scvi$inherited_celltype_lvl_5_extended == "Lower isthmus",
  "Basal JZ_4",
  extended_psu_scvi$inherited_celltype_lvl_5_extended
)

# ==== Visualize level 5 annotations ====
DimPlot(extended_psu_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(extended_ec_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_esg_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_muscle_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_immu_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_fibro_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_IFE_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)
DimPlot(extended_other_scvi, group.by = "inherited_celltype_lvl_5_extended", label = T)

# ==== Combine Seurat objects into a named list =====
seurat_objects <- list(
  "fb" = extended_fibro_scvi,
  "immu" = extended_immu_scvi,
  "muscle" = extended_muscle_scvi,
  "psu" = extended_psu_scvi,
  "EC" = extended_ec_scvi,
  "ESG" = extended_esg_scvi,
  "ife" = extended_IFE_scvi,
  "other" = extended_other_scvi
)

seurat_objects

# ===============================================
# Define hierarchical mapping of cell types
# Level 5 -> Level 1 Bottom-up mapping
# ===============================================

hierarchy <- list(
  "A1_4" = c("A1", "Fibro A", "Fibroblast", "Stroma"),
  "A2_4" = c("A2", "Fibro A", "Fibroblast", "Stroma"),
  "A3_4" = c("A3", "Fibro A", "Fibroblast", "Stroma"),
  "A4_4" = c("A4", "Fibro A", "Fibroblast", "Stroma"),
  "PDZRN4+ FB" = c("A1", "Fibro A", "Fibroblast", "Stroma"),
  "B1_4" = c("B1", "Fibro B", "Fibroblast", "Stroma"),
  "B2_4" = c("B2", "Fibro B", "Fibroblast", "Stroma"),
  "B3_4" = c("B3", "Fibro B", "Fibroblast", "Stroma"),
  "B4_4" = c("B4", "Fibro B", "Fibroblast", "Stroma"),
  "D1_4" = c("D1", "Fibro D", "Fibroblast", "Stroma"),
  "D2_4" = c("D2", "Fibro D", "Fibroblast", "Stroma"),
  "Dermal sheath (C1)_4" = c("Dermal sheath (C1)", "Fibro C", "Fibroblast", "Stroma"),
  "Outer bulge DP (C2)_4" = c("Outer bulge DP (C2)", "Fibro C", "Fibroblast", "Stroma"),
  "C3_4" = c("C3", "Fibro C", "Fibroblast", "Stroma"),
  "Anagen DP (C5)_4" = c("Anagen DP (C5)", "Fibro C", "Fibroblast", "Stroma"),
  "RAMP1+ Fibro (E1)_4" = c("RAMP1+ Fibro (E1)", "Fibro E", "Fibroblast", "Stroma"),
  "CAF1_4" = c("CAF1", "Fibro diseased", "Fibroblast", "Stroma"),
  "CAF2_4" = c("CAF2", "Fibro diseased", "Fibroblast", "Stroma"),

  "CD8+ T cell_4" = c("CD8+ T cell", "T cell", "Lymphoid", "Immune"),
  "naive CD8+ T cell" = c("CD8+ T cell", "T cell", "Lymphoid", "Immune"),
  "GD-T cell_4" = c("GD-T cell", "T cell", "Lymphoid", "Immune"),
  "CD4+ T cell_4" = c("CD4+ T cell", "T cell", "Lymphoid", "Immune"),
  "naive CD4+ T cell" = c("CD4+ T cell", "T cell", "Lymphoid", "Immune"),
  "Reg. T cell_4" = c("Reg. T cell", "T cell", "Lymphoid", "Immune"),
  "Prolif. T cell_4" = c("Prolif. T cell", "T cell", "Lymphoid", "Immune"),
  "Naive B cell_4" = c("Naive B cell", "B cell", "Lymphoid", "Immune"),
  "Plasma cell_4" = c("Plasma cell", "B cell", "Lymphoid", "Immune"),
  "FGFBP2+ NK_4" = c("FGFBP2+ NK", "NK", "Lymphoid", "Immune"),
  "SPINK2+ NK" = c("XCL2+ NK", "NK", "Lymphoid", "Immune"),
  "GZMK+ NK" = c("XCL2+ NK", "NK", "Lymphoid", "Immune"),

  "LC_4" = c("LC", "DC", "Myeloid", "Immune"),
  "Prolif. DC_4" = c("Prolif. DC", "DC", "Myeloid", "Immune"),
  "cDC1" = c("cDC", "DC", "Myeloid", "Immune"),
  "cDC2" = c("cDC", "DC", "Myeloid", "Immune"),
  "pDC" = c("cDC", "DC", "Myeloid", "Immune"),
  "Mature DC_4" = c("Mature DC", "DC", "Myeloid", "Immune"),
  "Inflammatory DC_4" = c("Inflammatory DC", "DC", "Myeloid", "Immune"),
  "Anti-inflammatory Mph_4" = c("Anti-inflammatory Mph", "Mph", "Myeloid", "Immune"),
  "Inflammatory Mph_4" = c("Inflammatory Mph", "Mph", "Myeloid", "Immune"),
  "C3+ Mph" = c("TREM2+ Mph", "Mph", "Myeloid", "Immune"),
  "LPL+ Mph" = c("TREM2+ Mph", "Mph", "Myeloid", "Immune"),
  "EREG+ Mono" = c("cMono", "Monocyte", "Myeloid", "Immune"),
  "cMono_4" = c("cMono", "Monocyte", "Myeloid", "Immune"),
  "ncMono_4" = c("ncMono", "Monocyte", "Myeloid", "Immune"),
  "Monocyte_3" = c("Monocyte_3", "Monocyte", "Myeloid", "Immune"),
  "Neutrophil_3" = c("Neutrophil_3", "Neutrophil", "Myeloid", "Immune"),
  "Mast cell_3" = c("Mast cell_3", "Mast cell", "Myeloid", "Immune"),

  "RSG5+ EC_4" = c("RSG5+ EC", "Venous EC", "Vascular EC", "Endothelial"),
  "Venous 1 EC_4" = c("Venous 1 EC", "Venous EC", "Vascular EC", "Endothelial"),
  "Venous 2 EC_4" = c("Venous 2 EC", "Venous EC", "Vascular EC", "Endothelial"),
  "Capillary EC_3" = c("Capillary EC_3", "Capillary EC", "Vascular EC", "Endothelial"),
  "Arterial EC_3" = c("Arterial EC_3", "Arterial EC", "Vascular EC", "Endothelial"),
  "SCG3+ LEC_4" = c("SCG3+ LEC", "Lymphatic EC", "Lymphatic EC", "Endothelial"),
  "NEO1+ LEC_4" = c("NEO1+ LEC", "Lymphatic EC", "Lymphatic EC", "Endothelial"),
  "LYVE1+ LEC_4" = c("LYVE1+ LEC", "Lymphatic EC", "Lymphatic EC", "Endothelial"),

  "Skeletal muscle_3" = c("Skeletal muscle_3", "Skeletal muscle", "Muscle", "Stroma"),
  "Muscle progenitor_4" = c("Muscle progenitor", "SMC", "Muscle", "Stroma"),
  "DES+ SMC_4" = c("DES+ SMC", "SMC", "Muscle", "Stroma"),
  "RERGL+ SMC_4" = c("RERGL+ SMC", "SMC", "Muscle", "Stroma"),
  "CYP26B1+ SMC" = c("STEAP4+ SMC", "SMC", "Muscle", "Stroma"),
  "FRMD3+ SMC" = c("STEAP4+ SMC", "SMC", "Muscle", "Stroma"),
  "TM4SF1+ SMC" = c("STEAP4+ SMC", "SMC", "Muscle", "Stroma"),

  "Dark cell_4" = c("Dark cell", "Coil", "Eccrine sweat Gland", "Cutaneous Epithelial"),
  "Myoepithelial_4" = c("Myoepithelial", "Coil", "Eccrine sweat Gland", "Cutaneous Epithelial"),
  "Clear cell 1_4" = c("Clear cell 1", "Coil", "Eccrine sweat Gland", "Cutaneous Epithelial"),
  "Clear cell 2_4" = c("Clear cell 2", "Coil", "Eccrine sweat Gland", "Cutaneous Epithelial"),
  "S100P+ luminal_4" = c("S100P+ luminal", "Duct", "Eccrine sweat Gland", "Cutaneous Epithelial"),
  "S100P- luminal_4" = c("S100P- luminal", "Duct", "Eccrine sweat Gland", "Cutaneous Epithelial"),
  "Basal luminal_4" = c("Basal luminal", "Duct", "Eccrine sweat Gland", "Cutaneous Epithelial"),

  "Infund. KC_4" = c("Infund. KC", "Infundibulum", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Supr. Infund. KC_4" = c("Supr. Infund. KC", "Infundibulum", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Basal SD_4" = c("Basal SD", "Isthmus", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Basal JZ_4" = c("Basal JZ", "Isthmus", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Supr. JZ & SD_4" = c("Supr. JZ & SD", "Isthmus", "Pilosebaceous Unit", "Cutaneous Epithelial"),

  "SEB-2L" = c("Supr. SG", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "SEB-2" = c("Supr. SG", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "SEB-1" = c("Supr. SG", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "SEB-T" = c("Supr. SG", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "SEB-B" = c("Basal SG", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Prolif. SEB-B" = c("Basal SG", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "SG progenitors_4" = c("SG progenitors", "SG", "Pilosebaceous Unit", "Cutaneous Epithelial"),

  "Matrix (Early)" = c("Matrix", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Matrix (Late) & Henle's layer (lower)" = c("Matrix", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "LEF1+ Cortex" = c("Hair shaft", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "KRT83+ Cortex" = c("Hair shaft", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Cuticle hair shaft" = c("Hair shaft", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Cuticle IRS" = c("IRS", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Henle's layer (upper)" = c("IRS", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Companion layer" = c("ORS", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Basal ORS" = c("ORS", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Supr. ORS" = c("ORS", "Bulb", "Pilosebaceous Unit", "Cutaneous Epithelial"),

  "SHG" = c("Outer bulge", "Bulge", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "Outer bulge" = c("Outer bulge", "Bulge", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "WFDC3+ Inner bulge" = c("Inner bulge", "Bulge", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "CST6+ Inner bulge" = c("Inner bulge", "Bulge", "Pilosebaceous Unit", "Cutaneous Epithelial"),
  "MUCL1+ Inner bulge" = c("Inner bulge", "Bulge", "Pilosebaceous Unit", "Cutaneous Epithelial"),

  "Basal KC_3" = c("Basal KC_3", "Basal KC", "Interfollicular Epidermis", "Cutaneous Epithelial"),
  "Prolif. KC_3" = c("Prolif. KC_3", "Prolif. KC", "Interfollicular Epidermis", "Cutaneous Epithelial"),
  "Spinous KC_3" = c("Spinous KC_3", "Spinous KC", "Interfollicular Epidermis", "Cutaneous Epithelial"),
  "Granular KC_3" = c("Granular KC_3", "Granular KC", "Interfollicular Epidermis", "Cutaneous Epithelial"),
  "Cornified KC_3" = c("Cornified KC_3", "Cornified KC", "Interfollicular Epidermis", "Cutaneous Epithelial"),
  "Merkel cell_3" = c("Merkel cell_3", "Merkel cell", "Interfollicular Epidermis", "Cutaneous Epithelial"),

  "Adipocyte_2" = c("Adipocyte_2", "Adipocyte_2", "Adipocyte", "Stroma"),
  "Melanocyte_2" = c("Melanocyte_2", "Melanocyte_2", "Melanocyte", "Neural lineage"),
  "Neuron_2" = c("Neuron_2", "Neuron_2", "Neuron", "Neural lineage"),
  "Neural progenitor_2" = c("Neural progenitor_2", "Neural progenitor_2", "Neural progenitor", "Neural lineage"),
  "Schwann cell_2" = c("Schwann cell_2", "Schwann cell_2", "Schwann cell", "Neural lineage"),
  "CNTNAP2+ Neuron_4" = c("CNTNAP2+ Neuron", "Sensory neuron", "Neuron", "Neural lineage"),

  "Chondrocyte_2" = c("Chondrocyte_2", "Chondrocyte_2", "Chondrocyte", "Other"),
  "Erythrocyte_2" = c("Erythrocyte_2", "Erythrocyte_2", "Erythrocyte", "Other"),
  "Schwann cell (Ex vivo)_2" = c("Schwann cell (Ex vivo)_2", "Schwann cell (Ex vivo)_2", "Schwann cell (Ex vivo)", "Other")
)

extended_esg_scvi$inherited_celltype_lvl_4_extended %>% table()

DimPlot(extended_immu_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(extended_esg_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T)
DimPlot(extended_muscle_scvi, group.by = "inherited_celltype_lvl_4_extended", label = T)

# ===============================================
# Update metadata for all objects based on hierarchy
# Propagate lvl5 -> lvl4 -> lvl3 -> lvl2 -> lvl1
# ===============================================

for (obj_name in names(seurat_objects)) {
  print(obj_name)
  seurat_obj <- seurat_objects[[obj_name]]

  if ("inherited_celltype_lvl_5_extended" %in% colnames(seurat_obj@meta.data)) {
    print("inherited_celltype_lvl_5_extended")
    seurat_obj@meta.data[["inherited_celltype_lvl_4_extended"]] <- unname(
      sapply(seurat_obj@meta.data[["inherited_celltype_lvl_5_extended"]], function(x)
        ifelse(x %in% names(hierarchy), hierarchy[[x]][1], NA))
    )
  }


  if ("inherited_celltype_lvl_4_extended" %in% colnames(seurat_obj@meta.data)) {
    print("inherited_celltype_lvl_4_extended")
    seurat_obj@meta.data[["inherited_celltype_lvl_3_extended"]] <- unname(
      sapply(seurat_obj@meta.data[["inherited_celltype_lvl_5_extended"]], function(x)
        ifelse(x %in% names(hierarchy), hierarchy[[x]][2], NA))
    )
  }

  if ("inherited_celltype_lvl_3_extended" %in% colnames(seurat_obj@meta.data)) {
    print("inherited_celltype_lvl_3_extended")
    seurat_obj@meta.data[["inherited_celltype_lvl_2_extended"]] <- unname(
      sapply(seurat_obj@meta.data[["inherited_celltype_lvl_5_extended"]], function(x)
        ifelse(x %in% names(hierarchy), hierarchy[[x]][3], NA))
    )
  }

  if ("inherited_celltype_lvl_2_extended" %in% colnames(seurat_obj@meta.data)) {
    print("inherited_celltype_lvl_2_extended")
    seurat_obj@meta.data[["inherited_celltype_lvl_1_extended"]] <- unname(
      sapply(seurat_obj@meta.data[["inherited_celltype_lvl_5_extended"]], function(x)
        ifelse(x %in% names(hierarchy), hierarchy[[x]][4], NA))
    )
  }

  seurat_objects[[obj_name]] <- seurat_obj
}

print(names(seurat_objects))

seurat_obj$inherited_celltype_lvl_1_extended

# ==== Visualize cell type hierarchy across subsets ====
DimPlot(seurat_objects[["psu"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["psu"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(seurat_objects[["fb"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["fb"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(seurat_objects[["immu"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["immu"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(seurat_objects[["muscle"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["muscle"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(seurat_objects[["EC"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["EC"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(seurat_objects[["ESG"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["ESG"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)
DimPlot(seurat_objects[["ife"]], group.by = "inherited_celltype_lvl_1_extended", label = T, raster = F)
DimPlot(seurat_objects[["ife"]], group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(seurat_objects[["other"]], group.by = "inherited_celltype_lvl_4_extended", label = T, raster = F)
DimPlot(seurat_objects[["other"]], group.by = "inherited_celltype_lvl_5_extended", label = T, raster = F)

# ==== Save filtered Seurat objects ====
saveRDS(seurat_objects[["psu"]], "./Extended_subsets/final/psu_extended.rds")
saveRDS(seurat_objects[["fb"]], "./Extended_subsets/final/fb_extended.rds")
saveRDS(seurat_objects[["immu"]], "./Extended_subsets/final/immune_extended.rds")
saveRDS(seurat_objects[["muscle"]], "./Extended_subsets/final/muscle_extended.rds")
saveRDS(seurat_objects[["EC"]], "./Extended_subsets/final/ec_extended.rds")
saveRDS(seurat_objects[["ESG"]], "./Extended_subsets/final/esg_extended.rds")
saveRDS(seurat_objects[["ife"]], "./Extended_subsets/final/ife_extended.rds")
saveRDS(seurat_objects[["other"]], "./Extended_subsets/final/other_cells_extended.rds")


# ===============================================
# Transfer the cell type labels from the individual objects to the
# merged final HSCA core (do this for each cell type level)
# ===============================================

for (obj_name in names(seurat_objects)) {
  print(paste("Compare:", obj_name))

  seurat_obj <- seurat_objects[[obj_name]]

  # Identify common cells
  common_cells <- intersect(rownames(seurat_obj@meta.data), rownames(extended_scvi@meta.data))
  print(paste("Common cells found:", length(common_cells)))

  # Transfer cell type levels if shared cells exist
  if (length(common_cells) > 0) {
    for (level in 1:5) {
      level_col <- paste0("inherited_celltype_lvl_", level, "_extended")

      if (level_col %in% colnames(seurat_obj@meta.data)) {
        extended_scvi@meta.data[common_cells, level_col] <-
          seurat_obj@meta.data[common_cells, level_col]
      }
    }
  }
}


sum(is.na(extended_scvi@meta.data$inherited_celltype_lvl_1_extended))

DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_1_extended", raster = T, label = T)
DimPlot(extended_scvi, group.by = "Dataset", raster = F, label =T)

# Keep only cells with assigned labels (low-quality/unlabeled cells removed)
extended_scvi <- subset(extended_scvi, subset = !is.na(inherited_celltype_lvl_1_extended))

extended_scvi
sum(is.na(extended_scvi$inherited_celltype_lvl_5_extended))
DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_1_extended", raster = T, label = T)
DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_2_extended", raster = T, label = T)
DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_3_extended", raster = T, label = T)
DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_4_extended", raster = T, label = T)
DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_5_extended", raster = T, label = T)

extended_scvi$inherited_celltype_lvl_4_extended %>% table()

set.seed(8)
extended_scvi <- RunUMAP(extended_scvi, reduction = "X_scvi_emb",reduction.name = "umap", return.model = T, dims = 1:10)

DimPlot(extended_scvi, group.by = "inherited_celltype_lvl_1_extended", raster = T, label = T)

saveRDS(extended_scvi, "./final_atlas/HSCA_extended.rds")
