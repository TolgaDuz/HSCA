library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(hdf5r)
library(ggplot2)
library(gridExtra)
library(pals)
library(data.table)
library(tidyr)
library(harmony)

setwd("~/HSCA_code/Figures/Figure_4/")

visiumHD_full <- readRDS("/path/to/figshare/VisiumHD_full.rds")

# Inspect
DimPlot(visiumHD_full, group.by = "celltype", label = TRUE)
DimPlot(visiumHD_full, group.by = "seurat_clusters", label = TRUE)

Idents(visiumHD_full) <- "seurat_clusters"

# ==== check spatial location of clusters ====
cells <- CellsByIdentities(visiumHD_full, idents = c("7"))

SpatialDimPlot(
  visiumHD_full,
  images = "slice1.008um",
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = TRUE, combine = TRUE
) + NoLegend()


SpatialDimPlot(
  visiumHD_full,
  images = "slice1.008um.2",
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = TRUE, combine = TRUE
) + NoLegend()


# subset to lower hair follicle cell types
Idents(visiumHD_full) <- "celltype"
lower_hairfollicle <- subset(visiumHD_full, idents = c(
  "ORS",
  "Companion layer",
  "lower catagen bulge",
  "Outer bulge",
  "KRT6A+ Inner bulge",
  "CST6+ Inner bulge",
  "MUCL1+ Inner bulge",
  "Dermal sheath",
  "IRS",
  "Matrix",
  "Hair shaft"
))

lower_hairfollicle@reductions <- list()
lower_hairfollicle@assays[["BANKSY"]] <- NULL

Idents(lower_hairfollicle) <- "seurat_clusters"

lower_hairfollicle
DefaultAssay(lower_hairfollicle) <- "Spatial.008um"

# Process Seurat object
lower_hairfollicle <- SCTransform(lower_hairfollicle, assay = "Spatial.008um")
set.seed(42)

# ==== Prepare data for multi-sample Banksy analysis ====
coords1 <- lower_hairfollicle@images[["slice1.008um"]]@boundaries[["centroids"]]@coords
coords2 <- lower_hairfollicle@images[["slice1.008um.2"]]@boundaries[["centroids"]]@coords

head(coords1, 5)

head(colnames(lower_hairfollicle), 5)

# Retrieve cell names for each slice (the column names start with the slice name)
cells_slice1 <- grep("^D1", colnames(lower_hairfollicle), value = TRUE)
cells_slice2 <- grep("^D2", colnames(lower_hairfollicle), value = TRUE)

# Check whether the number of cells matches the number of coordinates
length(cells_slice1) == nrow(coords1)  # must equal TRUE
length(cells_slice2) == nrow(coords2)  # must euqal TRUE

# Insert coordinates
lower_hairfollicle@meta.data[cells_slice1, "sdimx"] <- coords1[, "x"]
lower_hairfollicle@meta.data[cells_slice1, "sdimy"] <- coords1[, "y"]

lower_hairfollicle@meta.data[cells_slice2, "sdimx"] <- coords2[, "x"]
lower_hairfollicle@meta.data[cells_slice2, "sdimy"] <- coords2[, "y"]

head(lower_hairfollicle@meta.data[c(cells_slice1, cells_slice2), c("sdimx", "sdimy")])

# Run Banksy embedding
lower_hairfollicle_banksy <- RunBanksy(
  lower_hairfollicle,
  lambda = 0.15,
  assay = "SCT",
  slot = "data",
  dimx = "sdimx",
  dimy = "sdimy", # column names in metadata for coord x and coord y
  features = "variable",
  group = "Sample",
  split.scale = TRUE,
  k_geom = 5
)

DefaultAssay(lower_hairfollicle_banksy) <- "BANKSY"

Idents(lower_hairfollicle_banksy) <- lower_hairfollicle_banksy$Sample

# ==== Dimensionality reduction and batch correction ====
# PCA on Banksy embedding
lower_hairfollicle_banksy <- RunPCA(
  lower_hairfollicle_banksy,
  assay = "BANKSY",
  reduction.name = "pca.banksy",
  features = rownames(lower_hairfollicle_banksy),
  npcs = 30
)

# Harmony integration across samples
set.seed(42)
lower_hairfollicle_har <- RunHarmony(
  lower_hairfollicle_banksy,
  reduction.use = "pca.banksy",
  group.by.vars = "Sample"
)

# UMAP visualization
set.seed(42)
lower_hairfollicle_har <- RunUMAP(
  lower_hairfollicle_har,
  reduction = "harmony",
  dims = 1:30
)

# Neighbor graph
lower_hairfollicle_har <- FindNeighbors(
  lower_hairfollicle_har,
  reduction = "harmony",
  dims = 1:30
)

# Clustering
set.seed(42)
DefaultAssay(lower_hairfollicle_har) <- "BANKSY"
lower_hairfollicle_har <- FindClusters(
  lower_hairfollicle_har,
  cluster.name = "banksy_cluster",
  resolution = 6.5
)

DimPlot(lower_hairfollicle_har, group.by = "banksy_cluster", label = TRUE)
DimPlot(lower_hairfollicle_har, group.by = "celltype", label = TRUE)

DefaultAssay(lower_hairfollicle_har) <- "SCT"

Idents(visiumHD_full) <- "seurat_clusters"
FeaturePlot(lower_hairfollicle_har, c("percent_mito", "nFeature_Spatial.008um"))

# remove LQ and VEC
lower_hairfollicle_har <- subset(lower_hairfollicle_har, idents = c(44, 24), invert = TRUE)
DimPlot(lower_hairfollicle_har, group.by = "banksy_cluster", label = TRUE)
DimPlot(lower_hairfollicle_har, group.by = "celltype", label = TRUE)

# Huxleys layer
FeaturePlot(lower_hairfollicle_har, c("TGFA", "CDSN", "SPINK5", "CLDN4"))
# Henle's layer (lower)
FeaturePlot(lower_hairfollicle_har, c("SLPI", "RHCG", "KRT71", "CAPN8"))
# Henle's layer (middle)
FeaturePlot(lower_hairfollicle_har, c("CRCT1", "CPA4", "CNFN", "ASPRV1"))
# Henle's layer (upper)
# Caution: comp and selenop also in ORS
FeaturePlot(lower_hairfollicle_har, c("CHI3L1", "KRT75", "COMP", "SELENOP"))
# Late Matrix
FeaturePlot(lower_hairfollicle_har, c("LPAR6", "MYCN", "FABP5", "SCEL"))
# Early Matrix
FeaturePlot(lower_hairfollicle_har, c("HIST1H1B", "TOP2A", "MYCN", "HIST1H1D"))
# Medulla
FeaturePlot(lower_hairfollicle_har, c("KRTAP10-2", "SCYGR4", "KRTAP12-1", "SCYGR8"))
FeaturePlot(lower_hairfollicle_har, c("VSIG8", "KRT86", "KRT83", "KRTAP5-8"))
# Cortex
FeaturePlot(lower_hairfollicle_har, c("SELENBP1", "KRT85", "LY6G6D", "GPNMB")) # LEF1+ Cortex
FeaturePlot(lower_hairfollicle_har, c("LEF1")) # LEF1+ Cortex
FeaturePlot(lower_hairfollicle_har, c("KRT81", "KRT83", "KRT86")) # KRT83+ Cortex
# Melanocytes
FeaturePlot(lower_hairfollicle_har, c("TYRP1", "MLANA", "TYR", "QPCT"))
# DPC
FeaturePlot(lower_hairfollicle_har, c("IGFBP3", "DKK2", "RSPO4", "WIF1"))
# Cuticle
# IRS cuticle
FeaturePlot(lower_hairfollicle_har, c("KRT25", "KRT26", "KRT27", "KRT28", "KRT71", "KRT72", "KRT73"))
# Cuticle hair shaft
FeaturePlot(lower_hairfollicle_har, c("KRT82", "PROCR", "VSIG8", "ACTBL2"))
FeaturePlot(lower_hairfollicle_har, c("KRT32", "SOX21", "CYP26B1", "PPP2R1B", "PADI3"))
FeaturePlot(lower_hairfollicle_har, c("KRT32", "KRT35", "KRT39", "KRT40", "KRT82", "KRT85"))
# CL
FeaturePlot(lower_hairfollicle_har, c("GAL", "GABRP", "SERPINA1", "CALB2"))
FeaturePlot(lower_hairfollicle_har, c("KRT6A", "SH3D21", "NES", "ADAMTS1"))
FeaturePlot(lower_hairfollicle_har, c("VEGFA", "CA9", "MEST", "ADAMTS1"))
# Basal ORS
FeaturePlot(lower_hairfollicle_har, c("LGR5", "COMP", "SELENOP", "ANGPTL7"))
FeaturePlot(lower_hairfollicle_har, c("FGF18", "CRLF1", "ALAD", "RNF152"))
# Bulge
FeaturePlot(lower_hairfollicle_har, c("KRT6A", "KRT16", "S100A2", "SERPINA3"))
FeaturePlot(lower_hairfollicle_har, c("DSG3", "DSP", "PKP1", "FABP5"))
FeaturePlot(lower_hairfollicle_har, c("CST6", "KRT6A"))
FeaturePlot(lower_hairfollicle_har, c("MUCL1"))
DFeaturePlot(lower_hairfollicle_har, c("FABP5"))
FeaturePlot(lower_hairfollicle_har, c("DIO2", "PDK4", "CXCL14"))
FeaturePlot(lower_hairfollicle_har, c("COL6A1"))
FeaturePlot(lower_hairfollicle_har, c("ANGPTL7"))
# lower catagen
FeaturePlot(lower_hairfollicle_har, c("SULT1E1", "NMB", "THBS1", "CLU", "PAPLN"))
# Dermal Sheath
FeaturePlot(lower_hairfollicle_har, c("COL1A2", "COL11A1", "SPARC", "MGP"))
# SHG (not present in visium HD dataset)
FeaturePlot(lower_hairfollicle_har, c(
  "KREMEN2", "SLC2A3", "EPCAM", "TNF", "SOX11", "CASC15",
  "SPRY4", "SPX5", "ROBO1", "CUX2", "RARB"
))

DefaultAssay(lower_hairfollicle_har) <- "Spatial.008um"
lower_hairfollicle_har <- JoinLayers(lower_hairfollicle_har)

# ==== DGEA with SCT ====
lower_hairfollicle@assays$SCT@SCTModel.list

# Update UMI assay name and set SCT as default
slot(object = lower_hairfollicle_har@assays$SCT@SCTModel.list[[2]], name = "umi.assay") <- "Spatial.008um"
DefaultAssay(lower_hairfollicle_har) <- "SCT"

# Run differential gene expression
lower_hairfollicle_har <- PrepSCTFindMarkers(lower_hairfollicle_har)
mark_spatial <- FindMarkers(
  lower_hairfollicle_har,
  ident.1 = c("31"),
  only.pos = TRUE,
  assay = "SCT"
)
head(mark_spatial, 15)

FeaturePlot(lower_hairfollicle_har, rownames(mark_spatial)[1:12])

# ==== DGEA with log-normalized data ====m
DefaultAssay(lower_hairfollicle_har) <- "Spatial.008um"
lower_hairfollicle_har <- NormalizeData(lower_hairfollicle_har)
mark_spatial <- FindMarkers(
  lower_hairfollicle_har,
  ident.1 = c("76"),
  only.pos = TRUE,
  assay = "Spatial.008um"
)
head(mark_spatial, 15)
FeaturePlot(lower_hairfollicle_har, rownames(mark_spatial)[1:12])
###

DefaultAssay(lower_hairfollicle_har) <- "nFeature_Spatial.008um"
lower_hairfollicle_har <- NormalizeData(lower_hairfollicle_har)

DimPlot(lower_hairfollicle_har.subset, group.by = "banksy_cluster", label = TRUE)
DefaultAssay(lower_hairfollicle_har) <- "SCT"

DimPlot(lower_hairfollicle_har, label = TRUE)
FeaturePlot(lower_hairfollicle_har, c(rownames(markers)[1:12]), raster = FALSE)

# ==== Assign cell types ====
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(63, 56), "Huxley's layer", lower_hairfollicle_har$seurat_clusters)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(41), "Henle's layer (lower)", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(64), "Henle's layer (middle)", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(55, 19, 8), "Henle's layer (uppper)", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(32), "Melanocyte", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(65), "Cortex/Medulla", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(40), "Cuticle", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(33), "LEF1+ Cortex", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(60), "KRT83+ Cortex", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(45, 15, 12, 70, 9, 7, 2, 68, 42, 37), "CST6+ Inner bulge", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(66,21,13,35), "MUCL1+ Inner bulge", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(51), "WFDC3+ Inner bulge", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(26, 22, 6, 46, 39, 17, 38,73, 20, 11), "Basal ORS", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(4,52, 23,49,48,34,31,43), "Supr. ORS", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(0,36,47,30,3,62,54,27,5,71,10,25,28, 72), "Outer bulge", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(61), "Matrix (Late)", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(16), "Matrix (Early)", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(1,53,26,14,50,67), "Dermal sheath", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(69), "Dermal papilla", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(18,29,59,57), "Companion layer", lower_hairfollicle_har$celltype)
lower_hairfollicle_har$celltype <- ifelse(lower_hairfollicle_har$seurat_clusters %in% c(58), "lower catagen bulge", lower_hairfollicle_har$celltype)

DimPlot(lower_hairfollicle_har, group.by = "celltype", label = TRUE)

saveRDS(lower_hairfollicle_har, "./visiumHD_lower_hairfollicle.rds")
