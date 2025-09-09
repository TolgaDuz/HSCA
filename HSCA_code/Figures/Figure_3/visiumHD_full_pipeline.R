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

spaceranger_out_D1 <- "~/HSCA_data/VisiumHD/Skin/D1/outs/"
spaceranger_out_D2 <- "~/HSCA_data/VisiumHD/Skin/D2/outs/"

# Inspect output files
list.files(spaceranger_out_D1)

object_D1 <- Load10X_Spatial(data.dir = spaceranger_out_D1)
object_D2 <- Load10X_Spatial(data.dir = spaceranger_out_D2)

DefaultAssay(object_D1) <- "Spatial.008um"
DefaultAssay(object_D2) <- "Spatial.008um"

object_D1$Sample <- "D1"
object_D2$Sample <- "D2"

object_D1[["percent_mito"]] <- PercentageFeatureSet(object_D1, pattern = "^MT-")
object_D2[["percent_mito"]] <- PercentageFeatureSet(object_D2, pattern = "^MT-")

# Inspect mitochondrial read content across the tissue
SpatialFeaturePlot(object_D1, features = "percent_mito") +
  theme(legend.position = "right")

seurat_merge <- merge(object_D1, y = object_D2, add.cell.ids = c("D1", "D2"))

# ==== QC filtering ====
object_D1.qc <- subset(object_D1, subset = nFeature_Spatial.008um > 100 & percent_mito < 25)
object_D2.qc <- subset(object_D2, subset = nFeature_Spatial.008um > 100 & percent_mito < 25)

# Inspect
object_D1.qc
object_D2.qc

SpatialFeaturePlot(object_D1.qc, features = "percent_mito") + theme(legend.position = "right")
SpatialFeaturePlot(object_D1.qc, features = "nFeature_Spatial.008um") + theme(legend.position = "right")

# Merge to one seurat object
seurat_merge <- merge(
  object_D1.qc,
  y = object_D2.qc,
  add.cell.ids = c("D1", "D2")
)

DefaultAssay(seurat_merge) <- "Spatial.008um"

# Process Seurat object
seurat_merge <- SCTransform(seurat_merge, assay = "Spatial.008um")
set.seed(42)

# ==== Prepare data for multi-sample Banksy analysis ====
set.seed(42)
coords1 <- seurat_merge@images[["slice1.008um"]]@boundaries[["centroids"]]@coords
coords2 <- seurat_merge@images[["slice1.008um.2"]]@boundaries[["centroids"]]@coords

head(coords1, 5)

head(colnames(seurat_merge), 5)
# Retrieve cell names for each slice (the column names start with the slice name)
cells_slice1 <- grep("^D1", colnames(seurat_merge), value = TRUE)
cells_slice2 <- grep("^D2", colnames(seurat_merge), value = TRUE)

# Check whether the number of cells matches the number of coordinates
length(cells_slice1) == nrow(coords1)  # must equal TRUE
length(cells_slice2) == nrow(coords2)  # must equal TRUE

# Insert coordinates
seurat_merge@meta.data[cells_slice1, "sdimx"] <- coords1[, "x"]
seurat_merge@meta.data[cells_slice1, "sdimy"] <- coords1[, "y"]

seurat_merge@meta.data[cells_slice2, "sdimx"] <- coords2[, "x"]
seurat_merge@meta.data[cells_slice2, "sdimy"] <- coords2[, "y"]

head(seurat_merge@meta.data[c(cells_slice1, cells_slice2), c("sdimx", "sdimy")])

# Run Banksy embedding
DefaultAssay(seurat_merge) <- "SCT"
seurat_merge <- RunBanksy(
  seurat_merge, lambda = 0.2,
  assay = "SCT",
  slot = "data",
  dimx = "sdimx",
  dimy = "sdimy", # column names in metadata for coord x and coord y
  features = "variable",
  group = "Sample",
  split.scale = TRUE,
  k_geom = 15
)

DefaultAssay(seurat_merge) <- "BANKSY"
Idents(seurat_merge) <- seurat_merge$Sample

# ==== Dimensionality reduction and batch correction ====
# PCA on Banksy embedding
seurat_merge <- RunPCA(
  seurat_merge,
  assay = "BANKSY",
  reduction.name = "pca.banksy",
  features = rownames(seurat_merge),
  npcs = 30
)

# Harmony integration across samples
seurat_merge_har <- RunHarmony(
  seurat_merge,
  reduction.use = "pca.banksy",
  group.by.vars = "Sample"
)

# UMAP visualization
set.seed(42)
seurat_merge_har <- RunUMAP(
  seurat_merge_har,
  reduction = "harmony",
  dims = 1:30
)

# Neighbor graph
seurat_merge_har <- FindNeighbors(
  seurat_merge_har,
  reduction = "harmony",
  dims = 1:30
)

# Clustering
set.seed(42)
DefaultAssay(seurat_merge_har) <- "BANKSY"
seurat_merge_har <- FindClusters(
  seurat_merge_har,
  cluster.name = "banksy_cluster",
  resolution = 2.8
)

DimPlot(seurat_merge_har, group.by = "banksy_cluster", label = TRUE)
DimPlot(seurat_merge_har, group.by = "seurat_clusters", label = TRUE)

Idents(seurat_merge_har) <- "banksy_cluster"

DefaultAssay(seurat_merge_har) <- "Spatial.008um"

seurat_merge_har <- JoinLayers(seurat_merge_har)

# ==== DGEA with SCT ====
seurat_merge@assays$SCT@SCTModel.list

# Update UMI assay name and set SCT as default
slot(object = seurat_merge_har@assays$SCT@SCTModel.list[[2]], name = "umi.assay") <- "Spatial.008um"
DefaultAssay(seurat_merge_har) <- "SCT"

# Run differential gene expression
seurat_merge_har <- PrepSCTFindMarkers(seurat_merge_har)
mark_spatial <- FindMarkers(
  seurat_merge_har,
  ident.1 = c("57"),
  only.pos = TRUE,
  assay = "SCT"
)
head(mark_spatial, 15)

FeaturePlot(seurat_merge_har, rownames(mark_spatial)[1:12])


# ==== DGEA with log-normalized data ====
DefaultAssay(seurat_merge_har) <- "Spatial.008um"
seurat_merge_har <- NormalizeData(seurat_merge_har)
mark_spatial <- FindMarkers(
  seurat_merge_har,
  ident.1 = c("36"),
  only.pos = TRUE,
  assay = "Spatial.008um"
)
head(mark_spatial, 15)

FeaturePlot(seurat_merge_har, rownames(mark_spatial)[1:12])

# Additional plots
SpatialFeaturePlot(seurat_merge_har, features = "FASN") + theme(legend.position = "right")
FeaturePlot(seurat_merge_har, features = c("S100A8")) + theme(legend.position = "right")
DimPlot(seurat_merge_har, group.by = "Sample", label = TRUE)

Idents(seurat_merge_har) <- "seurat_clusters"

# Inspect cell clusters spatially
cells <- CellsByIdentities(seurat_merge_har, idents = c("41"))

SpatialDimPlot(
  seurat_merge_har,
  images = "slice1.008um",
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = TRUE, combine = TRUE
) + NoLegend()


SpatialDimPlot(
  seurat_merge_har,
  images = "slice1.008um.2",
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = TRUE, combine = TRUE
) + NoLegend()

DefaultAssay(seurat_merge_har) <- "SCT"
FeaturePlot(seurat_merge_har, features = c("percent_mito", "nFeature_Spatial.008um"))
FeaturePlot(seurat_merge_har, features = c("TOB1", "MMP7", "IFI27"))
FeaturePlot(seurat_merge_har, features = c("APOD", "DCN", "CFD")) # APOD+ FB
FeaturePlot(seurat_merge_har, features = c("TNXB", "FBN1", "DCN")) # DCN+ FB
FeaturePlot(seurat_merge_har, features = c("COL1A1", "COL3A1", "COL1A2", "COL6A3")) # COL1A1+ FB
FeaturePlot(seurat_merge_har, features = c("COL6A2", "COL6A1", "SPARC", "MMP2")) # COL1A1+ FB
FeaturePlot(seurat_merge_har, features = c("FABP4", "PLIN4", "PLIN1", "ADIPOQ")) # Adipocyte
FeaturePlot(seurat_merge_har, features = c("CD74", "VIM", "TMSB10", "CST3")) # Immune Cells
FeaturePlot(seurat_merge_har, features = c("CCL19", "CXCR4", "LYZ", "CLU")) # Immune Cells
FeaturePlot(seurat_merge_har, features = c("KRT5", "PTN", "C1QTNF12", "LMO1", "ANGPTL2")) # Isthmus
FeaturePlot(seurat_merge_har, features = c("TYRP1", "DCT", "COL17A1", "MLANA")) # Melanocyte
FeaturePlot(seurat_merge_har, features = c("KRT10", "KRT1", "AQP3", "KRTDAP")) # Spinous KC
FeaturePlot(seurat_merge_har, features = c("KRT15", "KRT5", "COL17A1", "POSTN")) # Basal KC
FeaturePlot(seurat_merge_har, features = c("S100A8", "S100A9", "KRT10", "KRTDAP")) # Infund. KC
FeaturePlot(seurat_merge_har, features = c("KRTDAP", "LOR", "FLG", "CALML5")) # Granular KC
FeaturePlot(seurat_merge_har, features = c("CST6")) # CST6+ Inner Bulge
FeaturePlot(seurat_merge_har, features = c("MUCL1", "CTSV", "KRT16", "SPINK5")) # MUCL1+ Inner Bulge
FeaturePlot(seurat_merge_har, features = c("KRT6A", "KRT16", "KRT6B", "SERPINA3")) # MUCL1+ Inner Bulge
FeaturePlot(seurat_merge_har, features = c("PDK4", "CXCL14", "DIO2", "DST")) # Outer Bulge
FeaturePlot(seurat_merge_har, features = c("TCHH", "KRT71", "KRT25", "KRT27")) # IRS
FeaturePlot(seurat_merge_har, features = c("TCHH", "KRT71", "KRT25", "KRT27")) # IRS
FeaturePlot(seurat_merge_har, features = c("S100A3", "CDSN", "PSORS1C2", "SPINK5", "KRT73")) # IRS
FeaturePlot(seurat_merge_har, features = c("ID3", "GPNMB", "KRT35", "MYCN", "HIST1H1B", "TOP2A")) # IRS
FeaturePlot(seurat_merge_har, features = c("KRT31", "KRT33A", "KRT33B", "KRT34",  "KRT36")) # CORTEX
FeaturePlot(seurat_merge_har, features = c("KRT37", "KRT38", "KRT39", "KRT81", "KRT83", "KRT86")) # CORTEX
FeaturePlot(seurat_merge_har, features = c("KRT7", "WFDC2", "IL1R2", "NNAT", "PPARG", "FASN", "PLIN5")) # CORTEX

# remove LQ
seurat_merge_har_filt <- subset(seurat_merge_har, subset = seurat_clusters %in% c(57, 42, 44), invert = TRUE)
DimPlot(seurat_merge_har_filt, label = TRUE)

# ==== Assign cell types ====
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(26), "SMC", seurat_merge_har_filt$seurat_clusters)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(27), "APM", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(13, 4), "Endothelial", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(47, 51, 3, 9), "ESG coil", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(19, 2), "ESG duct", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(45), "Papillary FB", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(6), "COL1A1+ FB", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(58), "Adipocyte", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(14, 50), "LMT cell", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(5), "Isthmus", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(32), "Melanocyte", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(15), "Basal KC", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(0), "Early spinous KC", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(1), "Late spinous KC", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(8), "Inflamed IFE", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(55), "Granular KC", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(30, 29, 33), "Infundibular KC", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(20, 34, 18), "CST6+ Inner bulge", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(43, 28), "MUCL1+ Inner bulge", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(25, 21), "KRT6A+ Inner bulge", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(16, 17, 24, 31), "ORS", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(11, 12), "Outer bulge", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(38, 40, 49), "IRS", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(46, 52), "Hair shaft", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(41), "Matrix", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(35), "Bursted seb.", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(36), "Late SEB-2.", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(10), "SEB-1", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(39, 22), "SEB-B", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(23), "Sebaceous duct & JZ", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(53), "activated inflamed FB", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(48), "SG associated FB", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(7, 54), "Dermal sheath", seurat_merge_har_filt$celltype)
seurat_merge_har_filt$celltype <- ifelse(seurat_merge_har_filt$seurat_clusters %in% c(56), "lower catagen bulge", seurat_merge_har_filt$celltype)
seurat_merge$celltype <- ifelse(seurat_merge$seurat_clusters %in% c(37), "Companion layer", seurat_merge$celltype)

DimPlot(seurat_merge_har_filt, group.by = "celltype", label = TRUE)

saveRDS(seurat_merge_har_filt, "./VisiumHD_full.rds")
