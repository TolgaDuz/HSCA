library(SCpubr)
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
library(CellChat)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

set.seed(42)
setwd("~/HSCA_code/Figures/Figure_4/")

visiumHD_lower_hf <- readRDS("/path/to/figshare/visiumHD_lower_hairfollicle.rds")

Idents(visiumHD_lower_hf) <- visiumHD_lower_hf$celltype
# Inspect
table(Idents(visiumHD_lower_hf))
DimPlot(visiumHD_lower_hf, group.by = "celltype", label = T)
DimPlot(visiumHD_lower_hf, group.by = "seurat_clusters", label = T)
visiumHD_lower_hf$celltype %>% table()

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
  "KRT83- Cortex",
  "KRT83+ Cortex",
  "Cortex/Medulla",
  "Melanocyte",
  "Matrix (Late)",
  "Matrix (Early)",
  "Dermal papilla",
  "Dermal sheath",
  "lower catagen bulge"
)

levels(visiumHD_lower_hf) <- new_order_ident
DimPlot(visiumHD_lower_hf, group.by = "celltype", label = TRUE)

### Adapted from the CellChat tutorial (https://github.com/jinworks/CellChat) ###
### Original authors: Suoqin Jin ###
### Modifications: adapted for Visium HD data ###
DefaultAssay(visiumHD_lower_hf) <- "SCT"
data.input <- visiumHD_lower_hf@assays[["SCT"]]@data

visiumHD_lower_hf$labels <- Idents(visiumHD_lower_hf)
visiumHD_lower_hf$samples <- visiumHD_lower_hf$Sample
meta <- visiumHD_lower_hf[[]]

meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
unique(meta$samples) # check the sample labels

coords1 <- visiumHD_lower_hf@images[["slice1.008um"]]@boundaries[["centroids"]]@coords
coords2 <- visiumHD_lower_hf@images[["slice1.008um.2"]]@boundaries[["centroids"]]@coords

head(coords1, 5)

head(colnames(visiumHD_lower_hf), 5)

# 1. Identify cells per slice
# 2. Retrieve cell names for each slice (the colnames start with the slice name)
cells_slice1 <- grep("^D1", colnames(visiumHD_lower_hf), value = TRUE)
cells_slice2 <- grep("^D2", colnames(visiumHD_lower_hf), value = TRUE)

# 3. Check whether the number of cells matches the number of coordinates
length(cells_slice1) == nrow(coords1)  # sollte TRUE sein
length(cells_slice2) == nrow(coords2)  # sollte TRUE sein

# 4. Insert coordinates
visiumHD_lower_hf@meta.data[cells_slice1, "x"] <- coords1[, "x"]
visiumHD_lower_hf@meta.data[cells_slice1, "y"] <- coords1[, "y"]

visiumHD_lower_hf@meta.data[cells_slice2, "x"] <- coords2[, "x"]
visiumHD_lower_hf@meta.data[cells_slice2, "y"] <- coords2[, "y"]

head(visiumHD_lower_hf@meta.data[
  c(cells_slice1, cells_slice2),
  c("sdimx", "sdimy")
])

spatial.locs <- as.matrix(cbind(
  visiumHD_lower_hf@meta.data[["y"]],
  visiumHD_lower_hf@meta.data[["x"]]
))

spot.size <- 8 # The theoretical spot size (um) in 10X Visium HD
conversion.factor <- spot.size/12.411839512911587 # section D1
conversion.factor <- spot.size/10.438706322430242 # section D2
conversion.factor

rownames(spatial.locs) <- rownames(visiumHD_lower_hf@meta.data)
colnames(spatial.locs) <- c("imagerow", "imagecol")

spatial.factors1 <- data.frame(ratio = conversion.factor, tol = spot.size/2)
spatial.factors2 <- data.frame(ratio = conversion.factor, tol = spot.size/2)

spatial.factors <- rbind(spatial.factors1, spatial.factors2)
rownames(spatial.factors) <- c("D1", "D2")
spatial.factors

# Compute center-to-cetner distance
d.spatial <- computeCellDistance(
  coordinates = spatial.locs,
  ratio = spatial.factors$ratio,
  tol = spatial.factors$tol
)

min(d.spatial[d.spatial != 0])

# Create cellchat obejct
cellchat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "labels",
  datatype = "spatial",
  coordinates = spatial.locs,
  spatial.factors = spatial.factors
)

cellchat@idents
cellchat

CellChatDB <- CellChatDB.human

interaction_input <- CellChatDB$interaction
search <- c(
  "Secreted Signaling",
  "ECM-Receptor",
  # "Cell-Cell Contact",
  "Non-protein Signaling"
)
interaction_input <- interaction_input[interaction_input[["annotation"]] %in%
                                         search, ]

CellChatDB$interaction <- interaction_input
CellChatDB.use <- subsetDB(CellChatDB, non_protein = T)
showDatabaseCategory(CellChatDB.use)

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = T)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.10,
  distance.use = TRUE,
  interaction.range = 250,
  scale.distance = 0.1,
  contact.dependent = TRUE,
  contact.range = 25
)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)

ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  width = 5,
  height = 13,
  font.size = 6,
  color.use =  c(
    "Outer bulge" = "#8B4513",
    "MUCL1+ Inner bulge" = "#8B80C4",
    "CST6+ Inner bulge" = "#C6C6C6",
    "Basal ORS" = "#e9d8a6",
    "Supr. ORS" = "#556b2f",
    "Companion layer" = "#3CA0EC",
    "Henle's layer (upper)" = "#a30000",
    "Henle's layer (middle)" = "#ff2a1f",
    "Henle's layer (lower)" = "#fb8804",
    "Huxley's layer" = "#ffb8b8",
    "Cuticle" = "#D6cf00",
    "LEF1+ Cortex" = "#Fff70f",
    "KRT83+ Cortex" = "#02007a",
    "Cortex/Medulla" = "#E0115F",
    "Melanocyte" = "#FFFFFF",
    "Matrix (Late)" = "#196fb0",
    "Matrix (Early)" = "#084943",
    "Dermal papilla" = "#bfffad",
    "Dermal sheath" = "#94d2db",
    "lower catagen bulge" = "#00fbff"
  ),
  color.heatmap = "YlGnBu",
  cluster.rows = F
)


ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  width = 5,
  height = 13,
  font.size = 6,
  color.use =  c(
    "Outer bulge" = "#8B4513",
    "MUCL1+ Inner bulge" = "#8B80C4",
    "CST6+ Inner bulge" = "#C6C6C6",
    "Basal ORS" = "#e9d8a6",
    "Supr. ORS" = "#556b2f",
    "Companion layer" = "#3CA0EC",
    "Henle's layer (upper)" = "#a30000",
    "Henle's layer (middle)" = "#ff2a1f",
    "Henle's layer (lower)" = "#fb8804",
    "Huxley's layer" = "#ffb8b8",
    "Cuticle" = "#D6cf00",
    "LEF1+ Cortex" = "#Fff70f",
    "KRT83+ Cortex" = "#02007a",
    "Cortex/Medulla" = "#E0115F",
    "Melanocyte" = "#FFFFFF",
    "Matrix (Late)" = "#196fb0",
    "Matrix (Early)" = "#084943",
    "Dermal papilla" = "#bfffad",
    "Dermal sheath" = "#94d2db",
    "lower catagen bulge" = "#00fbff"
  ),
  color.heatmap = "YlGnBu",
  cluster.rows = F
)

# Create Figure 4g
ht1 + ht2

pdf("heatmap_CCC.pdf", width = 9, height = 12)
draw(ht1 + ht2)
dev.off()
