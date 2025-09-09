library(Seurat)
library(SoupX)

get_soup_groups <- function(sobj) {
  set.seed(42)
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(object = sobj)
  sobj <- ScaleData(sobj)
  sobj <- RunPCA(sobj, npcs = 30)
  sobj <- FindNeighbors(sobj, dims = 1:30)
  sobj <- FindClusters(sobj)

  sobj@meta.data[["seurat_clusters"]]
}

add_soup_groups <- function(sobj) {
  sobj$soup_group <- get_soup_groups(sobj)
  sobj
}

make_soup <- function(sobj) {
  sample_id <- unique(sobj$sample)
  print(sample_id)
  path <- paste0("./", sample_id, "/outs/raw_feature_bc_matrix/")
  raw <- Read10X(data.dir = path)

  counts_sobj <- GetAssayData(sobj, assay = "RNA", layer = "counts")
  set.seed(42)
  sc <- SoupChannel(raw, counts_sobj)
  sc <- setClusters(sc, sobj$soup_group)
  accession_source <- unique(sobj$Accession_source)

  if (accession_source == "GSE191067" || accession_source == "GSE274955") {
    sc <- setContaminationFraction(sc, contFrac = 0.1, forceAccept = TRUE)
    message("Processing Dataset with contFrac 0.1: ", accession_source)
  } else {
    sc <- autoEstCont(sc, doPlot = TRUE)
  }
  out <- adjustCounts(sc, roundToInt = TRUE)

  # optional: save original counts
  sobj[["original.counts"]] <- CreateAssayObject(counts = counts_sobj)

  sobj <- SetAssayData(sobj, assay = "RNA", slot = "counts", new.data = out)

  sobj
}

process_soupx <- function(data_split) {
  data_split <- lapply(data_split, add_soup_groups)
  data_split <- lapply(data_split, make_soup)
  data_split
}
