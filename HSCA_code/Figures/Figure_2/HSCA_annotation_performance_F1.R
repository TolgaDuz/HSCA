library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SCpubr)
library(RColorBrewer)

setwd("~/HSCA_code/Figures/Figure_2/")

# Input the HSCA core object (available on Figshare)
core <- readRDS("/path/to/figshare/HSCA_core.rds")

# Inspect
core$Dataset %>% table() %>% sort(decreasing = TRUE)
core$inherited_celltype_lvl_5 %>% unique()

# Remove the inherited suffix from celltypes
core$goldstandard_celltype_lvl_3 <-
  gsub("_[0-9]+$", "", core$inherited_celltype_lvl_3)

# Inspect
core$orig_celltype_lvl_3 %>% table()
core$goldstandard_celltype_lvl_3 %>% table()

# Harmonize cell type names between the original dataset annotations
# and the atlas re-annotations
core$orig_celltype_lvl_3 <- ifelse(
  core$orig_celltype_lvl_3 == "Neuron_CNTNAP2+",
  "Sensory neuron",
  core$orig_celltype_lvl_3
)

core$orig_celltype_lvl_3 <- ifelse(
  core$orig_celltype_lvl_3 == "NK cell",
  "NK",
  core$orig_celltype_lvl_3
)

core$orig_celltype_lvl_3 <- ifelse(
  core$orig_celltype_lvl_3 == "RGS5+ EC",
  "Venous EC",
  core$orig_celltype_lvl_3
)

# Extract annotations
df <- FetchData(core, vars = c("orig_celltype_lvl_3", "goldstandard_celltype_lvl_3"))

dim(df)
# Filter out irrelevant cell types for the comparison
df <- df %>%
  filter(
    !grepl("^Other", orig_celltype_lvl_3),
    !goldstandard_celltype_lvl_3 %in% c("Fibro E", "Fibro D", "Merkel cell")
  )

dim(df)

# Clean up
df <- df %>%
  filter(!is.na(orig_celltype_lvl_3), !is.na(goldstandard_celltype_lvl_3))

dim(df)

# Get all unique cell types from both predicted and gold standard
celltypes <- union(unique(df$orig_celltype_lvl_3), unique(df$goldstandard_celltype_lvl_3))
celltypes

# Initialize result data frame
results <- data.frame(
  CellType = character(),
  TP = integer(),
  FP = integer(),
  FN = integer(),
  TN = integer(),
  Precision = numeric(),
  Recall = numeric(),
  F1 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each cell type (ct)
for (ct in celltypes) {
  # True Positive: predicted == ct and actual == ct
  TP <- sum(df$orig_celltype_lvl_3 == ct & df$goldstandard_celltype_lvl_3 == ct)

  # False Positive: predicted == ct but actual != ct
  FP <- sum(df$orig_celltype_lvl_3 == ct & df$goldstandard_celltype_lvl_3 != ct)

  # False Negative: predicted != ct but actual == ct
  FN <- sum(df$orig_celltype_lvl_3 != ct & df$goldstandard_celltype_lvl_3 == ct)

  # True Negative: predicted != ct and actual != ct
  TN <- sum(df$orig_celltype_lvl_3 != ct & df$goldstandard_celltype_lvl_3 != ct)

  # Precision, Recall, F1
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA,
               2 * precision * recall / (precision + recall))

  # Append to results
  results <- rbind(results, data.frame(
    CellType = ct,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    Precision = round(precision, 3),
    Recall = round(recall, 3),
    F1 = round(f1, 3)
  ))
}

# View results
print(results)


p1 <- ggplot(results, aes(x = F1, y = reorder(CellType, F1), fill = F1)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(100)) +
  theme_minimal() +
  geom_text(aes(label = round(F1, 2)), hjust = -0.1, size = 3) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p1


ggsave(
  filename = "./HSCA_annot_performance.pdf",
  plot = p1,
  device = "pdf",
  width = 7,
  height = 9
)
