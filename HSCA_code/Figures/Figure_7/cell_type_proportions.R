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

setwd("~/HSCA_code/Figures/Figure_7")

# Load the HSCA extended atlas
extended_atlas <- readRDS("~/HSCA_data/extended/processed_seurat_objects/final_atlas/HSCA_extended.rds")

# Inspect
extended_atlas
DimPlot(extended_atlas, group.by = "inherited_celltype_lvl_1_extended", label = TRUE)
extended_atlas$inherited_celltype_lvl_4_extended %>% table()

# Extract metadata for further analysis
df <- extended_atlas@meta.data

df

df$Core %>% table()

# Create a new category: 'core' vs 'extended_only'
df <- df %>%
  mutate(dataset_group = ifelse(Core == "Yes", "core", "extended_only"))

df$dataset_group %>% table()

# ===============================================
# Percentage of cells per cell type
# ===============================================

counts <- df %>%
  group_by(inherited_celltype_lvl_4_extended, dataset_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(inherited_celltype_lvl_4_extended) %>%
  mutate(percent = n / sum(n) * 100) # Compute percentage within each cell type

# Set factor levels for cell types in descending order of core proportion
core_order <- counts %>%
  filter(dataset_group == "core") %>%
  arrange(desc(percent)) %>%
  pull(inherited_celltype_lvl_4_extended)

counts$inherited_celltype_lvl_4_extended <-
  factor(counts$inherited_celltype_lvl_4_extended, levels = core_order)

# Set factor levels for stacking: Core at the bottom, Extended-only on top
counts$dataset_group <- factor(counts$dataset_group, levels = c("core", "extended_only"))

# Create percentage bar plot (Figure 7f)
p <- ggplot(counts, aes(x = inherited_celltype_lvl_4_extended, y = percent, fill = dataset_group)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.6) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  labs(y = "Percentage of cells", x = "Cell type", fill = "Dataset") +
  scale_fill_manual(values = c("core" = "#00B4D8", "extended_only" = "#363636"))

# Display the plot
p

# Save Figure 7f
ggsave("./cell_type_proportions.pdf", plot = p, width = 14, height = 4)
