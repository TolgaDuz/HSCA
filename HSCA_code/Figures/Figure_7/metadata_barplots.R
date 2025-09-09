# setting working dir
setwd("~/HSCA_code/Figures/Figure_7/")

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
library(scales)

# Load the final, filtered, and annotated atlas object
atlas <-  readRDS("~/HSCA_data/extended/processed_seurat_objects/HSCA_extended.rds")

DimPlot(atlas, label = T)

# ========================================================
# Metadata inspection
# ========================================================

atlas$sample %>% table()
atlas$sample %>% unique() %>% length()
atlas$sex %>% table() # zeigen
atlas$subject_ID %>% table()
atlas$subject_ID %>% unique() %>% length()
atlas$age_years %>% table()
atlas$age_range %>% table()
atlas$ethnicity %>% table()
atlas$anatomical_region_level1 %>% table()
atlas$anatomical_region_level2 %>% table()
atlas$anatomical_region_level3 %>% table()
atlas$sequencing_platform %>% table()
atlas$reference_genome_compact %>% table()
atlas$cell_ranger_version_compact %>% table()
atlas$single_cell_platform %>% table()
atlas$single_cell_platform_specific %>% table()
atlas$tissue_sampling_type %>% table()
atlas$Accession_source %>% table()
atlas$Dataset %>% table()
atlas$Dataset %>% unique() %>% length()
atlas$Condition %>% table()
atlas$Core %>% table()

# ========================================================
# Define variables of interest for barplots
# ========================================================

cols_to_plot <- c(
  "sex",
  "ethnicity",
  "anatomical_region_level1",
  "anatomical_region_level3",
  "sequencing_platform",
  "single_cell_platform_specific",
  "tissue_sampling_type"
)

# Collect all possible levels for each metadata field
# - Exclude true NA values
# - Ensure "nan" string is always present as an explicit level
all_levels_list <- lapply(cols_to_plot, function(col) {
  vals <- unique(atlas@meta.data[[col]])
  vals <- vals[!is.na(vals)]
  vals <- as.character(vals)
  vals <- unique(c(vals, "nan"))
  sort(vals)
})

names(all_levels_list) <- cols_to_plot


# ========================================================
# Function to plot metadata distributions
# ========================================================

plot_metadata_barplot <- function(df, col_name, all_levels = NULL) {
  # Count values, excluding NA but keeping "nan"
  df_plot <- df %>%
    mutate(value = .data[[col_name]]) %>%
    filter(!is.na(value)) %>%
    mutate(value = ifelse(value == "nan", "nan", as.character(value))) %>%
    count(value)

  # Order categories by frequency (nan always last)
  ordered_levels <- df_plot %>%
    filter(value != "nan") %>%
    arrange(desc(n)) %>%
    pull(value)
  final_levels <- ordered_levels
  if ("nan" %in% df_plot$value) {
    final_levels <- c(final_levels, "nan")
  }

  # If full set of levels provided, add missing ones with count = 0
  if (!is.null(all_levels)) {
    missing_levels <- setdiff(all_levels, final_levels)
    if (length(missing_levels) > 0) {
      df_plot <- bind_rows(df_plot, tibble(value = missing_levels, n = 0))
      final_levels <- all_levels
    }
  }

  # Convert to factor with defined order
  df_plot$value <- factor(df_plot$value, levels = final_levels)

  # Create barplot
  p <- ggplot(df_plot, aes(x = value, y = n)) +
    geom_col(fill = "#1e347b", width = 0.7) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = comma
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(x = col_name, y = "Count", title = paste("Counts per", col_name))

  return(p)
}

# ========================================================
# Generate and save barplots
# ========================================================


bar_target_width_cm <- 0.7

for (col in cols_to_plot) {
  p <- plot_metadata_barplot(atlas@meta.data, col, all_levels_list[[col]])
  print(p)

  n_levels <- length(levels(p$data$value))
  plot_width_cm <- n_levels * (bar_target_width_cm / 0.7)

  plot_height_cm <- if (grepl("anatomical_region_level1", col)) 9 else 6

  filename <- paste0("./barplot_", col, ".pdf")
  ggsave(
    filename,
    plot = p,
    device = "pdf",
    width = plot_width_cm,
    height = plot_height_cm
  )
}

# ========================================================
# Special case: Age range barplot
# ========================================================

# Define natural order of age groups
age_levels <- c("0-1", "13-25", "26-40", "41-60", "40-80", "61-80", "nan")

plot_age_range <- function(df, age_levels) {
  df_plot <- df %>%
    mutate(age_range = .data[["age_range"]]) %>%
    filter(!is.na(age_range)) %>%
    mutate(age_range = ifelse(age_range == "nan", "nan", as.character(age_range))) %>%
    count(age_range) %>%
    complete(age_range = age_levels, fill = list(n = 0))

  df_plot$age_range <- factor(df_plot$age_range, levels = age_levels)

  df_plot$fill_color <- ifelse(df_plot$age_range == "nan", "gray", "#1e347b")

  p <- ggplot(df_plot, aes(x = age_range, y = n, fill = fill_color)) +
    geom_col(width = 0.7) +
    scale_fill_identity() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = comma
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(x = "Age Range", y = "Count", title = "Counts per Age Range")

  return(p)
}

# Generate and save age range plot
p_age <- plot_age_range(atlas@meta.data, age_levels)
print(p_age)

bar_target_width_cm <- 0.7
n_levels <- length(levels(p_age$data$age_range))
plot_width_cm <- n_levels * (bar_target_width_cm / 0.7)

filename <- "./barplot_age_range.pdf"
ggsave(filename, plot = p_age, device = "pdf", width = plot_width_cm, height = 7)
