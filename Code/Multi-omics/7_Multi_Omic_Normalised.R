rm(list = ls())

# Load required libraries
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)

dir.create("Normalised2", showWarnings = FALSE)

# =========================================================
# Helper function: plot density of all samples in a dataset
# =========================================================
plot_density <- function(data_long, title, filename_prefix) {
  p <- ggplot(data_long, aes(x = Value, group = Sample)) +
    geom_density(alpha = 0.05, linewidth = 0.2, color = "steelblue") +
    theme_minimal() +
    labs(title = title, x = "Value", y = "Density")
  
  ggsave(paste0("Normalised2/", filename_prefix, ".png"), p, width = 10, height = 7, dpi = 300)
  ggsave(paste0("Normalised2/", filename_prefix, ".pdf"), p, width = 10, height = 7)
  
  return(p)
}

# =========================================================
# Helper function: full pipeline for one dataset
# Applies LOG transform and INVERSE-NORMAL (rank-based) 
# transform SEPARATELY on the RAW data (not chained)
# =========================================================
process_dataset <- function(file_path, dataset_name) {
  
  cat("\n=== Processing:", dataset_name, "===\n")
  
  # ---- Load raw data ----
  data <- readRDS(file_path)
  cat("Dimensions:", dim(data), "\n")
  
  mat <- as.matrix(data)
  
  # ---- Long format of RAW data ----
  data_long_raw <- as.data.frame(mat) %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Sample", values_to = "Value")
  
  # ---- 1. BEFORE normalization: plot raw density ----
  plot_density(data_long_raw,
               paste0(dataset_name, ": Raw (before normalization)"),
               paste0(dataset_name, "_BeforeNorm_density"))
  
  # =========================================================
  # OPTION A: LOG2 transform (applied to RAW data directly)
  # =========================================================
  min_val <- min(mat, na.rm = TRUE)
  offset <- ifelse(min_val <= 0, abs(min_val) + 1e-3, 0)
  mat_log <- log2(mat + offset)
  
  data_long_log <- as.data.frame(mat_log) %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Sample", values_to = "Value")
  
  plot_density(data_long_log,
               paste0(dataset_name, ": Log2 transformed (from raw)"),
               paste0(dataset_name, "_Log2_density"))
  
  saveRDS(mat_log, file = paste0("Normalised2/", dataset_name, "_log2.rds"))
  
  # =========================================================
  # OPTION B: INVERSE NORMAL (rank-based quantile) transform
  # Applied to RAW data directly, per column (per sample)
  # =========================================================
  mat_invnorm <- apply(mat, 2, function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  })
  rownames(mat_invnorm) <- rownames(mat)
  colnames(mat_invnorm) <- colnames(mat)
  
  data_long_invnorm <- as.data.frame(mat_invnorm) %>%
    rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Sample", values_to = "Value")
  
  plot_density(data_long_invnorm,
               paste0(dataset_name, ": Inverse-normal transformed (from raw)"),
               paste0(dataset_name, "_InverseNormal_density"))
  
  saveRDS(mat_invnorm, file = paste0("Normalised2/", dataset_name, "_inverse_normal.rds"))
  
  cat("Saved log2 and inverse-normal transformed data for:", dataset_name, "\n")
  
  return(list(log2 = mat_log, invnorm = mat_invnorm))
}

# =========================================================
# Run pipeline for all three datasets
# =========================================================
files <- list(
  methylation = "Normalised/methylation_aligned.rds",
  blood       = "Normalised/blood_aligned.rds",
  proteomic   = "Normalised/proteomic_aligned.rds"
)

results <- list()

for (name in names(files)) {
  results[[name]] <- process_dataset(files[[name]], name)
}

cat("\n=== All datasets processed and saved in Normalised2/ ===\n")
list.files("Normalised2/")