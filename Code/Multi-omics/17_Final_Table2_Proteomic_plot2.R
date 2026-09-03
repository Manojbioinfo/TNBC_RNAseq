rm(list=ls())
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
library(ggrepel)

## ------------------------------------------------------------
## Paths
## ------------------------------------------------------------
base_dir <- "Z:/LAPTOP/Newmanuscript/Paper/2_Working/W2_Gopi_TNBC/Multiomic_Aug26_Complete_Data"
out_dir  <- file.path(base_dir, "Result_Combined")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Load data
## ------------------------------------------------------------
data <- readRDS(file.path(out_dir, "Combined_Proteomic_Results_clean.rds"))

cat("Data loaded. Dim:", dim(data), "\n")
cat("Columns:", paste(colnames(data), collapse = ", "), "\n")

data <- data %>%
  mutate(sig = case_when(
    is.na(FDR) ~ "",
    FDR < 0.05 ~ "*",
    TRUE ~ ""
  ))

## ------------------------------------------------------------
## Build matrix: Gene (rows) x Protein (columns) — Estimate values
## ------------------------------------------------------------
mat <- data %>%
  select(Gene, Protein, Estimate) %>%
  pivot_wider(names_from = Protein, values_from = Estimate) %>%
  as.data.frame()

rownames(mat) <- mat$Gene
mat$Gene <- NULL
mat <- as.matrix(mat)

## ------------------------------------------------------------
## Build matching matrix of significance stars
## ------------------------------------------------------------
sig_mat <- data %>%
  select(Gene, Protein, sig) %>%
  pivot_wider(names_from = Protein, values_from = sig) %>%
  as.data.frame()

rownames(sig_mat) <- sig_mat$Gene
sig_mat$Gene <- NULL
sig_mat <- sig_mat[rownames(mat), colnames(mat)]
sig_mat[is.na(sig_mat)] <- ""
sig_mat <- as.matrix(sig_mat)

cat("Initial mat dim:", dim(mat), " | NAs:", sum(is.na(mat)), "\n")

## ------------------------------------------------------------
## Remove all-NA rows/columns
## ------------------------------------------------------------
mat <- mat[rowSums(!is.na(mat)) > 1, colSums(!is.na(mat)) > 1, drop = FALSE]
sig_mat <- sig_mat[rownames(mat), colnames(mat)]

## ------------------------------------------------------------
## Remove zero-variance / constant rows & columns
## ------------------------------------------------------------
row_var <- apply(mat, 1, var, na.rm = TRUE)
col_var <- apply(mat, 2, var, na.rm = TRUE)

keep_rows <- !is.na(row_var) & row_var > 0
keep_cols <- !is.na(col_var) & col_var > 0

cat("Rows with zero/NA variance (removed):", sum(!keep_rows), "\n")
cat("Cols with zero/NA variance (removed):", sum(!keep_cols), "\n")

mat <- mat[keep_rows, keep_cols, drop = FALSE]
sig_mat <- sig_mat[keep_rows, keep_cols, drop = FALSE]

## ------------------------------------------------------------
## Cluster settings (guard against too-few rows/cols to cluster)
## ------------------------------------------------------------
row_clust <- nrow(mat) > 2
col_clust <- ncol(mat) > 2

n_rows <- nrow(mat)
n_cols <- ncol(mat)
cat("Final matrix used for heatmap: rows =", n_rows, " cols =", n_cols, "\n")

## ------------------------------------------------------------
## Transpose: CpGs/Proteins on Y axis, Genes on X axis
## ------------------------------------------------------------
mat_t     <- t(mat)   # rows = Protein/CpG, cols = Gene
sig_mat_t <- t(sig_mat)

## ------------------------------------------------------------
## Color scale
## ------------------------------------------------------------
max_abs <- max(abs(mat_t), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("navy", "white", "darkred"))

## ------------------------------------------------------------
## FIXED FIGURE SIZE
## ------------------------------------------------------------
fig_width  <- 7
fig_height <- 12

## ------------------------------------------------------------
## Dynamic font size for row labels (Protein/CpG names) based on count
## Smaller when many rows, but with a sensible floor
## ------------------------------------------------------------
row_fontsize <- max(4, min(9, 300 / nrow(mat_t)))
col_fontsize <- max(5, min(10, 200 / ncol(mat_t)))

cat("Row fontsize used:", row_fontsize, "\n")
cat("Col fontsize used:", col_fontsize, "\n")

## ------------------------------------------------------------
## Build heatmap - fixed width/height, non-overlapping row labels via
## ggrepel-style approach: we manually create repelled label positions
## using ggplot2 + ggrepel and overlay as annotation, OR use ComplexHeatmap
## row label repel via 'anno_text' along with 'link' annotation.
## ------------------------------------------------------------

## Step 1: create row index and repelled y-positions using ggrepel logic
row_labels <- rownames(mat_t)
n <- length(row_labels)

label_df <- data.frame(
  y = seq_len(n),
  label = row_labels
)

## Use ggrepel's internal repel algorithm via a dummy ggplot to compute
## non-overlapping y-positions, then map back to ComplexHeatmap row annotation
dummy_plot <- ggplot(label_df, aes(x = 1, y = y, label = label)) +
  geom_point(alpha = 0) +
  geom_text_repel(
    size = row_fontsize / .pt,
    direction = "y",
    max.overlaps = Inf,
    force = 2,
    box.padding = 0.1,
    segment.size = 0.2,
    segment.color = "grey50",
    hjust = 0
  ) +
  theme_void()

## ------------------------------------------------------------
## Use ComplexHeatmap's built-in row annotation with repel via
## anno_mark() -- this is the correct, robust way to get non-overlapping
## row labels with leader lines, fixed heatmap body size.
## ------------------------------------------------------------
row_anno <- rowAnnotation(
  labels = anno_mark(
    at = seq_len(n_rows) |> (\(x) which(rownames(mat_t) %in% rownames(mat_t)))(),
    labels = rownames(mat_t),
    labels_gp = gpar(fontsize = row_fontsize),
    padding = unit(1, "mm"),
    link_width = unit(4, "mm"),
    extend = unit(2, "mm")
  )
)

## ------------------------------------------------------------
## Build heatmap object (fixed body size independent of row/col count,
## since overall figure width/height is fixed at 7 x 16 inches)
## ------------------------------------------------------------
ht2 <- Heatmap(
  mat_t,
  name = "Estimate",
  col = col_fun,
  na_col = "grey80",
  
  cluster_rows = row_clust,
  cluster_columns = col_clust,
  
  show_row_names = FALSE,     # turned off; using anno_mark repel labels instead
  right_annotation = row_anno,
  
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = col_fontsize, fontface = "plain"),
  column_names_rot = 45,
  column_names_side = "bottom",
  column_names_max_height = unit(3, "cm"),
  
  column_title = "Gene",
  row_title = "Protein / CpG",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  row_title_gp = gpar(fontsize = 13, fontface = "bold"),
  
  width  = unit(fig_width - 3, "in"),   # leave room for repelled labels + legend
  height = unit(fig_height - 2, "in"),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat_t[i, j], x, y, gp = gpar(fontsize = 7, col = "black"))
  },
  
  heatmap_legend_param = list(
    title = "Estimate",
    direction = "horizontal",
    legend_width = unit(4, "cm")
  )
)

cat("Heatmap object 'ht2' built successfully. Class:", class(ht2), "\n")

## ------------------------------------------------------------
## Save at FIXED size: width = 7 in, height = 16 in
## ------------------------------------------------------------
title_txt <- "Estimate (Protein/CpG x Gene) — clustered\n* FDR < 0.05"

pdf(file.path(out_dir, "Heatmap_Protein_Gene_Estimate_Tall.pdf"),
    width = fig_width, height = fig_height)
draw(ht2, heatmap_legend_side = "top",
     column_title = title_txt,
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

jpeg(file.path(out_dir, "Heatmap_Protein_Gene_Estimate_Tall.jpeg"),
     width = fig_width, height = fig_height, units = "in", res = 300)
draw(ht2, heatmap_legend_side = "top",
     column_title = title_txt,
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

cat("Heatmap saved to:", out_dir, "\n")
cat("Fixed dimensions used: width =", fig_width, "| height =", fig_height, "\n")