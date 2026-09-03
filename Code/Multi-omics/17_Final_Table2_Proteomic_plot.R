rm(list=ls())
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

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

cat("Rows with zero/NA variance (removed):\n")
print(names(row_var)[row_var == 0 | is.na(row_var)])

cat("Columns with zero/NA variance (removed):\n")
print(names(col_var)[col_var == 0 | is.na(col_var)])

mat <- mat[!(row_var == 0 | is.na(row_var)), !(col_var == 0 | is.na(col_var)), drop = FALSE]
sig_mat <- sig_mat[rownames(mat), colnames(mat)]

cat("Final mat dim after filtering:", dim(mat), "\n")

## ------------------------------------------------------------
## Safe NA-tolerant distance function
## ------------------------------------------------------------
dist_na <- function(m) {
  cor_mat <- suppressWarnings(cor(t(m), use = "pairwise.complete.obs"))
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 1
  as.dist(1 - cor_mat)
}

row_clust <- hclust(dist_na(mat), method = "average")
col_clust <- hclust(dist_na(t(mat)), method = "average")

## ------------------------------------------------------------
## Color scale (symmetric around 0, based on actual data range)
## ------------------------------------------------------------
max_abs <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("navy", "white", "darkred"))

## ------------------------------------------------------------
## Build heatmap
## ------------------------------------------------------------
ht <- Heatmap(
  mat,
  name = "Estimate",
  col = col_fun,
  na_col = "grey80",
  cluster_rows = row_clust,
  cluster_columns = col_clust,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),      # smaller x-axis (column) label font
  column_names_rot = 90,                     # angle labels so they don't overlap
  column_title = "Protein",
  row_title = "Gene",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  row_title_gp = gpar(fontsize = 13, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
  },
  heatmap_legend_param = list(
    title = "Estimate",
    direction = "horizontal",
    legend_width = unit(4, "cm")
  )
)

cat("Heatmap object created successfully. Class:", class(ht), "\n")

## ------------------------------------------------------------
## Draw and save
## ------------------------------------------------------------
title_txt <- "Estimate (Gene x Protein) — clustered rows & columns\n* FDR < 0.05"

pdf(file.path(out_dir, "Heatmap_Gene_Protein_Estimate.pdf"), width = 35, height = 8)
draw(ht, heatmap_legend_side = "top",
     column_title = title_txt,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

jpeg(file.path(out_dir, "Heatmap_Gene_Protein_Estimate.jpeg"), width = 35, height = 8, units = "in", res = 300)
draw(ht, heatmap_legend_side = "top",
     column_title = title_txt,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

# draw(ht, heatmap_legend_side = "top",
#      column_title = title_txt,
#      column_title_gp = gpar(fontsize = 14, fontface = "bold"))

cat("Done. Heatmap saved to:", out_dir, "\n")

