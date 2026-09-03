rm(list=ls())


## ------------------------------------------------------------
## Libraries
## ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

out_dir <- "Result_Combined/"  # change to your output directory
data =readRDS("Result_Combined/Combined_Blood_Results_clean.rds")
## ------------------------------------------------------------
## Assume your data is in a data frame called `data`
## with columns: Gene, Celltype, Estimate, StdError, z-score, pval, FDR
## ------------------------------------------------------------

## Only mark FDR < 0.05 as "*"
data <- data %>%
  mutate(sig = case_when(
    FDR < 0.05 ~ "*",
    TRUE ~ ""
  ))

## ------------------------------------------------------------
## Build matrix of Estimates: Gene (rows) x Celltype (columns)
## ------------------------------------------------------------
mat <- data %>%
  select(Gene, Celltype, Estimate) %>%
  pivot_wider(names_from = Celltype, values_from = Estimate) %>%
  as.data.frame()

rownames(mat) <- mat$Gene
mat$Gene <- NULL
mat <- as.matrix(mat)

## ------------------------------------------------------------
## Build matching matrix of significance stars
## ------------------------------------------------------------
sig_mat <- data %>%
  select(Gene, Celltype, sig) %>%
  pivot_wider(names_from = Celltype, values_from = sig) %>%
  as.data.frame()

rownames(sig_mat) <- sig_mat$Gene
sig_mat$Gene <- NULL
sig_mat <- sig_mat[rownames(mat), colnames(mat)]
sig_mat[is.na(sig_mat)] <- ""
sig_mat <- as.matrix(sig_mat)

## ------------------------------------------------------------
## Color scale (based on Estimate values)
## ------------------------------------------------------------
max_abs <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("navy", "white", "darkred"))

## ------------------------------------------------------------
## Plot heatmap with legend on top
## ------------------------------------------------------------




ht <- Heatmap(
  mat,
  name = "Estimate",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_title = "Gene x Celltype Association (Estimate)\n* FDR < 0.05",
  heatmap_legend_param = list(
    title = "Estimate",
    direction = "horizontal",
    legend_width = unit(5, "cm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
  }
)

## ------------------------------------------------------------
## Save plots (Blood)
## ------------------------------------------------------------
pdf(file.path(out_dir, "Blood_GeneCelltype_heatmap.pdf"), width =12, height = 5)
draw(ht, heatmap_legend_side = "top")
dev.off()

jpeg(file.path(out_dir, "Blood_GeneCelltype_heatmap.jpeg"), width = 12, height = 5, units = "in", res = 300)
draw(ht, heatmap_legend_side = "top")
dev.off()

cat("Done. Blood heatmap saved to:", out_dir, "\n")
