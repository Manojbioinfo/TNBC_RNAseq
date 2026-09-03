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

data <- data %>%
  mutate(sig = case_when(
    is.na(FDR) ~ "",
    FDR < 0.05 ~ "*",
    TRUE ~ ""
  ))

## ------------------------------------------------------------
## Build matrix: Protein (rows) x Gene (columns) — Estimate
## ------------------------------------------------------------
mat_t <- data %>%
  select(Gene, Protein, Estimate) %>%
  pivot_wider(names_from = Gene, values_from = Estimate) %>%
  as.data.frame()
rownames(mat_t) <- mat_t$Protein
mat_t$Protein <- NULL
mat_t <- as.matrix(mat_t)

sig_mat <- data %>%
  select(Gene, Protein, sig) %>%
  pivot_wider(names_from = Gene, values_from = sig) %>%
  as.data.frame()
rownames(sig_mat) <- sig_mat$Protein
sig_mat$Protein <- NULL
sig_mat <- as.matrix(sig_mat)
sig_mat[is.na(sig_mat)] <- ""

## ------------------------------------------------------------
## >>> YOUR 6 TARGET GENES <<<
## ------------------------------------------------------------
six_genes_input <- c("NMU", "NEK2", "HORMAD1", "PRAME", "KIF18B", "MMP1")
six_genes <- intersect(six_genes_input, colnames(mat_t))

if (length(six_genes) == 0) {
  stop("None of the requested genes were found in colnames(mat_t). Check spelling/case.")
}

sub_mat <- mat_t[, six_genes, drop = FALSE]
sig_sub <- sig_mat[rownames(sub_mat), six_genes, drop = FALSE]

keep_rows <- rowSums(!is.na(sub_mat)) > 0
sub_mat <- sub_mat[keep_rows, , drop = FALSE]
sig_sub <- sig_sub[keep_rows, , drop = FALSE]

## ------------------------------------------------------------
## Classify each PROTEIN (row) by direction across the 6 genes
## ------------------------------------------------------------
row_direction <- apply(sub_mat, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("Mixed")
  if (all(x > 0)) return("Positive")
  if (all(x < 0)) return("Negative")
  return("Mixed")
})

proteins_neg   <- names(row_direction)[row_direction == "Negative"]
proteins_pos   <- names(row_direction)[row_direction == "Positive"]
proteins_mixed <- names(row_direction)[row_direction == "Mixed"]

## ------------------------------------------------------------
## Transpose so Gene = rows, Protein = columns (for plotting)
## ------------------------------------------------------------
mat_neg_t   <- t(sub_mat[proteins_neg, , drop = FALSE])
mat_pos_t   <- t(sub_mat[proteins_pos, , drop = FALSE])
mat_mixed_t <- t(sub_mat[proteins_mixed, , drop = FALSE])

sig_neg_t   <- t(sig_sub[proteins_neg, , drop = FALSE])
sig_pos_t   <- t(sig_sub[proteins_pos, , drop = FALSE])
sig_mixed_t <- t(sig_sub[proteins_mixed, , drop = FALSE])

## ------------------------------------------------------------
## Shared color scale across all 3 panels
## ------------------------------------------------------------
all_vals <- c(mat_neg_t, mat_pos_t, mat_mixed_t)
all_vals <- all_vals[!is.na(all_vals)]
max_abs  <- max(abs(all_vals))

col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("navy", "white", "darkred"))

## ------------------------------------------------------------
## Cell function factory (adds significance asterisks)
## ------------------------------------------------------------
make_cell_fun <- function(sig_mat_t) {
  function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(fill = fill, col = "grey80", lwd = 0.5))
    lab <- sig_mat_t[i, j]
    if (!is.na(lab) && lab != "") {
      grid.text(lab, x = x, y = y, gp = gpar(fontsize = 9, col = "black"))
    }
  }
}

## ------------------------------------------------------------
## Build each Heatmap object independently (own clustering)
## ------------------------------------------------------------
build_ht <- function(mat_t, sig_t, title) {
  if (nrow(mat_t) == 0 || ncol(mat_t) == 0) return(NULL)
  Heatmap(
    mat_t,
    name = "Estimate",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_dend = TRUE,
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 13, fontface = "bold"),
    row_title = paste0(title, " (n=", ncol(mat_t), ")"),
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    cell_fun = make_cell_fun(sig_t),
    show_heatmap_legend = FALSE,
    border = TRUE
  )
}

ht_neg   <- build_ht(mat_neg_t,   sig_neg_t,   "Negative")
ht_pos   <- build_ht(mat_pos_t,   sig_pos_t,   "Positive")
ht_mixed <- build_ht(mat_mixed_t, sig_mixed_t, "Mixed")

## ------------------------------------------------------------
## Shared legend
## ------------------------------------------------------------
shared_legend <- Legend(
  col_fun = col_fun,
  title = "Estimate",
  direction = "horizontal",
  legend_width = unit(4, "cm"),
  title_position = "topcenter"
)

## ------------------------------------------------------------
## Draw 3-row layout with EQUAL panel heights
## ------------------------------------------------------------
draw_three_row_layout <- function(title_txt) {
  grid.newpage()
  
  pushViewport(viewport(layout = grid.layout(
    nrow = 5, ncol = 1,
    heights = unit(c(1.2, 1, 1, 1, 1),
                   c("lines", "lines", "null", "null", "null"))
  )))
  
  # Title
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(title_txt, gp = gpar(fontsize = 13, fontface = "bold"))
  popViewport()
  
  # Shared legend
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(shared_legend, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
  popViewport()
  
  # Row 3: Negative
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  if (!is.null(ht_neg)) draw(ht_neg, newpage = FALSE)
  popViewport()
  
  # Row 4: Positive
  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  if (!is.null(ht_pos)) draw(ht_pos, newpage = FALSE)
  popViewport()
  
  # Row 5: Mixed
  pushViewport(viewport(layout.pos.row = 5, layout.pos.col = 1))
  if (!is.null(ht_mixed)) draw(ht_mixed, newpage = FALSE)
  popViewport()
  
  popViewport()
}

## ------------------------------------------------------------
## Figure size
## ------------------------------------------------------------
fig_width  <- 16
fig_height <- 10
title_txt <- "Estimate — Split by Direction across 6 Genes (* FDR < 0.05)"

## ------------------------------------------------------------
## Save PDF
## ------------------------------------------------------------
pdf(file.path(out_dir, "Heatmap_SixGene_ThreePanels_3Row_Equal.pdf"),
    width = fig_width, height = fig_height)
draw_three_row_layout(title_txt)
dev.off()

## ------------------------------------------------------------
## Save JPEG
## ------------------------------------------------------------
jpeg(file.path(out_dir, "Heatmap_SixGene_ThreePanels_3Row_Equal.jpeg"),
     width = fig_width, height = fig_height, units = "in", res = 300)
draw_three_row_layout(title_txt)
dev.off()

cat("\nSaved 3-row EQUAL-height layout to:", out_dir, "\n")
