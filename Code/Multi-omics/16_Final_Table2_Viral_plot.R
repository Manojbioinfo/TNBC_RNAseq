rm(list=ls())
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

base_dir <- "Z:/LAPTOP/Newmanuscript/Paper/2_Working/W2_Gopi_TNBC/Multiomic_Aug26_Complete_Data"
out_dir  <- file.path(base_dir, "Result_Combined")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## Load data
## ------------------------------------------------------------
data <- readRDS(file.path(out_dir, "Combined_Viral_Results_clean.rds"))

data <- data %>%
  mutate(sig = case_when(
    is.na(FDR) ~ "",
    FDR < 0.05 ~ "*",
    TRUE ~ ""
  ))

range(data$Estimate)

data=dplyr::filter(data,Virus!="HTLV")
data$Virus=gsub("bk_polyomavirus","BKP",data$Virus)
data$Virus=gsub("saimiriine_herpesvirus","SHV",data$Virus)


## ------------------------------------------------------------
## Build matrix: Gene (rows) x Virus (columns) — Estimate values
## ------------------------------------------------------------
mat_raw <- data %>%
  select(Gene, Virus, Estimate) %>%
  pivot_wider(names_from = Virus, values_from = Estimate) %>%
  as.data.frame()

rownames(mat_raw) <- mat_raw$Gene
mat_raw$Gene <- NULL
mat_raw <- as.matrix(mat_raw)

sig_mat <- data %>%
  select(Gene, Virus, sig) %>%
  pivot_wider(names_from = Virus, values_from = sig) %>%
  as.data.frame()

rownames(sig_mat) <- sig_mat$Gene
sig_mat$Gene <- NULL
sig_mat <- sig_mat[rownames(mat_raw), colnames(mat_raw)]
sig_mat[is.na(sig_mat)] <- ""
sig_mat <- as.matrix(sig_mat)

## ------------------------------------------------------------
## Remove all-NA rows/columns
## ------------------------------------------------------------
mat_raw <- mat_raw[rowSums(!is.na(mat_raw)) > 1, colSums(!is.na(mat_raw)) > 1]
sig_mat <- sig_mat[rownames(mat_raw), colnames(mat_raw)]

## ------------------------------------------------------------
## Sign-preserving log-scaling so tiny AND large values are both visible
## Formula: sign(x) * log10(1 + |x| / eps)  then rescale to [-0.1, 0.1]
## ------------------------------------------------------------
eps <- 1e-9  # floor scale, values below this look like ~0
mat_log <- sign(mat_raw) * log10(1 + abs(mat_raw) / eps)

max_abs_log <- max(abs(mat_log), na.rm = TRUE)
mat_scaled <- mat_log / max_abs_log * 0.1   # rescale into [-0.1, 0.1]

## ------------------------------------------------------------
## Safe NA-tolerant distance function
## ------------------------------------------------------------
dist_na <- function(m) {
  cor_mat <- suppressWarnings(cor(t(m), use = "pairwise.complete.obs"))
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 1
  as.dist(1 - cor_mat)
}

## ==============================================================
## Function to build + draw + save a heatmap
## ==============================================================
make_heatmap <- function(mat, sig_mat, label, legend_title, filename_base) {
  
  col_fun <- colorRamp2(c(-0.1, 0, 0.1), c("navy", "white", "darkred"))
  
  row_clust <- hclust(dist_na(mat), method = "average")
  col_clust <- hclust(dist_na(t(mat)), method = "average")
  
  ht <- Heatmap(
    mat,
    name = legend_title,
    col = col_fun,
    na_col = "grey80",
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title = "Virus",
    row_title = "Gene",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
    },
    heatmap_legend_param = list(
      title = legend_title,
      direction = "horizontal",
      legend_width = unit(4, "cm"),
      at = c(-0.1, -0.05, 0, 0.05, 0.1),
      labels = c("- large", "- small", "0", "+ small", "+ large")
    )
  )
  
  title_txt <- paste0(label, " (Gene x Virus) — clustered rows & columns\n* FDR < 0.05\nSign-preserving log-scaled to [-0.1, 0.1]")
  
  pdf(file.path(out_dir, paste0(filename_base, ".pdf")), width = 7, height = 5)
  draw(ht, heatmap_legend_side = "top",
       column_title = title_txt,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
  
  jpeg(file.path(out_dir, paste0(filename_base, ".jpeg")), width = 7, height = 5, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "top",
       column_title = title_txt,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
  
  draw(ht, heatmap_legend_side = "top",
       column_title = title_txt,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  
  return(ht)
}

## ==============================================================
## Build heatmap using log-scaled matrix
## ==============================================================
ht_scaled <- make_heatmap(
  mat = mat_scaled,
  sig_mat = sig_mat,
  label = "Estimate (log-scaled)",
  legend_title = "Scaled\nEstimate",
  filename_base = "Viral_heatmap_Gene_Virus_logscaled"
)

cat("Done. Log-scaled heatmap saved to:", out_dir, "\n")

