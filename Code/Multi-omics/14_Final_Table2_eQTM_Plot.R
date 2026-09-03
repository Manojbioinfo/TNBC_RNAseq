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
## Build combined_eQTM
## ------------------------------------------------------------
gene <- read.csv("Results_eQTM/Tables/gene_regions_250kb.csv")

in_dir <- file.path(base_dir, "Results_eQTM_Mval_0kb/", "Results_per_CpG")
files  <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
cat("eQTM files:", length(files), "\n")

list_eQTM     <- lapply(files, readRDS)
combined_eQTM <- do.call(rbind, list_eQTM)

combined_eQTM$NearbyGene <- combined_eQTM$Gene
combined_eQTM$DEGs       <- combined_eQTM$TargetGene

combined_eQTM <- combined_eQTM %>%
  filter(NearbyGene %in% gene$NAME) %>%
  select(DEGs, CpG, NearbyGene, t, pvalue)

combined_eQTM$FDR <- p.adjust(combined_eQTM$pvalue, method = "fdr")
combined_eQTM$NearbyGene <- factor(combined_eQTM$NearbyGene, levels = gene$NAME)
combined_eQTM <- combined_eQTM %>% arrange(NearbyGene)

cat("Total rows:", nrow(combined_eQTM), "\n")
cat("Unique CpGs:", length(unique(combined_eQTM$CpG)), "\n")

saveRDS(combined_eQTM, file.path(out_dir, "Combined_eQTM_Results_ALL.rds"))
combined_eQTM_sig <- combined_eQTM %>% filter(FDR < 0.05)
saveRDS(combined_eQTM_sig, file.path(out_dir, "Combined_eQTM_Results_sig.rds"))
cat("Significant (FDR<0.05):", nrow(combined_eQTM_sig), "\n")

## ------------------------------------------------------------
## Base data prep
## ------------------------------------------------------------
base_df <- combined_eQTM %>%
  select(CpG, NearbyGene, t, FDR) %>%
  distinct() %>%
  mutate(
    CpG_Gene = paste(CpG, NearbyGene, sep = "_"),
    sig = ifelse(FDR < 0.05, "*", "")
  )

## ==============================================================
## MATRIX 1: CpG (rows) x Gene (cols) — DENSE, full t-stat matrix
## ==============================================================
mat1 <- base_df %>%
  select(CpG, NearbyGene, t) %>%
  pivot_wider(names_from = NearbyGene, values_from = t) %>%
  as.data.frame()
rownames(mat1) <- mat1$CpG
mat1$CpG <- NULL
mat1 <- as.matrix(mat1)

sig1 <- base_df %>%
  select(CpG, NearbyGene, sig) %>%
  pivot_wider(names_from = NearbyGene, values_from = sig) %>%
  as.data.frame()
rownames(sig1) <- sig1$CpG
sig1$CpG <- NULL
sig1 <- sig1[rownames(mat1), colnames(mat1)]
sig1[is.na(sig1)] <- ""
sig1 <- as.matrix(sig1)

cat("mat1 dim:", dim(mat1), " | NAs:", sum(is.na(mat1)), "\n")

## ------------------------------------------------------------
## Color scale
## ------------------------------------------------------------
max_abs <- max(abs(mat1), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("navy", "white", "darkred"))

## ------------------------------------------------------------
## Transpose BOTH mat1 and sig1 together so indices stay aligned
## New orientation: rows = Gene, columns = CpG
## ------------------------------------------------------------
mat1 <- t(mat1)
sig1 <- t(sig1)

## ==============================================================
## PLOT 1: Gene x CpG — ComplexHeatmap, top legend, axis titles
## ==============================================================
ht1 <- Heatmap(
  mat1,
  name = "t-statistic",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_title = "CpG",
  row_title = "Gene",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  row_title_gp = gpar(fontsize = 13, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig1[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
  },
  heatmap_legend_param = list(
    title = "t-stat",
    direction = "horizontal",
    legend_width = unit(4, "cm")
  )
)

pdf(file.path(out_dir, "eQTM_heatmap_CpG_only2.pdf"), width = 6, height = 5)
draw(ht1,
     heatmap_legend_side = "top",
     column_title = "eQTM t-statistics (Gene x CpG)\n* FDR < 0.05",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

jpeg(file.path(out_dir, "eQTM_heatmap_CpG_only.jpeg"), width = 6, height = 5, units = "in", res = 300)
draw(ht1,
     heatmap_legend_side = "top",
     column_title = "eQTM t-statistics (Gene x CpG)\n* FDR < 0.05",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

## Draw once in current session too
draw(ht1,
     heatmap_legend_side = "top",
     column_title = "eQTM t-statistics (Gene x CpG)\n* FDR < 0.05",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))

cat("Done. Plot saved to:", out_dir, "\n")

## ------------------------------------------------------------
## Which gene(s) is each CpG associated with in the original data?
## ------------------------------------------------------------
cpg_gene_map <- combined_eQTM %>%
  select(CpG, NearbyGene) %>%
  distinct() %>%
  arrange(CpG)

cpg_gene_summary <- combined_eQTM %>%
  select(CpG, NearbyGene) %>%
  distinct() %>%
  group_by(CpG) %>%
  summarise(Genes = paste(NearbyGene, collapse = ", "), .groups = "drop")

write.csv(cpg_gene_map, file.path(out_dir, "CpG_Gene_Map.csv"), row.names = FALSE)
write.csv(cpg_gene_summary, file.path(out_dir, "CpG_Gene_Summary.csv"), row.names = FALSE)

print(cpg_gene_map, n = 50)
print(cpg_gene_summary, n = 30)
table(cpg_gene_summary$Genes)