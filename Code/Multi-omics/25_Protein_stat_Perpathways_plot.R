## =============================================================================
## Multi-Omic Directionality Report: Proteomic, Blood, and Viral Associations
## For candidate genes: NMU, NEK2, HORMAD1, PRAME, KIF18B, MMP1
## =============================================================================

library(dplyr)
library(tidyr)
library(openxlsx)

data1=readRDS("Result_Combined/Combined_Proteomic_Results_clean.rds")
data2=read.csv("Result_Combined/RPPA/RPPA_Protein_Pathway_Map_Clean.csv")
head(data1)
head(data2)
#############


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

## ============================================================
## 1. Merge data1 (results) with data2 (pathway map)
## ============================================================

merged <- data1 %>%
  left_join(data2 %>% select(Protein, Pathway), by = "Protein")

cat("Unmatched proteins:", sum(is.na(merged$Pathway)), "\n")
merged$Pathway[is.na(merged$Pathway)] <- "Unassigned"


write.xlsx(merged,"Result_Combined/Combined_Proteomic_Results.xlsx")

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

## ============================================================
## 1. Summarise mean Estimate per Gene per Pathway
## ============================================================

summary_df <- merged %>%
  group_by(Gene, Pathway) %>%
  summarise(mean_estimate = mean(Estimate, na.rm = TRUE), .groups = "drop")

## ============================================================
## 2. Pivot to wide matrix: rows = Gene, columns = Pathway
## ============================================================

mat <- summary_df %>%
  pivot_wider(names_from = Pathway, values_from = mean_estimate) %>%
  as.data.frame()

rownames(mat) <- mat$Gene
mat$Gene <- NULL
mat <- as.matrix(mat)

## Remove rows that are entirely NA
mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

## Replace NA with 0 for display/clustering
mat[is.na(mat)] <- 0

## ============================================================
## 3. Build heatmap
## ============================================================

## Symmetric color scale centered at 0
max_abs <- max(abs(mat), na.rm = TRUE)

col_fun <- colorRamp2(
  c(-max_abs, 0, max_abs),
  c("navy", "white", "darkred")
)

ht <- Heatmap(
  mat,
  name = "Mean Estimate",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 90,
  heatmap_legend_param = list(title = "Mean Estimate")
)

pdf("Result_Combined/RPPA/Heatmap_Gene_by_Pathway_Estimate.pdf", width = 12, height = 6)
draw(ht, heatmap_legend_side = "top")
dev.off()

draw(ht, heatmap_legend_side = "right")



merged %>%
  filter(Pathway == "Apoptosis") %>%
  distinct(Protein)
