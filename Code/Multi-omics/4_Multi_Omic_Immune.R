rm(list = ls())
# Load required libraries
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DESeq2)
library(readr)
library(dplyr)
rnaseq <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")

colnames(rnaseq)



# =========================================================
# Complete pipeline: Align RNA-seq (DESeqDataSet) and 
# Immune cell signature data by sample ID, then save as .rds
# =========================================================

library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------
# 1. Load immune cell data
# ---------------------------------------------------------
immune <- readRDS("Z:/1PERSONALprojectcomplete/TCGAData/Immune_cells/TCGA_immunecell.rds")

cat("Immune dimensions:", dim(immune), "\n")

# Create a single combined row-identifier (Source + SetName)
immune <- as.data.frame(immune)
immune$feature_id <- paste(immune$Source, immune$SetName, sep = "_")

# Move feature_id to front, drop original Source/SetName (keep if you want them later)
immune <- immune[, c("feature_id", setdiff(colnames(immune), c("Source", "SetName", "feature_id")))]

# ---------------------------------------------------------
# 2. Extract raw count matrix from DESeqDataSet
# ---------------------------------------------------------
rnaseq_counts <- counts(rnaseq)   # genes x samples matrix

cat("RNAseq dimensions:", dim(rnaseq_counts), "\n")

# ---------------------------------------------------------
# 3. Truncate barcodes to 15-character sample IDs
# ---------------------------------------------------------
colnames(rnaseq_counts)   <- substr(colnames(rnaseq_counts), 1, 15)
colnames(immune)[-1]      <- substr(colnames(immune)[-1], 1, 15)

# ---------------------------------------------------------
# 4. Collapse duplicate RNA-seq columns (same sample ID) by SUM
# ---------------------------------------------------------
rnaseq_t <- as.data.frame(t(rnaseq_counts))
rnaseq_t$sample_id <- rownames(rnaseq_t)

rnaseq_collapsed <- aggregate(. ~ sample_id, data = rnaseq_t, FUN = sum)
rownames(rnaseq_collapsed) <- rnaseq_collapsed$sample_id
rnaseq_collapsed$sample_id <- NULL

rnaseq_final <- as.data.frame(t(rnaseq_collapsed))  # genes x samples
colnames(rnaseq_final) <- gsub("\\.", "-", colnames(rnaseq_final))

# ---------------------------------------------------------
# 5. Collapse duplicate immune columns (same sample ID) by MEAN
# ---------------------------------------------------------
collapse_dup_cols_mean <- function(mat, id_col = NULL) {
  if (!is.null(id_col)) {
    ids <- mat[[id_col]]
    mat_vals <- mat[, setdiff(colnames(mat), id_col)]
  } else {
    mat_vals <- mat
  }
  
  ucols <- unique(colnames(mat_vals))
  
  collapsed <- sapply(ucols, function(cn) {
    cols <- which(colnames(mat_vals) == cn)
    if (length(cols) > 1) {
      rowMeans(mat_vals[, cols, drop = FALSE], na.rm = TRUE)
    } else {
      mat_vals[, cols]
    }
  })
  
  collapsed <- as.data.frame(collapsed)
  colnames(collapsed) <- ucols
  
  if (!is.null(id_col)) {
    collapsed <- cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), id_col), collapsed)
  }
  
  return(collapsed)
}

immune_final <- collapse_dup_cols_mean(immune, id_col = "feature_id")

# ---------------------------------------------------------
# 6. Find common samples between RNA-seq and immune data
# ---------------------------------------------------------
common_samples <- intersect(colnames(rnaseq_final), colnames(immune_final)[-1])
cat("Number of common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found between rnaseq and immune data. Check barcode formats.")
}

# ---------------------------------------------------------
# 7. Subset & reorder both datasets to match sample order
# ---------------------------------------------------------
rnaseq_final <- rnaseq_final[, common_samples, drop = FALSE]
immune_final <- immune_final[, c("feature_id", common_samples), drop = FALSE]

# ---------------------------------------------------------
# 8. Verify alignment
# ---------------------------------------------------------
stopifnot(identical(colnames(rnaseq_final), colnames(immune_final)[-1]))
cat("Sample order matches:", identical(colnames(rnaseq_final), colnames(immune_final)[-1]), "\n")
cat("Final RNAseq dim:", dim(rnaseq_final), "\n")
cat("Final Immune dim:", dim(immune_final), "\n")

head(rnaseq_final)
head(immune_final)

# ---------------------------------------------------------
# 9. Move feature_id into rownames, then drop it
# ---------------------------------------------------------
rownames(immune_final) <- immune_final$feature_id
immune_final <- immune_final[, -1]
rnaseq_final=rnaseq_final[-1,]


head(immune_final [1:5,] )
head(rnaseq_final[1:5,])


# ---------------------------------------------------------
# 10. Save both as .rds files
# ---------------------------------------------------------
dir.create("RNAseq_Immune", showWarnings = FALSE)

saveRDS(rnaseq_final, file = "RNAseq_Immune/rnaseq_aligned.rds")
saveRDS(immune_final, file = "RNAseq_Immune/immune_aligned.rds")

# Optional: combined list object
saveRDS(list(rnaseq = rnaseq_final, immune = immune_final),
        file = "RNAseq_Immune/RNAseq_Immune_combined.rds")

cat("Saved files in RNAseq_Immune/:\n")
cat(" - rnaseq_aligned.rds\n")
cat(" - immune_aligned.rds\n")
cat(" - RNAseq_Immune_combined.rds\n")
