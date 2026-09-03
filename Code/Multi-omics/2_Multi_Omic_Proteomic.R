rm(list = ls())
# Load required libraries
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)

rnaseq <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")

colnames(rnaseq)


# Read the data
# =========================================================
# Complete pipeline: Align RNA-seq (DESeqDataSet) and 
# Proteomic (RPPA) data by sample ID, then save as .rds
# =========================================================

library(DESeq2)
library(readr)

# ---------------------------------------------------------
# 1. Load proteomic (RPPA) data
# ---------------------------------------------------------
proteomic <- read_tsv("RPPA", show_col_types = FALSE)

cat("Proteomic dimensions:", dim(proteomic), "\n")

# ---------------------------------------------------------
# 2. Extract raw count matrix from DESeqDataSet
# ---------------------------------------------------------
rnaseq_counts <- counts(rnaseq)   # genes x samples matrix

cat("RNAseq dimensions:", dim(rnaseq_counts), "\n")

# ---------------------------------------------------------
# 3. Truncate barcodes to 15-character sample IDs
# ---------------------------------------------------------
colnames(rnaseq_counts)  <- substr(colnames(rnaseq_counts), 1, 15)
colnames(proteomic)[-1]  <- substr(colnames(proteomic)[-1], 1, 15)

# ---------------------------------------------------------
# 4. Collapse duplicate RNA-seq columns (same sample ID) by SUM
# ---------------------------------------------------------
rnaseq_t <- as.data.frame(t(rnaseq_counts))
rnaseq_t$sample_id <- rownames(rnaseq_t)

rnaseq_collapsed <- aggregate(. ~ sample_id, data = rnaseq_t, FUN = sum)
rownames(rnaseq_collapsed) <- rnaseq_collapsed$sample_id
rnaseq_collapsed$sample_id <- NULL

rnaseq_final <- as.data.frame(t(rnaseq_collapsed))  # genes x samples

# ---------------------------------------------------------
# 5. Collapse duplicate proteomic columns (same sample ID) by MEAN
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

proteomic_final <- collapse_dup_cols_mean(proteomic, id_col = "sample")

# ---------------------------------------------------------
# 6. Fix RNAseq colnames (dots -> dashes, if applicable)
# ---------------------------------------------------------
colnames(rnaseq_final) <- gsub("\\.", "-", colnames(rnaseq_final))

# ---------------------------------------------------------
# 7. Find common samples between RNA-seq and proteomic data
# ---------------------------------------------------------
common_samples <- intersect(colnames(rnaseq_final), colnames(proteomic_final)[-1])
cat("Number of common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found between rnaseq and proteomic data. Check barcode formats.")
}

# ---------------------------------------------------------
# 8. Subset & reorder both datasets to match sample order
# ---------------------------------------------------------
rnaseq_final    <- rnaseq_final[, common_samples, drop = FALSE]
proteomic_final <- proteomic_final[, c("sample", common_samples), drop = FALSE]

# ---------------------------------------------------------
# 9. Verify alignment
# ---------------------------------------------------------
stopifnot(identical(colnames(rnaseq_final), colnames(proteomic_final)[-1]))
cat("Sample order matches:", identical(colnames(rnaseq_final), colnames(proteomic_final)[-1]), "\n")
cat("Final RNAseq dim:", dim(rnaseq_final), "\n")
cat("Final Proteomic dim:", dim(proteomic_final), "\n")

head(rnaseq_final)
head(proteomic_final)

# ---------------------------------------------------------
# 10. Move "sample" column into rownames, then drop it
# ---------------------------------------------------------
rownames(proteomic_final) <- proteomic_final$sample
proteomic_final <- proteomic_final[, -1]


rnaseq_final=rnaseq_final[-1,]

head(proteomic_final[1:5,] )
head(rnaseq_final[1:5,])
# ---------------------------------------------------------
# 11. Save both as .rds files
# ---------------------------------------------------------
dir.create("RNAseq_Proteomic", showWarnings = FALSE)

saveRDS(rnaseq_final,    file = "RNAseq_Proteomic/rnaseq_aligned.rds")
saveRDS(proteomic_final, file = "RNAseq_Proteomic/proteomic_aligned.rds")

# Optional: combined list object
saveRDS(list(rnaseq = rnaseq_final, proteomic = proteomic_final),
        file = "RNAseq_Proteomic/RNAseq_Proteomic_combined.rds")

cat("Saved files in RNAseq_Proteomic/:\n")
cat(" - rnaseq_aligned.rds\n")
cat(" - proteomic_aligned.rds\n")
cat(" - RNAseq_Proteomic_combined.rds\n")

