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
data <- read_tsv("TCGA.BRCA.sampleMap_HumanMethylation450.gz")
colnames(data)


####################
# =========================================================
# Complete pipeline: Align RNA-seq (DESeqDataSet) and 
# Methylation data by sample ID, then save as .rds
# =========================================================

library(DESeq2)

# ---------------------------------------------------------
# 1. Extract raw count matrix from DESeqDataSet
# ---------------------------------------------------------
rnaseq_counts <- counts(rnaseq)   # genes x samples matrix

cat("RNAseq dimensions:", dim(rnaseq_counts), "\n")
cat("Methylation/data dimensions:", dim(data), "\n")

# ---------------------------------------------------------
# 2. Truncate barcodes to 15-character sample IDs
#    (TCGA-XX-XXXX-01A -> TCGA-XX-XXXX-01)
# ---------------------------------------------------------
colnames(rnaseq_counts) <- substr(colnames(rnaseq_counts), 1, 15)
colnames(data)[-1]      <- substr(colnames(data)[-1], 1, 15)  # assumes col1 = "sample"/probe ID column

# ---------------------------------------------------------
# 3. Collapse duplicate RNA-seq columns (same sample ID) by SUM
#    (raw counts -> summing replicate aliquots is valid)
# ---------------------------------------------------------
rnaseq_t <- as.data.frame(t(rnaseq_counts))
rnaseq_t$sample_id <- rownames(rnaseq_t)

rnaseq_collapsed <- aggregate(. ~ sample_id, data = rnaseq_t, FUN = sum)
rownames(rnaseq_collapsed) <- rnaseq_collapsed$sample_id
rnaseq_collapsed$sample_id <- NULL

rnaseq_final <- as.data.frame(t(rnaseq_collapsed))  # back to genes x samples

# ---------------------------------------------------------
# 4. Collapse duplicate methylation columns (same sample ID) by MEAN
#    (beta-values are proportions -> averaging is valid, not summing)
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

data_final <- collapse_dup_cols_mean(data, id_col = "sample")

# ---------------------------------------------------------
# 5. Find common samples between RNA-seq and methylation data
# ---------------------------------------------------------

colnames(rnaseq_final)=gsub("\\.","-",colnames(rnaseq_final))

common_samples <- intersect(colnames(rnaseq_final), colnames(data_final)[-1])
cat("Number of common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found between rnaseq and methylation data. Check barcode formats.")
}

# ---------------------------------------------------------
# 6. Subset & reorder both datasets to match sample order
# ---------------------------------------------------------
rnaseq_final <- rnaseq_final[, common_samples, drop = FALSE]
data_final   <- data_final[, c("sample", common_samples), drop = FALSE]

# ---------------------------------------------------------
# 7. Verify alignment
# ---------------------------------------------------------
stopifnot(identical(colnames(rnaseq_final), colnames(data_final)[-1]))
cat("Sample order matches:", identical(colnames(rnaseq_final), colnames(data_final)[-1]), "\n")
cat("Final RNAseq dim:", dim(rnaseq_final), "\n")
cat("Final Methylation dim:", dim(data_final), "\n")


head(rnaseq_final)
head(data_final)
rownames(data_final)=data_final$sample


data_final=data_final[,-1]
rnaseq_final=rnaseq_final[-1,]


head(data_final)
head(rnaseq_final)
# ---------------------------------------------------------
# 8. Save both as .rds files
# ---------------------------------------------------------
dir.create("RNAseq_DNAMeth")

saveRDS(rnaseq_final, file = "RNAseq_DNAMeth/rnaseq_aligned.rds")
saveRDS(data_final,   file = "RNAseq_DNAMeth/methylation_aligned.rds")

cat("Saved: rnaseq_aligned.rds and methylation_aligned.rds\n")
