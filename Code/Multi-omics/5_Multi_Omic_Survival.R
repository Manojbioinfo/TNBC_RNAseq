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

# =========================================================
# Complete pipeline: Align RNA-seq (DESeqDataSet) and 
# Survival data by sample ID, then save as .rds
# =========================================================

library(DESeq2)
library(readr)
library(dplyr)

# ---------------------------------------------------------
# 1. Load survival data (long format: samples as rows)
# ---------------------------------------------------------
survival <- read_tsv("survival_BRCA_survival.txt", show_col_types = FALSE)

cat("Survival dimensions:", dim(survival), "\n")

# ---------------------------------------------------------
# 2. Extract raw count matrix from DESeqDataSet
# ---------------------------------------------------------
rnaseq_counts <- counts(rnaseq)   # genes x samples matrix

cat("RNAseq dimensions:", dim(rnaseq_counts), "\n")

# ---------------------------------------------------------
# 3. Truncate barcodes to 15-character sample IDs
# ---------------------------------------------------------
colnames(rnaseq_counts) <- substr(colnames(rnaseq_counts), 1, 15)
survival$sample         <- substr(survival$sample, 1, 15)

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
# 5. Fix RNAseq colnames (dots -> dashes, if applicable)
# ---------------------------------------------------------
colnames(rnaseq_final) <- gsub("\\.", "-", colnames(rnaseq_final))

# ---------------------------------------------------------
# 6. Collapse duplicate survival rows (same sample ID)
#    Survival is per-patient, duplicates are rare, but if 
#    present, keep the first non-NA record per sample
# ---------------------------------------------------------
survival_collapsed <- survival %>%
  group_by(sample) %>%
  summarise(across(everything(), ~ dplyr::first(na.omit(.x))), .groups = "drop")

cat("Number of duplicate sample rows collapsed:", 
    nrow(survival) - nrow(survival_collapsed), "\n")

# ---------------------------------------------------------
# 7. Find common samples between RNA-seq and survival data
# ---------------------------------------------------------
common_samples <- intersect(colnames(rnaseq_final), survival_collapsed$sample)
cat("Number of common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found between rnaseq and survival data. Check barcode formats.")
}

# ---------------------------------------------------------
# 8. Subset & reorder both datasets to match sample order
# ---------------------------------------------------------
rnaseq_final <- rnaseq_final[, common_samples, drop = FALSE]

survival_final <- survival_collapsed %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

# ---------------------------------------------------------
# 9. Verify alignment
# ---------------------------------------------------------
stopifnot(identical(colnames(rnaseq_final), survival_final$sample))
cat("Sample order matches:", identical(colnames(rnaseq_final), survival_final$sample), "\n")
cat("Final RNAseq dim:", dim(rnaseq_final), "\n")
cat("Final Survival dim:", dim(survival_final), "\n")

head(rnaseq_final)
head(survival_final)

# ---------------------------------------------------------
# 10. Move "sample" column into rownames
# ---------------------------------------------------------
survival_final <- as.data.frame(survival_final)
rownames(survival_final) <- survival_final$sample
survival_final <- survival_final[, -1]
rnaseq_final=rnaseq_final[-1,]




# ---------------------------------------------------------
# 11. Save both as .rds files
# ---------------------------------------------------------
dir.create("RNAseq_Survival", showWarnings = FALSE)

saveRDS(rnaseq_final,   file = "RNAseq_Survival/rnaseq_aligned.rds")
saveRDS(survival_final, file = "RNAseq_Survival/survival_aligned.rds")

# Optional: combined list object
saveRDS(list(rnaseq = rnaseq_final, survival = survival_final),
        file = "RNAseq_Survival/RNAseq_Survival_combined.rds")

cat("Saved files in RNAseq_Survival/:\n")
cat(" - rnaseq_aligned.rds\n")
cat(" - survival_aligned.rds\n")
cat(" - RNAseq_Survival_combined.rds\n")