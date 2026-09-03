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
# Blood/Immune cell fraction (CIBERSORT-style) data by sample ID
# =========================================================

library(DESeq2)
library(readr)
library(dplyr)

# ---------------------------------------------------------
# 1. Load blood cell fraction data
# ---------------------------------------------------------
blood <- readRDS("Z:/1PERSONALprojectcomplete/TCGAData/Immune_cells/TCGA_Bloodcount.rds")
blood <- as.data.frame(blood)

cat("Blood dimensions:", dim(blood), "\n")

# ---------------------------------------------------------
# 2. Fix sample IDs: dots -> dashes, truncate to 15 chars
# ---------------------------------------------------------
blood$SampleID <- gsub("\\.", "-", blood$SampleID)
blood$SampleID <- substr(blood$SampleID, 1, 15)
head(blood)
# ---------------------------------------------------------
# 3. Extract raw count matrix from DESeqDataSet
# ---------------------------------------------------------
rnaseq_counts <- counts(rnaseq)   # genes x samples matrix
colnames(rnaseq_counts) <- substr(colnames(rnaseq_counts), 1, 15)

cat("RNAseq dimensions:", dim(rnaseq_counts), "\n")

# ---------------------------------------------------------
# 4. Collapse duplicate RNA-seq columns (same sample ID) by SUM
# ---------------------------------------------------------
rnaseq_t <- as.data.frame(t(rnaseq_counts))
rnaseq_t$sample_id <- rownames(rnaseq_t)

rnaseq_collapsed <- aggregate(. ~ sample_id, data = rnaseq_t, FUN = sum)
rownames(rnaseq_collapsed) <- rnaseq_collapsed$sample_id
#rnaseq_collapsed$sample_id <- NULL

rnaseq_final <- as.data.frame(t(rnaseq_collapsed))  # genes x samples
colnames(rnaseq_final) <- gsub("\\.", "-", colnames(rnaseq_final))

# ---------------------------------------------------------
# 5. Collapse duplicate blood rows (same SampleID) by MEAN
#    (numeric columns only; keep CancerType as first value)
# ---------------------------------------------------------
numeric_cols <- sapply(blood, is.numeric)

blood_numeric <- blood %>%
  group_by(SampleID) %>%
  summarise(across(all_of(names(numeric_cols)[numeric_cols]), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

blood_meta <- blood %>%
  group_by(SampleID) %>%
  summarise(CancerType = CancerType[1], .groups = "drop")

blood_collapsed <- left_join(blood_meta, blood_numeric, by = "SampleID")

# ---------------------------------------------------------
# 6. Find common samples between RNA-seq and blood data
# ---------------------------------------------------------
common_samples <- intersect(colnames(rnaseq_final), blood_collapsed$SampleID)
cat("Number of common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found between rnaseq and blood data. Check barcode formats.")
}

# ---------------------------------------------------------
# 7. Subset & reorder both datasets to match sample order
# ---------------------------------------------------------
rnaseq_final <- rnaseq_final[, common_samples, drop = FALSE]

blood_final <- blood_collapsed %>%
  filter(SampleID %in% common_samples) %>%
  arrange(match(SampleID, common_samples))

# ---------------------------------------------------------
# 8. Verify alignment
# ---------------------------------------------------------
stopifnot(identical(colnames(rnaseq_final), blood_final$SampleID))
cat("Sample order matches:", identical(colnames(rnaseq_final), blood_final$SampleID), "\n")
cat("Final RNAseq dim:", dim(rnaseq_final), "\n")
cat("Final Blood dim:", dim(blood_final), "\n")

head(rnaseq_final)
head(blood_final)

# ---------------------------------------------------------
# 9. Move SampleID into rownames
# ---------------------------------------------------------
blood_final <- as.data.frame(blood_final)
rownames(blood_final) <- blood_final$SampleID
blood_final$SampleID <- NULL

# ---------------------------------------------------------
# 10. Save both as .rds files
# ---------------------------------------------------------
dir.create("RNAseq_Blood", showWarnings = FALSE)

colnames(blood_final)

blood_final=blood_final[,-c(1,24:26)]

rnaseq_final=rnaseq_final[-1,]


blood_final1=as.data.frame(t(blood_final))

head(blood_final1[1:5,] )
head(rnaseq_final[1:5,])



saveRDS(rnaseq_final, file = "RNAseq_Blood/rnaseq_aligned.rds")
saveRDS(blood_final1, file = "RNAseq_Blood/blood_aligned.rds")

# Optional: combined list object
saveRDS(list(rnaseq = rnaseq_final, blood = blood_final1),
        file = "RNAseq_Blood/RNAseq_Blood_combined.rds")

cat("Saved files in RNAseq_Blood/:\n")
cat(" - rnaseq_aligned.rds\n")
cat(" - blood_aligned.rds\n")
cat(" - RNAseq_Blood_combined.rds\n")

