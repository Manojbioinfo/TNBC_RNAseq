rm(list = ls())

library(tidyverse)
library(DESeq2)

# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------
rnaseq <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")
viral  <- readRDS("TCGA_VIRAL.rds")   # adjust path as needed

cat("Viral dimensions:", dim(viral), "\n")

# ---------------------------------------------------------
# 2. Extract raw count matrix from DESeqDataSet
# ---------------------------------------------------------
rnaseq_counts <- counts(rnaseq)

# ---------------------------------------------------------
# 3. Truncate barcodes to 15-character sample IDs
# ---------------------------------------------------------
colnames(rnaseq_counts) <- substr(colnames(rnaseq_counts), 1, 15)
viral$AliquotBarcode    <- substr(viral$AliquotBarcode, 1, 15)

# ---------------------------------------------------------
# 4. Collapse duplicate RNA-seq columns (same sample ID) by SUM
# ---------------------------------------------------------
rnaseq_t <- as.data.frame(t(rnaseq_counts))
rnaseq_t$sample_id <- rownames(rnaseq_t)

rnaseq_collapsed <- aggregate(. ~ sample_id, data = rnaseq_t, FUN = sum)
rownames(rnaseq_collapsed) <- rnaseq_collapsed$sample_id
rnaseq_collapsed$sample_id <- NULL

rnaseq_final <- as.data.frame(t(rnaseq_collapsed))
colnames(rnaseq_final) <- gsub("\\.", "-", colnames(rnaseq_final))

# ---------------------------------------------------------
# 5. Prepare viral data: drop metadata (keep only sample + numeric viral cols)
# ---------------------------------------------------------
meta_cols_to_drop <- c("ParticipantBarcode", "SampleBarcode", "Study", "SampleTypeLetterCode")

head(viral)

# Force correct namespace
viral_numeric <- viral %>%
  dplyr::select(-any_of(c("ParticipantBarcode", "SampleBarcode", "Study", "SampleTypeLetterCode"))) %>%
  dplyr::rename(sample = AliquotBarcode)

head(viral_numeric)
# ---------------------------------------------------------
# 6. Collapse duplicate sample IDs (mean of numeric columns)
# ---------------------------------------------------------
viral_collapsed <- viral_numeric %>%
  group_by(sample) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# ---------------------------------------------------------
# 7. Transpose to features x samples (to match rnaseq_final orientation)
# ---------------------------------------------------------
viral_mat <- as.data.frame(t(viral_collapsed[,-1]))
colnames(viral_mat) <- viral_collapsed$sample
viral_mat$feature <- rownames(viral_mat)
viral_final <- viral_mat %>% select(feature, everything())

# ---------------------------------------------------------
# 8. Find common samples between RNA-seq and viral data
# ---------------------------------------------------------
common_samples <- intersect(colnames(rnaseq_final), colnames(viral_final)[-1])
cat("Number of common samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found between rnaseq and viral data. Check barcode formats.")
}

# ---------------------------------------------------------
# 9. Subset & reorder both datasets to match sample order
# ---------------------------------------------------------
rnaseq_final <- rnaseq_final[, common_samples, drop = FALSE]
viral_final  <- viral_final[, c("feature", common_samples), drop = FALSE]

# ---------------------------------------------------------
# 10. Verify alignment
# ---------------------------------------------------------
stopifnot(identical(colnames(rnaseq_final), colnames(viral_final)[-1]))
cat("Sample order matches:", identical(colnames(rnaseq_final), colnames(viral_final)[-1]), "\n")
cat("Final RNAseq dim:", dim(rnaseq_final), "\n")
cat("Final Viral dim:", dim(viral_final), "\n")

# ---------------------------------------------------------
# 11. Move "feature" column into rownames, then drop it
# ---------------------------------------------------------
rownames(viral_final) <- viral_final$feature
viral_final <- viral_final[, -1]

head(rnaseq_final[1:5, 1:5])
head(viral_final[1:5,])

# ---------------------------------------------------------
# 12. Save both as .rds files
# ---------------------------------------------------------
dir.create("RNAseq_Viral", showWarnings = FALSE)

saveRDS(rnaseq_final, file = "RNAseq_Viral/rnaseq_aligned.rds")
saveRDS(viral_final,  file = "RNAseq_Viral/viral_aligned.rds")
saveRDS(list(rnaseq = rnaseq_final, viral = viral_final),
        file = "RNAseq_Viral/RNAseq_Viral_combined.rds")

cat("Saved files in RNAseq_Viral/:\n")
cat(" - rnaseq_aligned.rds\n")
cat(" - viral_aligned.rds\n")
cat(" - RNAseq_Viral_combined.rds\n")

