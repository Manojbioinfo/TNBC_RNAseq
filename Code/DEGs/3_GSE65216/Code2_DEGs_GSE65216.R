# ============================================================
# LOAD LIBRARIES
# ============================================================

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# ============================================================
# SETUP
# ============================================================

GEO_ID   <- "GSE65216"
out_dir  <- "."
data     <- readRDS("GSE65216_raw.rds")

coldata_pre <- colData(data)
coldata_pre_df <- as.data.frame(coldata_pre)
coldata_pre_df$sample_id <- rownames(coldata_pre_df)


# ============================================================
# 2 — VST Transformation
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(assay(data)),
  colData   = coldata_pre_df,
  design    = ~ condition
)


saveRDS(dds, file.path(out_dir, paste0(GEO_ID, "_dds_clean.rds")))

dds=readRDS("GSE65216_dds_clean.rds")
table(dds$condition)
