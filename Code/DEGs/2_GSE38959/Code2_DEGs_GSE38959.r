# ============================================================================
# LOAD LIBRARIES
# ============================================================================

library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(tibble)

# ============================================================================
# SETUP
# ============================================================================

GEO_ID  <- "GSE38959"
out_dir <- "."

# Load raw data
vsd_pre <- readRDS("GSE38959_raw.rds")

count_filtered <- assay(vsd_pre)
coldata_pre    <- colData(vsd_pre) |> as.data.frame()
coldata_pre$condition <- coldata_pre$sample_type
coldata_pre$sample_id <- rownames(coldata_pre)


# ============================================================================
# 2 — VST Transformation
# ============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(count_filtered),
  colData   = coldata_pre,
  design    = ~ condition
)


# ============================================================================
# 8 — SAVE CLEAN DESEQ2 OBJECT
# ============================================================================


saveRDS(dds, paste0(GEO_ID,"_clean.rds"))

dds=readRDS("GSE38959_clean.rds")
table(dds$condition)
range(sort(dds$age_numeric))
