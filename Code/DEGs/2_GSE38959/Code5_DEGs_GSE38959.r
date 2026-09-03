# ============================================================
# LOAD LIBRARIES
# ============================================================
Sys.setenv(http_proxy = "http://proxy.mh-hannover.de:8080")
Sys.setenv(https_proxy = "http://proxy.mh-hannover.de:8080")
options(download.file.method = "curl")
options(download.file.extra = "-L --proxy http://proxy.mh-hannover.de:8080")
library(DESeq2)
library(sva)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(gridExtra)

# ============================================================
# SETUP
# ============================================================

# =============================================================================
# Load libraries
# =============================================================================
library(DESeq2)

# =============================================================================
# Load data
# =============================================================================
dds <- readRDS("GSE38959_dds_clean_with_SVs.rds")

# =============================================================================
# Check condition levels
# =============================================================================
table(dds$condition)

# Set "healthy" as reference (adjust name if needed)
dds$condition <- relevel(dds$condition, ref = "Healthy")

# =============================================================================
# Set design formula (include 8 surrogate variables)
# =============================================================================
design(dds) <- ~ condition+SV1 + SV2 + SV3

# =============================================================================
# Run DESeq2
# =============================================================================
dds <- DESeq(dds)

# =============================================================================
# Extract results: TNBC vs healthy
# =============================================================================
res <- results(dds, contrast = c("condition", "TNBC", "Healthy"))

# =============================================================================
# Order results by adjusted p-value
# =============================================================================
res <- res[order(res$padj), ]

# =============================================================================
# Summary of results
# =============================================================================
summary(res)

# =============================================================================
# Filter significant DEGs
# =============================================================================
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
res=as.data.frame(res)
res$Gene=rownames(res)
res_sig=as.data.frame(res_sig)
res_sig$Gene=rownames(res_sig)


# For significant results
res_sig$direction <- ifelse(res_sig$log2FoldChange > 0, "Up",
                            ifelse(res_sig$log2FoldChange < 0, "Down", "NoChange"))

table(res_sig$direction)


# =============================================================================
# Save results
# =============================================================================
write.csv(res, "DEGs_GSE38959_vs_Healthy_with_SV.csv")
write.csv(as.data.frame(res_sig), "DEGs_GSE38959_vs_Healthy_with_SV_significant.csv")

# =============================================================================
# Optional checks
# =============================================================================
resultsNames(dds)
head(res)
