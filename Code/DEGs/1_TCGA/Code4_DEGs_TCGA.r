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


dds<- readRDS("TCGA_dds_clean.rds")


flagged_samples=readRDS("TCGA.flagged_samples.rds")
dds
head(flagged_samples)


# Get list of samples to remove
samples_to_remove <- flagged_samples$sample_id

# Keep only samples NOT in flagged list
dds<- dds[, !(colnames(dds) %in% samples_to_remove)]



######## sva


# Build design matrices for SVA
mod  <- model.matrix(~ condition, data = colData(dds))
mod0 <- model.matrix(~ 1, data = colData(dds))

# Transform counts
vst_obj <- vst(dds, blind = TRUE)

# Extract expression matrix
vst_matrix <- assay(vst_obj)

# ============================================================================
# 3 — Run SVA
# ============================================================================
svobj <- sva(as.matrix(vst_matrix), mod, mod0, n.sv = NULL)
n.sv <- svobj$n.sv
cat("Estimated surrogate variables:", n.sv, "\n")
GEO_ID="TCGA"
# Save SV matrix
sv_df <- as.data.frame(svobj$sv)

colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
rownames(sv_df) <- colnames(vst_matrix)

all(rownames(sv_df) == colnames(dds))


head(rownames(sv_df))
head(colnames(dds))

colData(dds) <- cbind(colData(dds), sv_df)

dim(sv_df)              # should be samples × SVs
head(sv_df[,1:3])
any(is.na(sv_df))       # should be FALSE



saveRDS(dds, file.path( paste0(GEO_ID, "_dds_clean_with_SVs.rds")))
