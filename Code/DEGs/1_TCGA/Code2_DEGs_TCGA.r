# ============================================================
# LOAD LIBRARIES
# ============================================================
Sys.setenv(http_proxy = "http://proxy.mh-hannover.de:8080")
Sys.setenv(https_proxy = "http://proxy.mh-hannover.de:8080")
options(download.file.method = "curl")
options(download.file.extra = "-L --proxy http://proxy.mh-hannover.de:8080")

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


vsd_pre <- readRDS("TCGA_raw.rds")
##### make sure we have the gene_symbol in the rowname
vsd_pre

# ============================================================
# Convert Ensembl IDs to gene symbols
# ============================================================

library(biomaRt)
library(tibble)

# Extract Ensembl IDs without version
ensembl_ids <- gsub("\\..*$", "", rownames(vsd_pre))

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene symbols
gene_map <- getBM(
  filters    = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values     = ensembl_ids,
  mart       = mart
)

# Match order to original Ensembl IDs
gene_map <- gene_map[match(ensembl_ids, gene_map$ensembl_gene_id), ]


# Get gene symbols
gene_symbols <- gene_map$hgnc_symbol

# Keep only rows with valid gene symbols
valid_rows <- !is.na(gene_symbols) & gene_symbols != ""
vsd_pre <- vsd_pre[valid_rows, ]

# Update rownames to gene symbols
rownames(vsd_pre) <- gene_symbols[valid_rows]

# Make rownames unique
rownames(vsd_pre) <- make.unique(rownames(vsd_pre))



# ============================================================
# 1 — Sample metadata
# ============================================================

sample_info_pre <- as.data.frame(colData(vsd_pre))
sample_info_pre$sample_id <- rownames(sample_info_pre)
sample_info_pre$age <- sample_info_pre$age_at_diagnosis / 365

# identify males
male_samples <- sample_info_pre$sample_id[sample_info_pre$gender == "male"]

# ============================================================
# 2 — Remove male samples from VSD object
# ============================================================

vsd_pre_female <- vsd_pre[, !(colnames(vsd_pre) %in% male_samples)]

# ============================================================
# 3 — Filter metadata to match remaining samples
# ============================================================

sample_info_pre_female <- sample_info_pre[
  sample_info_pre$sample_id %in% colnames(vsd_pre_female),
]

# enforce correct order (VERY IMPORTANT)
sample_info_pre_female <- sample_info_pre_female[
  match(colnames(vsd_pre_female), sample_info_pre_female$sample_id),
]

# ============================================================
# 4 — Create clean DESeq2 object (from VST if needed)
#    ⚠️ BEST PRACTICE: use raw counts if available
# ============================================================

dds_clean <- DESeqDataSetFromMatrix(
  countData = round(assay(vsd_pre_female)),
  colData   = sample_info_pre_female,
  design    = ~ condition
)

# ============================================================
# 5 — Save
# ============================================================
GEO_ID="TCGA"
saveRDS(dds_clean, file.path(paste0(GEO_ID, "_dds_clean.rds")))

dds_clean=readRDS("TCGA_dds_clean.rds")
table(dds_clean$condition)
range(sort(dds_clean$age))
