
rm(list = ls())

# =============================================================================
# Load required libraries
# =============================================================================
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(future)
library(gtools)
library(matrixStats)
library(data.table)
library(MASS)
library(sandwich)
library(lmtest)
library(parallel)
library(R.utils)
library(readxl)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(DESeq2)

# =============================================================================
# Set up directories
# =============================================================================
main_dir    <- "Result_viral"
result_dir  <- file.path(main_dir, "Results_per_viral")
plot_dir    <- file.path(main_dir, "Plots")
table_dir   <- file.path(main_dir, "Tables")

dir.create(main_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,  showWarnings = FALSE, recursive = TRUE)


genes=read.csv("Results_eQTM/Tables/gene_regions_250kb.csv")


# =============================================================================
# Load datasets
# =============================================================================
dds     <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")
rnaseq  <- readRDS("RNAseq_Viral/rnaseq_aligned.rds")
metab   <- readRDS("RNAseq_Viral/viral_aligned.rds")

head(genes)
head(rnaseq)
head(metab)

itk=which(rownames(rnaseq )%in%genes$NAME)
rnaseq=rnaseq[itk,]

pheno <- colData(dds) %>% data.frame()

cat("Original dimensions:\n")
cat("rnaseq:", dim(rnaseq), "\n")
cat("metab :", dim(metab), "\n")
cat("pheno :", dim(pheno), "\n")


table(pheno$condition)

# =============================================================================
# Ensure metab has samples as COLUMNS (genes/proteins as rows) for now,
# consistent with rnaseq, so we can use the same dedupe/truncate logic
# =============================================================================
# If metab is samples x proteins (samples as rows), transpose so proteins x samples
if (all(rownames(pheno) %in% colnames(metab)) == FALSE &&
    all(rownames(pheno) %in% rownames(metab))) {
  metab <- t(metab)
}

cat("\nBarcode length check (before truncation):\n")
cat("rnaseq colnames nchar:", unique(nchar(colnames(rnaseq))), "\n")
cat("metab  colnames nchar:", unique(nchar(colnames(metab))), "\n")
cat("pheno  rownames nchar:", unique(nchar(rownames(pheno))), "\n")

# =============================================================================
# Truncate all barcodes to 15 characters (TCGA-XX-XXXX-01) and deduplicate
# =============================================================================
truncate_barcode <- function(x, len = 15) substr(x, 1, len)

dedupe_by_id <- function(mat_or_df, ids, margin = 2) {
  dup <- duplicated(ids)
  n_dup <- sum(dup)
  if (n_dup > 0) {
    cat("  Found", n_dup, "duplicated IDs after truncation. Keeping first occurrence of each.\n")
  }
  keep <- !dup
  if (margin == 2) {
    mat_or_df <- mat_or_df[, keep, drop = FALSE]
    colnames(mat_or_df) <- ids[keep]
  } else {
    mat_or_df <- mat_or_df[keep, , drop = FALSE]
    rownames(mat_or_df) <- ids[keep]
  }
  mat_or_df
}

rnaseq_ids_short <- truncate_barcode(colnames(rnaseq))
metab_ids_short  <- truncate_barcode(colnames(metab))
pheno_ids_short  <- truncate_barcode(rownames(pheno))

cat("\n--- Deduplicating rnaseq (columns) ---\n")
rnaseq <- dedupe_by_id(rnaseq, rnaseq_ids_short, margin = 2)

cat("--- Deduplicating metab (columns) ---\n")
metab <- dedupe_by_id(metab, metab_ids_short, margin = 2)

cat("--- Deduplicating pheno (rows) ---\n")
pheno <- dedupe_by_id(pheno, pheno_ids_short, margin = 1)

# =============================================================================
# Align samples across all datasets
# =============================================================================
common_samples <- Reduce(intersect, list(colnames(rnaseq), colnames(metab), rownames(pheno)))
common_samples <- sort(common_samples)

cat("\nNumber of common samples found:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found even after truncation. Check barcode formats manually.")
}

rnaseq <- rnaseq[, common_samples, drop = FALSE]
metab  <- metab[, common_samples, drop = FALSE]
pheno  <- pheno[common_samples, , drop = FALSE]

stopifnot(all(colnames(rnaseq) == colnames(metab)))
stopifnot(all(colnames(rnaseq) == rownames(pheno)))

cat("\nFinal aligned dimensions:\n")
cat("rnaseq:", dim(rnaseq), "\n")
cat("metab :", dim(metab), "\n")
cat("pheno :", dim(pheno), "\n")

# =============================================================================
# Covariates
# =============================================================================
cov1 <- c("SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8")

rlm_formula <- as.formula(
  paste("as.numeric(cyto1) ~ expression +", paste(cov1, collapse = " + "))
)

# =============================================================================
# Prepare inputs for the loop
# =============================================================================
pheno3 <- pheno
pheno3$ID_idat <- rownames(pheno3)
pheno3$ID      <- rownames(pheno3)

df <- t(rnaseq) %>% as.data.frame()   # samples x genes, for indexing df[, gene]
genes <- colnames(df)

metab_df <- t(metab) %>% as.data.frame()  # samples x proteins
metab_df <- tibble::rownames_to_column(metab_df, var = "Row.names")

pheno.T0 <- pheno3

# =============================================================================
# Main loop: RNAseq gene x Protein association testing (rlm robust regression)
# =============================================================================
for (ikt in 1:(ncol(metab_df) - 1)) {
  
  print(ikt)
  
  id1 <- ikt + 1  # protein columns start from 2nd column (1st is Row.names)
  
  result <- pheno.T0
  result$SAM  <- as.numeric(metab_df[, id1])
  result$SAM1 <- log(result$SAM, 2)
  result$SAM2 <- qnorm((rank(result$SAM, na.last = "keep") - 0.5) / sum(!is.na(result$SAM)))
  
  TrImmRes.T0 <- result
  
  cor_test <- function(gene) {
    sub <- df[, gene] %>% data.frame()
    ind <- which(is.na(sub) == TRUE)
    colnames(sub) <- 'expression'
    sub$expression <- as.numeric(as.character(sub$expression))   # <-- force numeric
    
    cyto1 <- TrImmRes.T0$SAM2
    
    test.df <- cbind(pheno.T0, sub)
    test.df$cyto1 <- cyto1
    test.df$expression <- as.numeric(as.character(test.df$expression))  # <-- double-safe
    
    bad <- as.numeric(rep(NA, 4))
    names(bad) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    result <- bad
    
    if (length(ind) <= 50) {
      tryCatch(
        {
          ML <- rlm(rlm_formula, data = test.df)
          cf <- try(coeftest(ML, vcov = vcovHC(ML, type = "HC0")))
          
          x <<- x + 1
          cat("I am at ", x, " of total=", length(P_value_Index), ".\n", sep = "")
          cat("Left =", length(P_value_Index) - x, ".\n", sep = "")
          
          if (class(cf)[1] == "try-error") {
            result <- bad
          } else if (!"expression" %in% rownames(cf)) {
            # Safety net: if expression still becomes categorical, bail out
            result <- bad
          } else {
            result <- cf["expression", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
          }
          
          if (length(result) != 4) {
            result <- bad
          }
        },
        error = function(error_message) {
          message("This is my custom message.")
          message("And below is the error message from R:")
          message(error_message)
          return(bad)
        }
      )
    }
    result
  }
  
  
  P_value_Index <- genes
  x <- 0
  Result <- vapply(P_value_Index, cor_test, numeric(4))
  
  cor_info <- as.data.frame(t(Result))
  colnames(cor_info) <- c("Estimate", "StdError", "z-score", "pval")
  cor_info$Gene <- rownames(cor_info)
  cor_info$FDR <- p.adjust(cor_info$pval, method = "fdr")
  cor_info$Protein <- colnames(metab_df)[id1]
  
  output <- paste0(result_dir, "/cor_info.protein", ikt, ".rds")
  saveRDS(cor_info, output, compress = "xz")
  
  #rm(cor_info)
  gc()
}

cat("RNAseq-Protein rlm association analysis complete. Results saved in", result_dir, "\n")