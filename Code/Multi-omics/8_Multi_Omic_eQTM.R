# =============================================================================
# Load required libraries
# =============================================================================
library(limma)
library(edgeR)
library(bacon)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # switch to 450k anno if needed

# =============================================================================
# Set up output directory structure
# =============================================================================
#250000 

flank <- 250000  
flank1 <- "10kb"


main_dir    <- paste0("Results_eQTM")
result_dir  <- file.path(main_dir, "Results_per_CpG")
plot_dir    <- file.path(main_dir, "Plots")
table_dir   <- file.path(main_dir, "Tables")

dir.create(main_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,  showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load datasets
# =============================================================================
dds     <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")
rnaseq  <- readRDS("RNAseq_DNAMeth/rnaseq_aligned.rds")
methyl  <- readRDS("RNAseq_DNAMeth/methylation_aligned.rds")

pheno <- colData(dds) %>% data.frame()

# Remove all-NA CpG rows
methyl <- methyl[rowSums(!is.na(methyl)) > 0, ]

cat("Original dimensions:\n")
cat("rnaseq:", dim(rnaseq), "\n")
cat("methyl:", dim(methyl), "\n")
cat("pheno :", dim(pheno), "\n")

cat("\nBarcode length check (before truncation):\n")
cat("rnaseq colnames nchar:", unique(nchar(colnames(rnaseq))), "\n")
cat("methyl colnames nchar:", unique(nchar(colnames(methyl))), "\n")
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
methyl_ids_short <- truncate_barcode(colnames(methyl))
pheno_ids_short  <- truncate_barcode(rownames(pheno))

cat("\n--- Deduplicating rnaseq (columns) ---\n")
rnaseq <- dedupe_by_id(rnaseq, rnaseq_ids_short, margin = 2)

cat("--- Deduplicating methyl (columns) ---\n")
methyl <- dedupe_by_id(methyl, methyl_ids_short, margin = 2)

cat("--- Deduplicating pheno (rows) ---\n")
pheno <- dedupe_by_id(pheno, pheno_ids_short, margin = 1)

# =============================================================================
# Align samples across all datasets
# =============================================================================
common_samples <- Reduce(intersect, list(colnames(rnaseq), colnames(methyl), rownames(pheno)))
common_samples <- sort(common_samples)

cat("\nNumber of common samples found:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common samples found even after truncation. Check barcode formats manually.")
}

rnaseq <- rnaseq[, common_samples, drop = FALSE]
methyl <- methyl[, common_samples, drop = FALSE]
pheno  <- pheno[common_samples, , drop = FALSE]

stopifnot(all(colnames(rnaseq) == colnames(methyl)))
stopifnot(all(colnames(rnaseq) == rownames(pheno)))

cat("\nFinal aligned dimensions:\n")
cat("rnaseq:", dim(rnaseq), "\n")
cat("methyl:", dim(methyl), "\n")
cat("pheno :", dim(pheno), "\n")

# =============================================================================
# Check NA rates in methylation data
# =============================================================================
na_rates <- rowMeans(is.na(methyl))
cat("\nSummary of NA rate per CpG:\n")
print(summary(na_rates))

png(file.path(plot_dir, "CpG_NA_rate_histogram.png"), width = 800, height = 600)
hist(na_rates, breaks = 50, main = "NA rate per CpG")
dev.off()

cat("Number of CpGs with >50% missing:", sum(na_rates > 0.5), "\n")
cat("Number of CpGs with >20% missing:", sum(na_rates > 0.2), "\n")

# =============================================================================
# Save aligned datasets
# =============================================================================
saveRDS(pheno,  file.path(main_dir, "aligned_pheno.rds"))
saveRDS(rnaseq, file.path(main_dir, "aligned_rnaseq.rds"))
saveRDS(methyl, file.path(main_dir, "aligned_methyl.rds"))

cat("\nSaved aligned_pheno.rds, aligned_rnaseq.rds, aligned_methyl.rds to", main_dir, "\n")

# =============================================================================
# Gene vector and 250kb region definition
# =============================================================================
gene_vector <- c("MMP1", "PRAME", "HORMAD1", "NEK2", "KIF18B", "NMU")

gene_loc <- read_tsv("NCBI37.3.gene.loc", show_col_types = FALSE)

 

gene_regions <- gene_loc %>%
  filter(NAME %in% gene_vector) %>%
  mutate(
    REGION_START = pmax(0, START - flank),
    REGION_END   = END + flank
  ) %>%
  select(NAME, CHR, START, END, REGION_START, REGION_END, STRAND)

write.csv(gene_regions, file.path(table_dir, "gene_regions_250kb.csv"), row.names = FALSE)

# =============================================================================
# Map CpGs to gene regions using EPIC annotation
# =============================================================================
ann850k   <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- data.frame(ann850k[c("chr", "pos")])
annotdata$CpG <- rownames(annotdata)
annotdata$CHR <- gsub("chr", "", annotdata$chr)

# Restrict to CpGs present in methyl object
annotdata <- annotdata %>% filter(CpG %in% rownames(methyl))

cpg_gene_map <- lapply(seq_len(nrow(gene_regions)), function(i) {
  g <- gene_regions[i, ]
  annotdata %>%
    filter(CHR == g$CHR, pos >= g$REGION_START, pos <= g$REGION_END) %>%
    mutate(GENE = g$NAME)
}) %>% bind_rows()

write.csv(cpg_gene_map, file.path(table_dir, "cpg_gene_map.csv"), row.names = FALSE)

cat("Number of CpGs mapped per gene:\n")
print(table(cpg_gene_map$GENE))

# =============================================================================
# Covariates
# =============================================================================
cov1 <- c("SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8")

# Check covariates for NAs
cat("\nMissing values per covariate in pheno:\n")
print(sapply(pheno[, cov1], function(x) sum(is.na(x))))

# =============================================================================
# Bacon correction function
# =============================================================================
correction <- function(toptable) {
  bc <- bacon(effectsizes = as.matrix(toptable$logFC),
              standarderrors = as.matrix(toptable$logFC / toptable$t))
  toptable$logFC.cor     <- es(bc)
  toptable$t.cor         <- tstat(bc)
  toptable$P.Value.cor   <- pval(bc)
  toptable$adj.P.Val.cor <- p.adjust(toptable$P.Value.cor, method = "BH")
  return(toptable)
}

# =============================================================================
# eQTM function
# =============================================================================
do.twas <- function(cpg, covariates, sample_data, counts_data, dnam, filter = 0.8, min_frac = 0.7) {
  
  sample_data <- as.data.frame(sample_data)
  
  if (!cpg %in% rownames(dnam)) {
    stop(paste("CpG", cpg, "not found in methylation data"))
  }
  
  meth_vals <- dnam[cpg, ]
  meth_vals <- meth_vals[match(rownames(sample_data), names(meth_vals))]
  
  sample_data[[cpg]] <- as.numeric(meth_vals)
  
  vars <- sample_data[, c(cpg, covariates), drop = FALSE]
  vars <- na.omit(vars)
  
  min_required <- max(10, floor(min_frac * nrow(sample_data)))
  if (nrow(vars) < min_required) {
    stop(paste0("Too few samples with complete data for CpG ", cpg,
                " (", nrow(vars), "/", nrow(sample_data), ")"))
  }
  
  counts_sub <- counts_data[, match(rownames(vars), colnames(counts_data)), drop = FALSE]
  
  design <- model.matrix(~ ., vars)
  
  y <- DGEList(counts = counts_sub)
  y <- y[rowSums(y$counts > 0) > filter * ncol(y), ]
  y <- calcNormFactors(y)
  v <- voom(y, design)
  
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Explicit namespace to avoid masking issues
  results <- limma::topTable(fit, coef = 2, n = Inf)
  results$SYMBOL <- rownames(results)
  return(results)
}
# =============================================================================
# Run eQTM analysis for all mapped CpGs
# =============================================================================
cpg_list <- unique(cpg_gene_map$CpG)

cat("Running eQTM for", length(cpg_list), "CpGs...\n")

for (i in seq_along(cpg_list)) {
  cpg_i <- cpg_list[i]
  
  twas <- tryCatch(
    do.twas(cpg_i, cov1, pheno, rnaseq, methyl, filter = 0.8, min_frac = 0.7),
    error = function(e) {
      message("Error for CpG ", cpg_i, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(twas)) next
  
  twas <- correction(twas)
  
  gene_hit <- cpg_gene_map$GENE[cpg_gene_map$CpG == cpg_i]
  
  tdata <- data.frame(
    CpG        = cpg_i,
    Gene       = twas$SYMBOL,
    t          = twas$t,
    pvalue     = twas$P.Value,
    adj.pvalue = twas$adj.P.Val,
    t.cor      = twas$t.cor,
    pvalue.cor = twas$P.Value.cor,
    adj.pvalue.cor = twas$adj.P.Val.cor,
    TargetGene = paste(gene_hit, collapse = ";"),
    num        = nrow(pheno)
  )
  
  saveRDS(tdata, file.path(result_dir, paste0(cpg_i, ".rds")))
  
  if (i %% 50 == 0) cat("Processed", i, "/", length(cpg_list), "CpGs\n")
}

cat("eQTM analysis complete.\n")

# =============================================================================
# Combine all results into one table
# =============================================================================
result_files <- list.files(result_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No result files were generated. Check the errors above (likely all CpGs failed the sample-completeness filter).")
}

all_results  <- lapply(result_files, readRDS) %>% bind_rows()

saveRDS(all_results, file.path(main_dir, "all_eqtm_results_combined.rds"))
write.csv(all_results, file.path(table_dir, "all_eqtm_results_combined.csv"), row.names = FALSE)

# =============================================================================
# Filter to target genes of interest (cis-eQTM: CpG's target gene expression)
# =============================================================================
final_results <- all_results %>%
  filter(Gene %in% gene_vector)

final_results <- final_results %>% arrange(pvalue)

write.csv(final_results, file.path(table_dir, "eqtm_results_target_genes.csv"), row.names = FALSE)
saveRDS(final_results, file.path(main_dir, "eqtm_results_target_genes.rds"))

cat("Top eQTM hits for target genes:\n")
print(head(final_results, 20))

# =============================================================================
# Summary table: number of significant CpGs per gene
# =============================================================================
summary_table <- final_results %>%
  group_by(Gene) %>%
  summarise(
    n_CpGs_tested = n(),
    n_sig_p05     = sum(pvalue < 0.05, na.rm = TRUE),
    n_sig_fdr05   = sum(adj.pvalue < 0.05, na.rm = TRUE),
    n_sig_p05_cor = sum(pvalue.cor < 0.05, na.rm = TRUE)
  )

write.csv(summary_table, file.path(table_dir, "eqtm_summary_by_gene.csv"), row.names = FALSE)
print(summary_table)

# =============================================================================
# Plots: -log10(P) vs CpG position for each target gene
# =============================================================================
plot_data <- final_results %>%
  left_join(cpg_gene_map, by = c("CpG" = "CpG")) %>%
  filter(!is.na(pos))

for (g in gene_vector) {
  df_g <- plot_data %>% filter(GENE == g)
  if (nrow(df_g) == 0) next
  
  p <- ggplot(df_g, aes(x = pos, y = -log10(pvalue))) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(
      title = paste("eQTM results for", g, "(±250kb region)"),
      x = "CpG position (hg19)",
      y = "-log10(P-value)"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(plot_dir, paste0("eQTM_", g, "_region_plot.png")),
    plot = p, width = 8, height = 5, dpi = 300
  )
}

cat("All plots and tables saved under:", main_dir, "\n")

