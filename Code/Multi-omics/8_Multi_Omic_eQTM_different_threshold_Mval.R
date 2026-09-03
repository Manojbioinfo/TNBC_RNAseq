# =============================================================================
# Wrapper: run eQTM pipeline across multiple flank sizes
# =============================================================================
library(limma)
library(edgeR)
library(bacon)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(lumi)
# =============================================================================
# Load datasets ONCE (outside the loop - these don't change across flank sizes)
# =============================================================================
dds     <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")
rnaseq  <- readRDS("RNAseq_DNAMeth/rnaseq_aligned.rds")
methyl  <- readRDS("RNAseq_DNAMeth/methylation_aligned.rds")

methyl =beta2m(methyl )




pheno <- colData(dds) %>% data.frame()

methyl <- methyl[rowSums(!is.na(methyl)) > 0, ]

truncate_barcode <- function(x, len = 15) substr(x, 1, len)

dedupe_by_id <- function(mat_or_df, ids, margin = 2) {
  dup <- duplicated(ids)
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

rnaseq <- dedupe_by_id(rnaseq, rnaseq_ids_short, margin = 2)
methyl <- dedupe_by_id(methyl, methyl_ids_short, margin = 2)
pheno  <- dedupe_by_id(pheno, pheno_ids_short, margin = 1)

common_samples <- Reduce(intersect, list(colnames(rnaseq), colnames(methyl), rownames(pheno)))
common_samples <- sort(common_samples)

if (length(common_samples) == 0) {
  stop("No common samples found even after truncation. Check barcode formats manually.")
}

rnaseq <- rnaseq[, common_samples, drop = FALSE]
methyl <- methyl[, common_samples, drop = FALSE]
pheno  <- pheno[common_samples, , drop = FALSE]

stopifnot(all(colnames(rnaseq) == colnames(methyl)))
stopifnot(all(colnames(rnaseq) == rownames(pheno)))

cat("Final aligned dimensions:\n")
cat("rnaseq:", dim(rnaseq), "\n")
cat("methyl:", dim(methyl), "\n")
cat("pheno :", dim(pheno), "\n")

gene_vector <- c("MMP1", "PRAME", "HORMAD1", "NEK2", "KIF18B", "NMU")
gene_loc <- read_tsv("NCBI37.3.gene.loc", show_col_types = FALSE)

ann850k   <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata <- data.frame(ann850k[c("chr", "pos")])
annotdata$CpG <- rownames(annotdata)
annotdata$CHR <- gsub("chr", "", annotdata$chr)
annotdata <- annotdata %>% filter(CpG %in% rownames(methyl))

cov1 <- c("SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8")

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
  
  results <- limma::topTable(fit, coef = 2, n = Inf)
  results$SYMBOL <- rownames(results)
  return(results)
}

# =============================================================================
# MAIN FUNCTION: run full pipeline for a given flank size
# =============================================================================
run_eqtm_for_flank <- function(flank, flank_label) {
  
  cat("\n\n=====================================================\n")
  cat("Running eQTM analysis for flank =", flank, "bp (", flank_label, ")\n")
  cat("=====================================================\n\n")
  
  main_dir    <- paste0("Results_eQTM_Mval_", flank_label)
  result_dir  <- file.path(main_dir, "Results_per_CpG")
  plot_dir    <- file.path(main_dir, "Plots")
  table_dir   <- file.path(main_dir, "Tables")
  
  dir.create(main_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir,  showWarnings = FALSE, recursive = TRUE)
  
  # Save aligned datasets (same across flanks, but saved per-folder for completeness)
  saveRDS(pheno,  file.path(main_dir, "aligned_pheno.rds"))
  saveRDS(rnaseq, file.path(main_dir, "aligned_rnaseq.rds"))
  saveRDS(methyl, file.path(main_dir, "aligned_methyl.rds"))
  
  # -----------------------------------------------------------
  # Gene region definition for this flank
  # -----------------------------------------------------------
  gene_regions <- gene_loc %>%
    filter(NAME %in% gene_vector) %>%
    mutate(
      REGION_START = pmax(0, START - flank),
      REGION_END   = END + flank
    ) %>%
    select(NAME, CHR, START, END, REGION_START, REGION_END, STRAND)
  
  write.csv(gene_regions, file.path(table_dir, paste0("gene_regions_", flank_label, ".csv")), row.names = FALSE)
  
  # -----------------------------------------------------------
  # Map CpGs to gene regions
  # -----------------------------------------------------------
  cpg_gene_map <- lapply(seq_len(nrow(gene_regions)), function(i) {
    g <- gene_regions[i, ]
    annotdata %>%
      filter(CHR == g$CHR, pos >= g$REGION_START, pos <= g$REGION_END) %>%
      mutate(GENE = g$NAME)
  }) %>% bind_rows()
  
  write.csv(cpg_gene_map, file.path(table_dir, "cpg_gene_map.csv"), row.names = FALSE)
  
  cat("Number of CpGs mapped per gene (flank =", flank_label, "):\n")
  print(table(cpg_gene_map$GENE))
  
  if (nrow(cpg_gene_map) == 0) {
    cat("No CpGs mapped for this flank size. Skipping eQTM run.\n")
    return(NULL)
  }
  
  # -----------------------------------------------------------
  # Run eQTM analysis for all mapped CpGs
  # -----------------------------------------------------------
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
  
  cat("eQTM analysis complete for flank =", flank_label, "\n")
  
  # -----------------------------------------------------------
  # Combine all results
  # -----------------------------------------------------------
  result_files <- list.files(result_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(result_files) == 0) {
    cat("No result files generated for flank =", flank_label, ". Skipping downstream steps.\n")
    return(NULL)
  }
  
  all_results  <- lapply(result_files, readRDS) %>% bind_rows()
  
  saveRDS(all_results, file.path(main_dir, "all_eqtm_results_combined.rds"))
  write.csv(all_results, file.path(table_dir, "all_eqtm_results_combined.csv"), row.names = FALSE)
  
  # -----------------------------------------------------------
  # Filter to target genes
  # -----------------------------------------------------------
  final_results <- all_results %>%
    filter(Gene %in% gene_vector) %>%
    arrange(pvalue)
  
  write.csv(final_results, file.path(table_dir, "eqtm_results_target_genes.csv"), row.names = FALSE)
  saveRDS(final_results, file.path(main_dir, "eqtm_results_target_genes.rds"))
  
  cat("Top eQTM hits for target genes (flank =", flank_label, "):\n")
  print(head(final_results, 20))
  
  # -----------------------------------------------------------
  # Summary table
  # -----------------------------------------------------------
  summary_table <- final_results %>%
    group_by(Gene) %>%
    summarise(
      n_CpGs_tested = n(),
      n_sig_p05     = sum(pvalue < 0.05, na.rm = TRUE),
      n_sig_fdr05   = sum(adj.pvalue < 0.05, na.rm = TRUE),
      n_sig_p05_cor = sum(pvalue.cor < 0.05, na.rm = TRUE)
    ) %>%
    mutate(flank = flank_label)
  
  write.csv(summary_table, file.path(table_dir, "eqtm_summary_by_gene.csv"), row.names = FALSE)
  print(summary_table)
  
  # -----------------------------------------------------------
  # Plots
  # -----------------------------------------------------------
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
        title = paste("eQTM results for", g, "(flank =", flank_label, ")"),
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
  
  return(summary_table)
}

# =============================================================================
# Run for all flank sizes: 0kb, 1kb, 10kb, 100kb, 250kb
# =============================================================================
flank_settings <- list(
  "0kb"   = 0,
  "1kb"   = 1000,
  "10kb"  = 10000,
  "100kb" = 100000,
  "250kb" = 250000
)

all_summaries <- list()

for (flank_label in names(flank_settings)) {
  flank_val <- flank_settings[[flank_label]]
  
  summary_result <- tryCatch(
    run_eqtm_for_flank(flank_val, flank_label),
    error = function(e) {
      message("Error running flank ", flank_label, ": ", e$message)
      return(NULL)
    }
  )
  
  all_summaries[[flank_label]] <- summary_result
}

# =============================================================================
# Combine summary across all flank sizes for comparison
# =============================================================================
combined_summary <- bind_rows(all_summaries)

write.csv(combined_summary, "eqtm_summary_all_flanks.csv", row.names = FALSE)

cat("\n\n=====================================================\n")
cat("ALL FLANK SIZES COMPLETE\n")
cat("=====================================================\n")
print(combined_summary)