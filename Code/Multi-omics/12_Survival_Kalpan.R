## ============================================================
## KAPLAN-MEIER CURVES FOR CONFIRMED TRUE-TNBC COHORT
## (uses the exact aligned surv_mat / expr_mat_numeric / df logic
##  from the verified TNBC-restriction script)
## Produces: combined panel, individual plots, full summary table
## ============================================================
rm(list = ls())

library(survival)
library(survminer)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(gridExtra)
library(openxlsx)

## ------------------------------------------------------------
## 1. Load everything (identical to verified TNBC-restriction script)
## ------------------------------------------------------------
dds      <- readRDS("../Analysed_Dataset_May2026/1_TCGA/TCGA_dds_clean_with_SVs.rds")
surv_mat <- readRDS("RNAseq_Survival/survival_aligned.rds")
expr_mat <- readRDS("RNAseq_Survival/rnaseq_aligned.rds")

pheno <- colData(dds) %>% as.data.frame()

truncate_barcode <- function(x, len = 15) substr(x, 1, len)

pheno_ids_short <- truncate_barcode(rownames(pheno))
surv_ids_short  <- truncate_barcode(rownames(surv_mat))
expr_ids_short  <- truncate_barcode(colnames(expr_mat))

pheno    <- pheno[!duplicated(pheno_ids_short), , drop = FALSE]
rownames(pheno) <- pheno_ids_short[!duplicated(pheno_ids_short)]

surv_mat <- surv_mat[!duplicated(surv_ids_short), , drop = FALSE]
rownames(surv_mat) <- surv_ids_short[!duplicated(surv_ids_short)]

expr_mat <- expr_mat[, !duplicated(expr_ids_short), drop = FALSE]
colnames(expr_mat) <- expr_ids_short[!duplicated(expr_ids_short)]

## ------------------------------------------------------------
## 2. Restrict to TRUE TNBC + survival + expr overlap
## ------------------------------------------------------------
tnbc_ids <- rownames(pheno)[pheno$condition == "TNBC"]

common_samples <- Reduce(intersect, list(
  rownames(surv_mat),
  colnames(expr_mat),
  tnbc_ids
))
common_samples <- sort(common_samples)

if (length(common_samples) == 0) {
  stop("No common samples found. Check barcode formats / condition column values manually.")
}

surv_mat <- surv_mat[common_samples, , drop = FALSE]
expr_mat <- expr_mat[, common_samples, drop = FALSE]
stopifnot(all(rownames(surv_mat) == colnames(expr_mat)))

cat("Final TRUE-TNBC cohort size:", length(common_samples), "\n")

expr_mat_numeric <- as.matrix(expr_mat)
mode(expr_mat_numeric) <- "numeric"

genes <- c("HORMAD1", "NEK2", "MMP1", "KIF18B", "PRAME", "NMU")
genes <- genes[genes %in% rownames(expr_mat_numeric)]
cat("Genes available:", paste(genes, collapse = ", "), "\n")

## ------------------------------------------------------------
## 3. Output directory
## ------------------------------------------------------------
out_dir <- "Survival_Analysis_TRUE_TNBC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 4. Loop: build KM fit + plot + summary stats per gene
## ------------------------------------------------------------
km_plots     <- list()
km_summary   <- list()
cox_results  <- list()

for (g in genes) {
  
  cat("\n=====================================\n")
  cat("Processing gene:", g, "\n")
  
  expr_vals <- expr_mat_numeric[g, ]
  
  df <- data.frame(
    sample    = names(expr_vals),
    expr_cont = as.numeric(expr_vals),
    stringsAsFactors = FALSE
  )
  df <- merge(df, surv_mat, by.x = "sample", by.y = "row.names")
  df <- df[!is.na(df$OS) & !is.na(df$OS.time) & !is.na(df$expr_cont), ]
  
  if (nrow(df) < 10) {
    cat("Too few samples for", g, "- skipping.\n")
    next
  }
  
  med_val  <- median(df$expr_cont, na.rm = TRUE)
  df$group <- ifelse(df$expr_cont > med_val, "High", "Low")
  df$group <- factor(df$group, levels = c("Low", "High"))
  
  df$OS      <- as.numeric(df$OS)
  df$OS.time <- as.numeric(df$OS.time)
  
  ## ---- Cox model (categorical) for HR/CI/p annotation ----
  cox_cat <- tryCatch(
    coxph(Surv(OS.time, OS) ~ group, data = df),
    error = function(e) NULL
  )
  if (is.null(cox_cat)) {
    cat("Cox model failed for", g, "- skipping.\n")
    next
  }
  s_cat    <- summary(cox_cat)
  hr_cat   <- s_cat$conf.int["groupHigh", "exp(coef)"]
  ci_lower <- s_cat$conf.int["groupHigh", "lower .95"]
  ci_upper <- s_cat$conf.int["groupHigh", "upper .95"]
  p_cat    <- s_cat$coefficients["groupHigh", "Pr(>|z|)"]
  
  cox_results[[g]] <- data.frame(
    Gene = g, HR = hr_cat, CI_lower = ci_lower, CI_upper = ci_upper, p_value = p_cat
  )
  
  hr_label <- paste0("HR = ", sprintf("%.2f", hr_cat),
                     " (", sprintf("%.2f", ci_lower), "-", sprintf("%.2f", ci_upper), ")")
  
  ## ---- KM fit ----
  fit <- survfit(Surv(OS.time, OS) ~ group, data = df)
  
  n_low  <- sum(df$group == "Low")
  n_high <- sum(df$group == "High")
  
  ## ---- Individual KM plot (with risk table) ----
  p <- ggsurvplot(
    fit, data = df,
    pval          = TRUE,
    conf.int      = FALSE,
    risk.table    = TRUE,
    risk.table.height = 0.25,
    palette       = c("#1F77B4", "#D62728"),
    legend.title  = g,
    legend.labs   = c(paste0("Low (n=", n_low, ")"), paste0("High (n=", n_high, ")")),
    title         = paste0(g, " -- Confirmed TNBC only (n=", nrow(df), ")"),
    subtitle      = hr_label,
    xlab          = "Time (days)",
    ylab          = "Overall survival probability",
    ggtheme       = theme_minimal(base_size = 12)
  )
  
  pdf(file.path(out_dir, paste0("KM_TRUE_TNBC_", g, ".pdf")), width = 6, height = 6)
  print(p)
  dev.off()
  
  png(file.path(out_dir, paste0("KM_TRUE_TNBC_", g, ".png")), width = 6, height = 6, units = "in", res = 300)
  print(p)
  dev.off()
  
  ## store simplified ggplot (no risk table) for combined panel
  km_plots[[g]] <- p$plot +
    labs(subtitle = hr_label) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9))
  
  ## ---- log-rank test (independent cross-check) ----
  lr   <- survdiff(Surv(OS.time, OS) ~ group, data = df)
  lr_p <- 1 - pchisq(lr$chisq, df = 1)
  
  km_summary[[g]] <- data.frame(
    Gene           = g,
    N_total        = nrow(df),
    N_Low          = n_low,
    N_High         = n_high,
    Median_OS_Low  = summary(fit)$table["group=Low", "median"],
    Median_OS_High = summary(fit)$table["group=High", "median"],
    LogRank_p      = round(lr_p, 5),
    Cox_HR         = round(hr_cat, 3),
    Cox_CI_lower   = round(ci_lower, 3),
    Cox_CI_upper   = round(ci_upper, 3),
    Cox_p          = round(p_cat, 4)
  )
}

cat("\nIndividual KM plots saved for:", paste(names(km_plots), collapse = ", "), "\n")

## ------------------------------------------------------------
## 5. Combined KM panel (2x3 grid, no risk tables, for overview figure)
## ------------------------------------------------------------
combined_km <- gridExtra::grid.arrange(grobs = km_plots, ncol = 2, nrow = 3)

ggsave(file.path(out_dir, "KM_COMBINED_6genes_TRUE_TNBC.pdf"),
       plot = combined_km, width = 10, height = 12)
ggsave(file.path(out_dir, "KM_COMBINED_6genes_TRUE_TNBC.png"),
       plot = combined_km, width = 10, height = 12, dpi = 300)

cat("Combined KM panel saved (PDF + PNG).\n")

## ------------------------------------------------------------
## 6. Save summary table (CSV + XLSX) - full downstream-ready table
## ------------------------------------------------------------
km_summary_df <- do.call(rbind, km_summary)
rownames(km_summary_df) <- NULL

km_summary_df$FDR_LogRank <- p.adjust(km_summary_df$LogRank_p, method = "fdr")
km_summary_df$FDR_Cox     <- p.adjust(km_summary_df$Cox_p, method = "fdr")

write.csv(km_summary_df,
          file.path(out_dir, "KM_Summary_TRUE_TNBC.csv"),
          row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "KM_Summary_TRUE_TNBC")
writeData(wb, "KM_Summary_TRUE_TNBC", km_summary_df)

headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
addStyle(wb, "KM_Summary_TRUE_TNBC", headerStyle, rows = 1,
         cols = 1:ncol(km_summary_df), gridExpand = TRUE)
setColWidths(wb, "KM_Summary_TRUE_TNBC", cols = 1:ncol(km_summary_df), widths = "auto")
freezePane(wb, "KM_Summary_TRUE_TNBC", firstRow = TRUE)

saveWorkbook(wb, file.path(out_dir, "KM_Summary_TRUE_TNBC.xlsx"), overwrite = TRUE)

cat("\nKM summary table saved as CSV and XLSX.\n")

## ------------------------------------------------------------
## 7. Final console report
## ------------------------------------------------------------
cat("\n=== KAPLAN-MEIER SUMMARY (CONFIRMED TRUE TNBC ONLY) ===\n")
print(km_summary_df)

cat("\nAll outputs saved in:", out_dir, "\n")
cat("Files generated:\n")
cat(" - KM_TRUE_TNBC_<GENE>.pdf/.png (x", length(km_plots), ", individual, with risk tables)\n", sep="")
cat(" - KM_COMBINED_6genes_TRUE_TNBC.pdf/.png (2x3 panel)\n")
cat(" - KM_Summary_TRUE_TNBC.csv/.xlsx\n")