## =============================================================================
## Multi-Omic Directionality Report: Proteomic, Blood, and Viral Associations
## For candidate genes: NMU, NEK2, HORMAD1, PRAME, KIF18B, MMP1
## =============================================================================

library(dplyr)
library(tidyr)

dir.create("Result_Combined", showWarnings = FALSE)

## -----------------------------------------------------------------------
## 1. Load all combined datasets
## -----------------------------------------------------------------------
proteomic <- readRDS("Result_Combined/Combined_Proteomic_Results_clean.rds")
blood     <- readRDS("Result_Combined/Combined_Blood_Results_clean.rds")
viral     <- readRDS("Result_Combined/Combined_Viral_Results_clean.rds")

cat("Proteomic columns:", paste(colnames(proteomic), collapse = ", "), "\n")
cat("Blood columns:    ", paste(colnames(blood), collapse = ", "), "\n")
cat("Viral columns:    ", paste(colnames(viral), collapse = ", "), "\n")

## -----------------------------------------------------------------------
## 2. Standardize each dataset - rename feature column to "Feature"
## -----------------------------------------------------------------------
add_direction <- function(df, omic_name, feature_col) {
  df %>%
    rename(Feature = all_of(feature_col)) %>%
    mutate(
      Omic       = omic_name,
      direction  = case_when(
        Estimate > 0 ~ "Positive",
        Estimate < 0 ~ "Negative",
        TRUE ~ "No effect"
      ),
      sig_status = ifelse(FDR < 0.05, "Significant", "Not Significant")
    ) %>%
    select(Omic, Gene, Feature, Estimate, StdError, `z-score`, pval, FDR,
           direction, sig_status)
}

proteomic2 <- add_direction(proteomic, "Proteomic", "Protein")
blood2     <- add_direction(blood,     "Blood",     "Celltype")
viral2     <- add_direction(viral,     "Viral",     "Virus")

all_data <- bind_rows(proteomic2, blood2, viral2)

## -----------------------------------------------------------------------
## 3. Per-gene, per-omic directionality summary (ALL pairs)
## -----------------------------------------------------------------------
gene_summary_all <- all_data %>%
  group_by(Omic, Gene) %>%
  summarise(
    n_total        = n(),
    n_positive     = sum(direction == "Positive"),
    n_negative     = sum(direction == "Negative"),
    pct_positive   = round(100 * n_positive / n_total, 1),
    pct_negative   = round(100 * n_negative / n_total, 1),
    n_sig          = sum(sig_status == "Significant"),
    n_sig_positive = sum(sig_status == "Significant" & direction == "Positive"),
    n_sig_negative = sum(sig_status == "Significant" & direction == "Negative"),
    .groups = "drop"
  )

print(gene_summary_all)

## -----------------------------------------------------------------------
## 4. Overall directionality per omic (collapsed across genes)
## -----------------------------------------------------------------------
overall_by_omic <- all_data %>%
  group_by(Omic) %>%
  summarise(
    n_total      = n(),
    n_positive   = sum(direction == "Positive"),
    n_negative   = sum(direction == "Negative"),
    pct_positive = round(100 * n_positive / n_total, 1),
    pct_negative = round(100 * n_negative / n_total, 1),
    n_sig        = sum(sig_status == "Significant"),
    pct_sig      = round(100 * n_sig / n_total, 1),
    .groups = "drop"
  )

print(overall_by_omic)

## -----------------------------------------------------------------------
## 5. Significant-only directionality per omic
## -----------------------------------------------------------------------
sig_by_omic <- all_data %>%
  filter(sig_status == "Significant") %>%
  group_by(Omic) %>%
  summarise(
    n_sig_total      = n(),
    n_sig_positive   = sum(direction == "Positive"),
    n_sig_negative   = sum(direction == "Negative"),
    pct_sig_positive = round(100 * n_sig_positive / n_sig_total, 1),
    pct_sig_negative = round(100 * n_sig_negative / n_sig_total, 1),
    .groups = "drop"
  )

print(sig_by_omic)

## -----------------------------------------------------------------------
## 6. Significant hits detail
## -----------------------------------------------------------------------
sig_hits_detail <- all_data %>%
  filter(sig_status == "Significant") %>%
  select(Omic, Gene, Feature, Estimate, `z-score`, pval, FDR, direction) %>%
  arrange(Omic, Gene, FDR)

## -----------------------------------------------------------------------
## 7. Binomial test per omic: is direction randomly split (50/50)?
## -----------------------------------------------------------------------
binom_by_omic <- all_data %>%
  group_by(Omic) %>%
  summarise(
    n_total    = n(),
    n_positive = sum(direction == "Positive"),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(binom_pval = binom.test(n_positive, n_total, p = 0.5)$p.value) %>%
  ungroup()

print(binom_by_omic)

## -----------------------------------------------------------------------
## 8. Binomial test per gene WITHIN each omic (which genes drive direction)
## -----------------------------------------------------------------------
binom_by_gene_omic <- all_data %>%
  group_by(Omic, Gene) %>%
  summarise(
    n_total    = n(),
    n_positive = sum(direction == "Positive"),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(binom_pval = binom.test(n_positive, n_total, p = 0.5)$p.value) %>%
  ungroup() %>%
  mutate(sig_skew = ifelse(binom_pval < 0.05, "*", ""))

print(binom_by_gene_omic)

## -----------------------------------------------------------------------
## 9. Chi-square: is direction associated with Gene identity (within omic)?
## -----------------------------------------------------------------------
chisq_by_omic <- list()
for (omic_name in c("Proteomic", "Blood", "Viral")) {
  sub <- all_data %>% filter(Omic == omic_name, direction %in% c("Positive", "Negative"))
  tab <- table(sub$Gene, sub$direction)
  if (all(dim(tab) > 1)) {
    chisq_by_omic[[omic_name]] <- suppressWarnings(chisq.test(tab))
  }
}

## -----------------------------------------------------------------------
## 10. Save summary tables
## -----------------------------------------------------------------------
saveRDS(gene_summary_all, "Result_Combined/GeneSummary_AllOmics.rds")
write.csv(gene_summary_all, "Result_Combined/GeneSummary_AllOmics.csv", row.names = FALSE)
write.csv(sig_hits_detail, "Result_Combined/SignificantHits_AllOmics.csv", row.names = FALSE)
write.csv(binom_by_gene_omic, "Result_Combined/BinomialTest_ByGeneOmic.csv", row.names = FALSE)

## =============================================================================
## 11. WRITE FULL TEXT REPORT
## =============================================================================
report_path <- "Result_Combined/MultiOmic_Directionality_Report.txt"
sink(report_path)

cat("=================================================================\n")
cat("MULTI-OMIC DIRECTIONALITY REPORT\n")
cat("Candidate Genes: NMU, NEK2, HORMAD1, PRAME, KIF18B, MMP1\n")
cat("Omics: Proteomic (protein arrays), Blood (immune cell types),\n")
cat("       Viral (viral seropositivity/load)\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=================================================================\n\n")

cat("-----------------------------------------------------------------\n")
cat("1. OVERALL DIRECTIONALITY BY OMIC LAYER (ALL PAIRS TESTED)\n")
cat("-----------------------------------------------------------------\n")
for (i in seq_len(nrow(overall_by_omic))) {
  r <- overall_by_omic[i, ]
  cat(sprintf("\n%s (n = %d tests):\n", r$Omic, r$n_total))
  cat(sprintf("  Positive association: %d (%.1f%%)\n", r$n_positive, r$pct_positive))
  cat(sprintf("  Negative association: %d (%.1f%%)\n", r$n_negative, r$pct_negative))
  cat(sprintf("  Significant (FDR<0.05): %d (%.1f%%)\n", r$n_sig, r$pct_sig))
}

cat("\n-----------------------------------------------------------------\n")
cat("2. DIRECTIONALITY AMONG SIGNIFICANT HITS ONLY (FDR < 0.05)\n")
cat("-----------------------------------------------------------------\n")
if (nrow(sig_by_omic) > 0) {
  for (i in seq_len(nrow(sig_by_omic))) {
    r <- sig_by_omic[i, ]
    cat(sprintf("\n%s (n significant = %d):\n", r$Omic, r$n_sig_total))
    cat(sprintf("  Positive: %d (%.1f%%)\n", r$n_sig_positive, r$pct_sig_positive))
    cat(sprintf("  Negative: %d (%.1f%%)\n", r$n_sig_negative, r$pct_sig_negative))
  }
} else {
  cat("\nNo significant hits found at FDR < 0.05 in any omic layer.\n")
}

cat("\n-----------------------------------------------------------------\n")
cat("3. BINOMIAL TEST (OVERALL PER OMIC): IS DIRECTION SKEWED FROM 50/50?\n")
cat("-----------------------------------------------------------------\n")
for (i in seq_len(nrow(binom_by_omic))) {
  r <- binom_by_omic[i, ]
  cat(sprintf("\n%s: %d/%d positive (%.1f%%), binomial p = %.3e %s\n",
              r$Omic, r$n_positive, r$n_total,
              100 * r$n_positive / r$n_total,
              r$binom_pval,
              ifelse(r$binom_pval < 0.05, "*** SIGNIFICANT SKEW ***", "(not significant)")))
}

cat("\n-----------------------------------------------------------------\n")
cat("4. BINOMIAL TEST PER GENE WITHIN EACH OMIC\n")
cat("   (identifies which genes show a consistent +/- direction)\n")
cat("-----------------------------------------------------------------\n")
for (omic_name in c("Proteomic", "Blood", "Viral")) {
  cat(sprintf("\n### %s ###\n", omic_name))
  sub <- binom_by_gene_omic %>% filter(Omic == omic_name) %>% arrange(binom_pval)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("  %-10s: %d/%d positive (%.1f%%), p = %.3e %s\n",
                r$Gene, r$n_positive, r$n_total,
                100 * r$n_positive / r$n_total, r$binom_pval, r$sig_skew))
  }
}

cat("\n-----------------------------------------------------------------\n")
cat("5. CHI-SQUARE: IS DIRECTION ASSOCIATED WITH GENE IDENTITY?\n")
cat("   (tests whether some genes are more consistently + or - than others)\n")
cat("-----------------------------------------------------------------\n")
for (omic_name in names(chisq_by_omic)) {
  ct <- chisq_by_omic[[omic_name]]
  cat(sprintf("\n%s: X-squared = %.2f, df = %d, p-value = %.3e %s\n",
              omic_name, ct$statistic, ct$parameter, ct$p.value,
              ifelse(ct$p.value < 0.05, "*** SIGNIFICANT ***", "(not significant)")))
}

cat("\n-----------------------------------------------------------------\n")
cat("6. PER-GENE DIRECTIONALITY SUMMARY (ALL PAIRS)\n")
cat("-----------------------------------------------------------------\n")
for (omic_name in c("Proteomic", "Blood", "Viral")) {
  cat(sprintf("\n### %s ###\n", omic_name))
  sub <- gene_summary_all %>% filter(Omic == omic_name) %>% arrange(desc(pct_negative))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("  %-10s: %d tests | +:%d (%.1f%%)  -:%d (%.1f%%) | Sig: %d (+:%d / -:%d)\n",
                r$Gene, r$n_total, r$n_positive, r$pct_positive,
                r$n_negative, r$pct_negative,
                r$n_sig, r$n_sig_positive, r$n_sig_negative))
  }
}

cat("\n-----------------------------------------------------------------\n")
cat("7. SIGNIFICANT HITS DETAIL (FDR < 0.05)\n")
cat("-----------------------------------------------------------------\n")
if (nrow(sig_hits_detail) > 0) {
  for (omic_name in c("Proteomic", "Blood", "Viral")) {
    cat(sprintf("\n### %s ###\n", omic_name))
    sub <- sig_hits_detail %>% filter(Omic == omic_name)
    if (nrow(sub) == 0) {
      cat("  No significant hits.\n")
      next
    }
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      cat(sprintf("  %-10s <-> %-25s : Estimate = %+.4e, FDR = %.3e, Direction = %s\n",
                  r$Gene, r$Feature, r$Estimate, r$FDR, r$direction))
    }
  }
} else {
  cat("\nNo significant hits found across any omic layer.\n")
}

cat("\n-----------------------------------------------------------------\n")
cat("8. BIOLOGICAL INTERPRETATION SUMMARY\n")
cat("-----------------------------------------------------------------\n")
for (i in seq_len(nrow(overall_by_omic))) {
  r <- overall_by_omic[i, ]
  dominant <- ifelse(r$pct_negative > r$pct_positive, "predominantly NEGATIVE",
                     ifelse(r$pct_positive > r$pct_negative, "predominantly POSITIVE", "balanced"))
  cat(sprintf("- %s associations are %s (%.1f%% negative vs %.1f%% positive across %d tests).\n",
              r$Omic, dominant, r$pct_negative, r$pct_positive, r$n_total))
}

cat("\nPer-gene highlights:\n")
for (g in unique(gene_summary_all$Gene)) {
  cat(sprintf("\n  %s:\n", g))
  sub <- gene_summary_all %>% filter(Gene == g)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("    - %s: %.1f%% negative, %.1f%% positive (%d sig hits)\n",
                r$Omic, r$pct_negative, r$pct_positive, r$n_sig))
  }
}

cat("\n=================================================================\n")
cat("END OF REPORT\n")
cat("=================================================================\n")

sink()

cat("\nReport successfully written to:", report_path, "\n")
cat("Supporting tables saved:\n")
cat("  - Result_Combined/GeneSummary_AllOmics.csv\n")
cat("  - Result_Combined/SignificantHits_AllOmics.csv\n")
cat("  - Result_Combined/BinomialTest_ByGeneOmic.csv\n")



## Which genes have the most FDR-significant hits, and in which omic?
sig_summary <- all_data %>%
  filter(FDR < 0.05) %>%
  group_by(Gene, Omic) %>%
  summarise(n_sig = n(), .groups = "drop") %>%
  pivot_wider(names_from = Omic, values_from = n_sig, values_fill = 0) %>%
  arrange(desc(rowSums(across(where(is.numeric)))))

print(sig_summary)

## Effect size magnitude comparison across omics (ignoring sign)
magnitude_summary <- all_data %>%
  mutate(abs_estimate = abs(Estimate)) %>%
  group_by(Omic) %>%
  summarise(
    median_abs_effect = median(abs_estimate),
    mean_abs_effect   = mean(abs_estimate),
    n_sig             = sum(FDR < 0.05),
    pct_sig           = 100 * n_sig / n(),
    .groups = "drop"
  )

print(magnitude_summary)