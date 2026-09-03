## ============================================================
## SIMPLE RBIND OF PER-FEATURE RDS FILES - ONE BLOCK PER OMICS LAYER
## ============================================================

base_dir <- "Z:/LAPTOP/Newmanuscript/Paper/2_Working/W2_Gopi_TNBC/Multiomic_Aug26_Complete_Data"
out_dir  <- file.path(base_dir, "Result_Combined")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



## ------------------------------------------------------------
## Blood
## ------------------------------------------------------------
in_dir <- file.path(base_dir, "Result_Blood", "Results_per_Blood")
files  <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Blood files:", length(files), "\n")

list_Blood <- lapply(files, readRDS)
combined_Blood <- do.call(rbind, list_Blood)



## ------------------------------------------------------------
## Blood - keep Estimate, StdError, z-score, pval (+ ID columns)
## ------------------------------------------------------------

## Keep only needed columns
combined_Blood <- combined_Blood[, c("Gene", "Protein", "Estimate", "StdError", "z-score", "pval")]

## Recompute global FDR across all rows
combined_Blood$FDR <- p.adjust(combined_Blood$pval, method = "fdr")
colnames(combined_Blood )[2]="Celltype"

head(combined_Blood)

## Filter significant results
combined_Blood_sig <- combined_Blood[combined_Blood$FDR < 0.05, ]

cat("Total rows:", nrow(combined_Blood), "\n")
cat("Significant (FDR<0.05):", nrow(combined_Blood_sig), "\n")

head(combined_Blood_sig)

## Save
saveRDS(combined_Blood,     file.path(out_dir, "Combined_Blood_Results_clean.rds"))
saveRDS(combined_Blood_sig, file.path(out_dir, "Combined_Blood_Results_sig.rds"))
## ------------------------------------------------------------
## Proteomic
## ------------------------------------------------------------
in_dir <- file.path(base_dir, "Result_proteomic", "Results_per_Protein")
files  <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Proteomic files:", length(files), "\n")

list_Proteomic <- lapply(files, readRDS)
combined_Proteomic <- do.call(rbind, list_Proteomic)

head(combined_Proteomic )

## ------------------------------------------------------------
## Proteomic - keep Estimate, StdError, z-score, pval (+ ID columns)
## ------------------------------------------------------------

## Keep only needed columns
combined_Proteomic <- combined_Proteomic[, c("Gene", "Protein", "Estimate", "StdError", "z-score", "pval")]

## Recompute global FDR across all rows
combined_Proteomic$FDR <- p.adjust(combined_Proteomic$pval, method = "fdr")

head(combined_Proteomic)

## Filter significant results
combined_Proteomic_sig <- combined_Proteomic[combined_Proteomic$FDR < 0.05, ]

cat("Total rows:", nrow(combined_Proteomic), "\n")
cat("Significant (FDR<0.05):", nrow(combined_Proteomic_sig), "\n")

head(combined_Proteomic_sig)

## Save
saveRDS(combined_Proteomic,     file.path(out_dir, "Combined_Proteomic_Results_clean.rds"))
saveRDS(combined_Proteomic_sig, file.path(out_dir, "Combined_Proteomic_Results_sig.rds"))

## ------------------------------------------------------------
## Viral
## ------------------------------------------------------------
in_dir <- file.path(base_dir, "Result_viral", "Results_per_viral")
files  <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Viral files:", length(files), "\n")

list_Viral <- lapply(files, readRDS)
combined_Viral <- do.call(rbind, list_Viral)

## ------------------------------------------------------------
## Viral - keep Estimate, StdError, z-score, pval (+ ID columns)
## ------------------------------------------------------------

## Keep only needed columns
combined_Viral <- combined_Viral[, c("Gene", "Protein", "Estimate", "StdError", "z-score", "pval")]

## Rename Protein -> Virus for clarity (optional, since it's viral not proteomic)
colnames(combined_Viral)[colnames(combined_Viral) == "Protein"] <- "Virus"

## Recompute global FDR across all rows
combined_Viral$FDR <- p.adjust(combined_Viral$pval, method = "fdr")

head(combined_Viral)

## Filter significant results
combined_Viral_sig <- combined_Viral[combined_Viral$FDR < 0.05, ]

cat("Total rows:", nrow(combined_Viral), "\n")
cat("Significant (FDR<0.05):", nrow(combined_Viral_sig), "\n")

head(combined_Viral_sig)

## Save
saveRDS(combined_Viral,     file.path(out_dir, "Combined_Viral_Results_clean.rds"))
saveRDS(combined_Viral_sig, file.path(out_dir, "Combined_Viral_Results_sig.rds"))

