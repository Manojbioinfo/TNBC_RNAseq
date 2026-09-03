## ============================================================
## FOREST PLOTS + SUMMARY TABLE FOR TRUE TNBC-ONLY COX RESULTS
## PUBLICATION-STYLE: forest panel + aligned numeric table panel
## ============================================================
rm(list = ls())

library(dplyr)
library(ggplot2)
library(openxlsx)
library(patchwork)   # for combining forest + table panels side-by-side

out_dir <- "Survival_Analysis_Results_CHECKED"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Load KM summary results
## ------------------------------------------------------------
summary_df <- read.csv("Survival_Analysis_TRUE_TNBC/KM_Summary_TRUE_TNBC.csv",
                       stringsAsFactors = FALSE)

cat("Columns available BEFORE rename:\n")
print(colnames(summary_df))

## ------------------------------------------------------------
## 2. Rename columns explicitly (avoid namespace clashes)
## ------------------------------------------------------------
summary_df <- dplyr::rename(summary_df,
                            HR_HighvsLow = Cox_HR,
                            CI_lower     = Cox_CI_lower,
                            CI_upper     = Cox_CI_upper,
                            p_value_cat  = Cox_p
)

stopifnot("HR_HighvsLow" %in% colnames(summary_df))
stopifnot("CI_lower"     %in% colnames(summary_df))
stopifnot("CI_upper"     %in% colnames(summary_df))
stopifnot("p_value_cat"  %in% colnames(summary_df))
cat("\nRename successful - proceeding.\n")

## ------------------------------------------------------------
## 3. Optional: merge in continuous HR/p if available
## ------------------------------------------------------------
cox_cont_path <- "Survival_Analysis_TRUE_TNBC/Cox_Survival_Results_TRUE_TNBC_ONLY.csv"

if (file.exists(cox_cont_path)) {
  cox_cont <- read.csv(cox_cont_path, stringsAsFactors = FALSE) %>%
    dplyr::select(Gene, HR_continuous, p_value_cont)
  summary_df <- dplyr::left_join(summary_df, cox_cont, by = "Gene")
} else {
  cat("\nNOTE: continuous Cox results file not found - HR_continuous/p_value_cont set to NA.\n")
  summary_df$HR_continuous <- NA_real_
  summary_df$p_value_cont  <- NA_real_
}

## ------------------------------------------------------------
## 4. Add derived columns
## ------------------------------------------------------------
summary_df <- summary_df %>%
  dplyr::mutate(
    Direction = ifelse(HR_HighvsLow > 1,
                       "High expression = WORSE survival (risk)",
                       "High expression = BETTER survival (protective)"),
    Significant_cat     = ifelse(p_value_cat < 0.05, "Yes", "No"),
    Significant_cont    = ifelse(!is.na(p_value_cont) & p_value_cont < 0.05, "Yes", "No"),
    Significant_FDR_cat = ifelse(FDR_Cox < 0.05, "Yes", "No"),
    log2HR_cat  = round(log2(HR_HighvsLow), 3),
    log2HR_cont = ifelse(!is.na(HR_continuous), round(log2(HR_continuous), 3), NA_real_),
    HR_label    = paste0(sprintf("%.2f", HR_HighvsLow),
                         " (", sprintf("%.2f", CI_lower),
                         "-", sprintf("%.2f", CI_upper), ")"),
    p_label     = ifelse(p_value_cat < 0.001, "p<0.001",
                         paste0("p=", sprintf("%.3f", p_value_cat))),
    p_display   = ifelse(p_value_cat < 0.001, "<0.001", sprintf("%.3f", p_value_cat)),
    FDR_display = ifelse(FDR_Cox < 0.001, "<0.001", sprintf("%.3f", FDR_Cox)),
    N_display   = paste0(N_total, " (", N_High, "H/", N_Low, "L)")
  ) %>%
  dplyr::arrange(HR_HighvsLow)

summary_df$Gene <- factor(summary_df$Gene, levels = summary_df$Gene)

cat("\nFinal summary_df structure:\n")
str(summary_df)

## ============================================================
## 5. PUBLICATION-STYLE COMBINED FOREST PLOT
##    Left: forest panel | Right: aligned numeric annotation table
## ============================================================

x_min <- min(summary_df$CI_lower, na.rm = TRUE) * 0.7
x_max <- max(summary_df$CI_upper, na.rm = TRUE) * 1.4

## ---- Left panel: forest plot ----
forest_panel <- ggplot(summary_df, aes(x = HR_HighvsLow, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Direction),
                 height = 0.22, linewidth = 1) +
  geom_point(aes(color = Direction, size = -log10(p_value_cat + 1e-6))) +
  scale_color_manual(values = c(
    "High expression = WORSE survival (risk)"        = "#D62728",
    "High expression = BETTER survival (protective)" = "#1F77B4"
  )) +
  scale_size_continuous(range = c(3, 6), guide = "none") +
  scale_x_log10(limits = c(x_min, x_max),
                breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  labs(
    title = "Cox Proportional Hazards: High vs Low Gene Expression",
    subtitle = paste0("Confirmed TNBC-only cohort (n = ", summary_df$N_total[1], ")  |  ",
                      "Point size \u221d statistical significance"),
    x = "Hazard Ratio (log scale)", y = NULL, color = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5, color = "grey30"),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold")
  )

## ---- Right panel: numeric annotation "table" as a ggplot ----
## Build a long-format df for aligned columns: N, HR (CI), p, FDR, Direction
tbl_df <- summary_df %>%
  dplyr::mutate(
    col_N   = N_display,
    col_HR  = HR_label,
    col_p   = p_display,
    col_FDR = FDR_display,
    col_dir = ifelse(HR_HighvsLow > 1, "Risk", "Protective")
  )

## reshape to long so each column of text sits at fixed x position
annot_long <- tbl_df %>%
  dplyr::select(Gene, col_N, col_HR, col_p, col_FDR, col_dir) %>%
  tidyr::pivot_longer(-Gene, names_to = "field", values_to = "value") %>%
  dplyr::mutate(
    x_pos = dplyr::case_when(
      field == "col_N"   ~ 1,
      field == "col_HR"  ~ 2,
      field == "col_p"   ~ 3,
      field == "col_FDR" ~ 4,
      field == "col_dir" ~ 5
    )
  )

header_labels <- data.frame(
  x_pos = 1:5,
  label = c("N (H/L)", "HR (95% CI)", "p-value", "FDR", "Direction")
)

table_panel <- ggplot(annot_long, aes(x = x_pos, y = Gene, label = value)) +
  geom_text(size = 3.6, hjust = 0.5) +
  geom_text(data = header_labels,
            aes(x = x_pos, y = length(levels(summary_df$Gene)) + 0.9, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 3.8) +
  scale_x_continuous(limits = c(0.5, 5.5), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.1))) +
  theme_void(base_size = 13) +
  theme(
    plot.margin = margin(t = 40, r = 10, b = 5, l = 10)
  )

## ---- Combine with patchwork ----
combined_forest <- forest_panel + table_panel +
  patchwork::plot_layout(widths = c(1.3, 1.6))

ggsave(file.path(out_dir, "Forest_Plot_COMBINED_6genes_TRUE_TNBC.pdf"),
       plot = combined_forest, width = 14, height = 6.5)
ggsave(file.path(out_dir, "Forest_Plot_COMBINED_6genes_TRUE_TNBC.png"),
       plot = combined_forest, width = 14, height = 6.5, dpi = 320)

cat("Publication-style combined forest plot (with table panel) saved.\n")

## ============================================================
## 6. INDIVIDUAL forest plots (enriched with stats box)
## ============================================================
for (g in levels(summary_df$Gene)) {
  
  row <- summary_df %>% dplyr::filter(Gene == g)
  col_main <- ifelse(row$HR_HighvsLow > 1, "#D62728", "#1F77B4")
  
  p_ind <- ggplot(row, aes(x = HR_HighvsLow, y = 1)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.12,
                   color = col_main, linewidth = 1.3) +
    geom_point(size = 5, color = col_main) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    annotate("label", x = row$CI_upper, y = 1.35,
             label = paste0(
               "HR = ", row$HR_label, "\n",
               "p = ", row$p_display, "  |  FDR = ", row$FDR_display, "\n",
               "N = ", row$N_total, " (High=", row$N_High, ", Low=", row$N_Low, ")"
             ),
             hjust = 1, size = 3.4, fill = "white", label.size = 0.3) +
    ylim(0.6, 1.6) +
    labs(
      title = paste0(g, " — Confirmed TNBC-only cohort"),
      subtitle = row$Direction,
      x = "Hazard Ratio (log scale)", y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10.5, color = col_main, face = "italic"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(out_dir, paste0("Forest_Plot_", g, "_TRUE_TNBC.pdf")),
         plot = p_ind, width = 6, height = 4)
  ggsave(file.path(out_dir, paste0("Forest_Plot_", g, "_TRUE_TNBC.png")),
         plot = p_ind, width = 6, height = 4, dpi = 320)
}

cat("Individual enriched forest plots saved for:",
    paste(levels(summary_df$Gene), collapse = ", "), "\n")

## ------------------------------------------------------------
## 7. Save master table (CSV + styled XLSX)
## ------------------------------------------------------------
write.csv(summary_df,
          file.path(out_dir, "Cox_Survival_MASTER_TABLE_TRUE_TNBC.csv"),
          row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Master_Table")
writeData(wb, "Master_Table", summary_df)

headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
addStyle(wb, "Master_Table", headerStyle, rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
setColWidths(wb, "Master_Table", cols = 1:ncol(summary_df), widths = "auto")
freezePane(wb, "Master_Table", firstRow = TRUE)

sigStyle <- createStyle(fgFill = "#FFF2CC")
sig_rows <- which(summary_df$Significant_cat == "Yes") + 1
if (length(sig_rows) > 0) {
  addStyle(wb, "Master_Table", sigStyle, rows = sig_rows,
           cols = 1:ncol(summary_df), gridExpand = TRUE)
}

saveWorkbook(wb, file.path(out_dir, "Cox_Survival_MASTER_TABLE_TRUE_TNBC.xlsx"), overwrite = TRUE)
cat("Master table saved as CSV and XLSX.\n")

## ------------------------------------------------------------
## 8. Final console report
## ------------------------------------------------------------
cat("\n=== FINAL SUMMARY TABLE (TRUE TNBC ONLY, n=", summary_df$N_total[1], ") ===\n")
print(summary_df %>%
        dplyr::select(Gene, HR_HighvsLow, CI_lower, CI_upper, p_value_cat,
                      HR_continuous, p_value_cont, FDR_Cox, Direction))

cat("\nAll outputs saved in:\n", out_dir, "\n")
cat("Files generated:\n")
cat(" - Cox_Survival_MASTER_TABLE_TRUE_TNBC.csv/.xlsx\n")
cat(" - Forest_Plot_COMBINED_6genes_TRUE_TNBC.pdf/.png  (forest + numeric table panel)\n")
cat(" - Forest_Plot_<GENE>_TRUE_TNBC.pdf/.png (x6, with stats annotation box)\n")