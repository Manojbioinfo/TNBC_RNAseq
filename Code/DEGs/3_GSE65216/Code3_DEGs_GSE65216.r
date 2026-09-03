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


dds<- readRDS("GSE65216_dds_clean.rds")
sample_info_pre=colData(dds)%>%as.data.frame()
vst_mat <- vst(dds, blind = TRUE)
expr_vst <- assay(vst_mat)

# PCA (top 5000 variable genes)
gene_var <- apply(expr_vst, 1, var)
top_probes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, length(gene_var))]
expr_top <- expr_vst[top_probes, ]
pca <- prcomp(t(expr_top), scale. = TRUE)
var_exp <- round(summary(pca)$importance[2, 1:5] * 100, 1)
pca_df <- as.data.frame(pca$x[,1:5])
pca_df$sample_id <- rownames(pca_df)
pca_df <- left_join(pca_df,
                    sample_info_pre %>% dplyr::select(sample_id, condition),
                    by = "sample_id")

GEO_ID="GSE65216"
# PCA plot colored by condition
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = condition, label = sample_id)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(size = 2.5, max.overlaps = 25) +
  scale_colour_manual(values = c("TNBC" = "darkred", "Healthy" = "navy")) +
  labs(title = paste0(GEO_ID, " PCA — Pre-QC"),
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Condition") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p_pca
ggsave(file.path(paste0("1_", GEO_ID, "_PCA_condition.png")),
       p_pca, width = 10, height = 7, dpi = 300)



###########################

# ============================================================
# Flag TNBC samples based on PC2 threshold
# Rule: TNBC with PC2 > 20
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

# -------------------------------
# 1. Define threshold (PC1)
#    Negative side = suspicious Healthy samples
# -------------------------------
pc1_threshold <- 10

# -------------------------------
# 2. Flag Healthy samples in negative PC1 region
# -------------------------------
pca_df <- pca_df %>%
  mutate(
    outlier_flag = condition == "Healthy" & PC1 < pc1_threshold
  )

# Check number of flagged samples
table(pca_df$outlier_flag)

# -------------------------------
# 3. PCA plot
# -------------------------------
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = condition)) +
  
  # All samples
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Highlight flagged Healthy samples
  geom_point(
    data = subset(pca_df, outlier_flag),
    shape = 21,
    size = 4,
    stroke = 1.2,
    color = "black"
  ) +
  
  # Label only flagged samples
  geom_text_repel(
    data = subset(pca_df, outlier_flag),
    aes(label = sample_id),
    size = 3,
    color = "black",
    max.overlaps = Inf
  ) +
  
  # Add vertical threshold line (PC1 cutoff)
  geom_vline(
    xintercept = pc1_threshold,
    linetype = "dashed",
    color = "black",
    size = 0.8
  ) +
  
  scale_colour_manual(values = c("TNBC" = "darkred", "Healthy" = "navy")) +
  
  labs(
    title = paste0(GEO_ID, " PCA — Healthy samples in negative PC1"),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    colour = "Condition"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Show plot
p_pca

# -------------------------------
# 4. Save plot
# -------------------------------
ggsave(
  paste0("1_", GEO_ID, "_PCA_PC1_negative_healthy.png"),
  p_pca,
  width = 10,
  height = 7,
  dpi = 300
)

saveRDS(p_pca, "GSE65216.PCA.rds")

# -------------------------------
# 5. Export flagged samples
# -------------------------------
flagged_samples <- pca_df %>%
  filter(outlier_flag) %>%
  select(sample_id, PC1, PC2)

flagged_samples

saveRDS(flagged_samples, "GSE65216.flagged_samples.rds")
