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


dds<- readRDS("TCGA_dds_clean.rds")
sample_info_pre=colData(dds)%>%as.data.frame()

sample_info_pre


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
                    sample_info_pre %>% dplyr::select(sample_id, condition, gender, age, ethnicity),
                    by = "sample_id")

GEO_ID="TCGA"
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
# 1. Define threshold
# -------------------------------
pc2_threshold <- 10

# -------------------------------
# 2. Flag samples
# -------------------------------
pca_df <- pca_df %>%
  mutate(
    outlier_flag = condition == "TNBC" & PC2 > pc2_threshold
  )

# Check number of flagged samples
table(pca_df$outlier_flag)

# -------------------------------
# 3. PCA plot
# -------------------------------
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = condition)) +
  
  # All samples
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Highlight flagged TNBC samples
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
  
  # Add threshold line (VERY useful for paper)
  geom_hline(
    yintercept = pc2_threshold,
    linetype = "dashed",
    color = "black",
    size = 0.8
  ) +
  
  scale_colour_manual(values = c("TNBC" = "darkred", "Healthy" = "navy")) +
  
  labs(
    title = paste0(GEO_ID, " PCA — PC2 threshold filtering"),
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
  paste0("1_", GEO_ID, "_PCA_PC2_threshold.png"),
  p_pca,
  width = 10,
  height = 7,
  dpi = 300
)

saveRDS(p_pca,"TCGA.PCA.rds")
# -------------------------------
# 5. Export flagged samples
# -------------------------------
flagged_samples <- pca_df %>%
  filter(outlier_flag) %>%
  select(sample_id, PC1, PC2)

flagged_samples

saveRDS(flagged_samples,"TCGA.flagged_samples.rds")





##########


p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = ethnicity)) +
  geom_point(size = 3.5, alpha = 0.85) +
  labs(
    title = paste0(GEO_ID, " PCA — Pre-QC"),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    colour = "Ethnicity"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_pca

saveRDS(p_pca,"TCGA.PCA.ethnicity.rds")




p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = age)) +
  geom_point(size = 3.5, alpha = 0.85) +
  labs(
    title = paste0(GEO_ID),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    colour = "Age"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_pca

saveRDS(p_pca, "TCGA.PCA.age.rds")