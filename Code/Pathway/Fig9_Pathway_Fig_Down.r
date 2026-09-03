# ============================================================================
# LOAD LIBRARIES AND DATA
# ============================================================================
# Clear workspace
rm(list = ls())

# Core packages
library(DESeq2)
library(sva)
library(dplyr)
library(tibble)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(cowplot)
library(jpeg)
#library(enrichR)
library(ggsci)
library(stringr)  # Added for str_extract

# Load data
final_pathway_data <- read.xlsx("Table2.all_pathways.xlsx")
remove=c("BioCarta_2016", "Elsevier_Pathway_Collection",  "WikiPathway_2023_Human") 

final_pathway_data=dplyr::filter(final_pathway_data,!Database%in%remove)


downregulated_pathways <- final_pathway_data %>%
  filter(grepl("Down", Dataset, fixed = TRUE))



cat("\n\n=== DOWNREGULATED PATHWAYS (TCGA_Down) ===\n")
print(downregulated_pathways)



final_pathway_data=downregulated_pathways

final_pathway_data$Term

# Prepare data for plotting
plot_data <- final_pathway_data %>%
  mutate(
    Gene_count = as.numeric(sub("/.*", "", Overlap)),  # Extract gene count
    Term = gsub("\\s*\\([^\\)]+\\)", "", as.character(Term)),  # Remove IDs
    Condition = paste(Dataset, Database, sep = "_")  # Unique key
  )

# Comprehensive cleaning of pathway terms
plot_data_top10 <- plot_data %>%
  mutate(
    # Remove all "Homo sapiens h ..." patterns
    #Term = str_remove(Term, "Homo sapiens h\\s*[A-Za-z0-9]+Pathway"),
    Term = str_remove(Term, "Homo sapiens h\\s*[A-Za-z0-9]+"),
    
    # Remove all WP identifiers (WP followed by numbers)
    Term = str_remove(Term, "\\s*WP\\d+$"),
    
    # Remove any remaining species identifiers
    Term = str_remove(Term, "Homo sapiens\\s*"),
    
    # Remove any trailing numbers or codes in parentheses
    Term = str_remove(Term, "\\s*\\(.*\\)$"),
    
    # Clean up extra spaces
    Term = str_trim(Term),
    
    # Remove any remaining trailing special characters
    Term = str_remove(Term, "[\\s\\.\\,]+$")
  ) 

plot_data_top10$Term

# Plot
p <- ggplot(plot_data_top10, aes(x = Condition, y = Term, color = FDR, size = Gene_count)) +
  geom_point() +
  scale_color_gradient(low = "#BC3C29FF", high = "#0072B5FF", guide = guide_colorbar(title.position = "top")) +
  scale_size_continuous(range = c(2, 6)) +
  theme_bw() +
  xlab("Cell Type with DEGs Direction") +
  ylab("Gene Ontology Term") +
  ggtitle("Top 10 GO Terms per Cell Type and Direction") +
  theme(
    text = element_text(size = 7, colour = "#45403f", face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 1, colour = "#45403f"),
    axis.text.y = element_text(colour = "#45403f"),
    legend.position = "top",
    legend.key.width = unit(1.5, "cm")
  )

p
# Save plots
ggsave("Fig2.Top10_Pathways_GO_down.pdf", plot = p, device = "pdf", width = 5, height = 7)

# Save JPEG version (300 DPI)
ggsave("Fig2.Top10_Pathways_GO_down.jpg", plot = p,
       device = "jpeg", width = 5, height = 7, units = "in",
       dpi = 300, quality = 90)

