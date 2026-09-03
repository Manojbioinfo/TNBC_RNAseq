# ============================================================================
# LOAD LIBRARIES AND DATA
# ============================================================================
# Clear workspace
rm(list = ls())
Sys.setenv(http_proxy = "http://proxy.mh-hannover.de:8080")
Sys.setenv(https_proxy = "http://proxy.mh-hannover.de:8080")
options(download.file.method = "curl")
options(download.file.extra = "-L --proxy http://proxy.mh-hannover.de:8080")
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
library(tidyverse)
library(enrichR)
library(ggsci)
library(stringr)  # Added for str_extract

# Load data
data1 <- read.csv("../../Analysed_Dataset_May2026/1_TCGA/DEGs_TNBC_vs_Healthy_with_SV_significant.csv")
data2 <- read.csv("../../Analysed_Dataset_May2026/2_GSE38959/DEGs_GSE38959_vs_Healthy_with_SV_significant.csv")
data3 <- read.csv("../../Analysed_Dataset_May2026/3_GSE65216/DEGs_GSE65216_vs_Healthy_with_SV_significant.csv")

head(data1)
common=read.csv("common_upregulated_genes.csv")




# ============================================================================
# FILTER SIGNIFICANT DEGs
# ============================================================================
# Function to extract significant DEGs
get_sig_degs <- function(data, dataset_name) {
  sig_degs <- data %>%
    filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
    mutate(
      Direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      Dataset = dataset_name
    ) %>%
    select(Gene, log2FoldChange, padj, Direction, Dataset)
  return(sig_degs)
}

# Apply to all datasets
sig_deg1 <- get_sig_degs(data1, "TCGA")
sig_deg2 <- get_sig_degs(data2, "GSE38959")
sig_deg3 <- get_sig_degs(data3, "GSE65216")

# Combine all significant DEGs
all_sig_degs <- bind_rows(sig_deg1, sig_deg2, sig_deg3)

# ============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# ============================================================================
# Define databases
databases = c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", 
             "GO_Cellular_Component_2023", "KEGG_2021_Human", 
             "BioCarta_2016", "WikiPathway_2023_Human", "Elsevier_Pathway_Collection")
abb = c("BP", "MF", "CC", "KEGG", "BioCarta", "WikiP", "Els")



########################
sig_deg1_up=dplyr::filter(sig_deg1,Direction=="Up")%>%pull(Gene)
sig_deg1_down=dplyr::filter(sig_deg1,Direction=="Down")%>%pull(Gene)

sig_deg2_up=dplyr::filter(sig_deg2,Direction=="Up")%>%pull(Gene)
sig_deg2_down=dplyr::filter(sig_deg2,Direction=="Down")%>%pull(Gene)


sig_deg3_up=dplyr::filter(sig_deg3,Direction=="Up")%>%pull(Gene)
sig_deg3_down=dplyr::filter(sig_deg3,Direction=="Down")%>%pull(Gene)


####################################
all_gene=list()
all_gene[[1]]=sig_deg1_up
all_gene[[2]]=sig_deg1_down
all_gene[[3]]=sig_deg2_up
all_gene[[4]]=sig_deg2_down
all_gene[[5]]=sig_deg3_up
all_gene[[6]]=sig_deg3_down

data_code=c("TCGA_Up","TCGA_Down","GSE38959_Up","GSE38959_Down","GSE65216_Up","GSE65216_Down")


i=1
j=1
# Store the results
all_pathways = data.frame()


for(i in 1:length(data_code))
{
  
  gene_list=all_gene[[i]]
  
  for (j in 1:length(databases)) {
   
    db=databases[j]
    
    enriched <- tryCatch({
      enrichr(gene_list, db)
    }, error = function(e) {
      message(paste("Error in", db, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(enriched) && length(enriched) > 0) {
      res_df <- enriched[[1]] %>%
        mutate(
          FDR = p.adjust(P.value, method = "fdr"),
          Dataset = data_code[i],
          #Direction = direction,
          Database = db,
          Gene_count = as.numeric(str_extract(Overlap, "^\\d+"))
        ) %>%
        filter(FDR < 0.05)
      
      if (nrow(res_df) > 0) {
        all_pathways <- rbind(all_pathways,res_df)
      }
    }
  }
  
}

head(all_pathways)

write.xlsx(all_pathways,"Table1.all_pathways.xlsx")

