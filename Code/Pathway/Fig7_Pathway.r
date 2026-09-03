# ============================================================================
# LOAD LIBRARIES AND DATA
# ============================================================================
# Clear workspace
# rm(list = ls())
# Sys.setenv(http_proxy = "http://proxy.mh-hannover.de:8080")
# Sys.setenv(https_proxy = "http://proxy.mh-hannover.de:8080")
# options(download.file.method = "curl")
# options(download.file.extra = "-L --proxy http://proxy.mh-hannover.de:8080")
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
#library(enrichR)
library(ggsci)
library(stringr)  # Added for str_extract


common=read.csv("common_upregulated_genes.csv")
data=read.xlsx("Table1.all_pathways.xlsx")
head(data)

remove=c("BioCarta_2016", "Elsevier_Pathway_Collection",  "WikiPathway_2023_Human") 

data=dplyr::filter(data,!Database%in%remove)



# Prepare data for plotting
plot_data <- data%>%
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

plot_data_top10$Term=gsub("DNA replication","DNA Replication",plot_data_top10$Term)
plot_data_top10$Term=gsub("DNA Unwinding Involved In DNA Replication" ,"DNA Duplex Unwinding" ,plot_data_top10$Term)
plot_data_top10$Term=gsub("DNA-templated DNA Replication" ,"DNA Replication",plot_data_top10$Term)
plot_data_top10$Term=gsub("Regulation Of DNA Replication"   ,"DNA Replication",plot_data_top10$Term)

plot_data_top10$Term=gsub("Alcohol Dehydrogenase Activity, Zinc-Dependent"  ,"Alcohol Dehydrogenase Activity" ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Positive Regulation Of DNA-directed DNA Polymerase Activity"  ,  "Regulation Of DNA-directed DNA Polymerase Activity"  ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Water Transmembrane Transporter Activity"    ,  "Water Channel Activity"    ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Mitotic Spindle Assembly"     ,  "Mitotic Spindle Organization"     ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Mitotic Spindle Organization Checkpoint Signaling"     ,  "Mitotic Spindle Checkpoint Signaling"     ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Ethanol Oxidation"     ,   "Ethanol Metabolic Process"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Regulation Of Chromosome Separation"     ,   "Regulation Of Chromosome Segregation"     ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Glycerol Channel Activity"     ,   "Glycerol Transmembrane Transport"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Glycerol Transmembrane Transporter Activity"        ,   "Glycerol Transmembrane Transport"      ,plot_data_top10$Term)


plot_data_top10$Term=gsub("Regulation Of Chromosome Separation"     ,   "Regulation Of Chromosome Segregation"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Mitotic Sister Chromatid Segregation"     ,   "Regulation Of Chromosome Segregation"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Positive Regulation Of Chromosome Segregation"     ,   "Regulation Of Chromosome Segregation"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Sister Chromatid Segregation"     ,   "Regulation Of Chromosome Segregation"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Oxidoreductase Activity, Acting On The CH-OH Group Of Donors, NAD Or NADP As Acceptor"     ,   "Oxidoreductase Activity"      ,plot_data_top10$Term)


plot_data_top10$Term=gsub("Cell cycle"     ,   "Cell Cycle"      ,plot_data_top10$Term)
plot_data_top10$Term=gsub("Cell Cycle Overiew"     ,  "Cell Cycle"     ,plot_data_top10$Term)


pp=as.data.frame(unique(sort(plot_data_top10$Term)))

write.xlsx(plot_data_top10,"Table2.all_pathways.xlsx")
