## =============================================================================
## Multi-Omic Directionality Report: Proteomic, Blood, and Viral Associations
## For candidate genes: NMU, NEK2, HORMAD1, PRAME, KIF18B, MMP1
## =============================================================================

library(dplyr)
library(tidyr)
library(openxlsx)
data1=readRDS("Result_Combined/Combined_Blood_Results_clean.rds")
write.xlsx(data1,"Result_Combined/Combined_Blood_Results_clean.xlsx")


data1=readRDS("Result_Combined/Combined_eQTM_Results_ALL.rds")
data2=readRDS("Result_Combined/Combined_eQTM_Results_ALL_annotated.rds")
write.xlsx(data2,"Result_Combined/Combined_eQTM_Results_ALL.xlsx")


data1=readRDS("Result_Combined/Combined_Proteomic_Results_clean.rds")
write.xlsx(data1,"Result_Combined/Combined_Proteomic_Results_clean.xlsx")


data1=readRDS("Result_Combined/Combined_Viral_Results_clean.rds")
write.xlsx(data1,"Result_Combined/Combined_Viral_Results_clean.xlsx")



rm(list=ls())

data1=readRDS("RNAseq_DNAMeth//rnaseq_aligned.rds")
