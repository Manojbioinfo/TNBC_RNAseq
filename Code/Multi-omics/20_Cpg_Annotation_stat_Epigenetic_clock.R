## =============================================================================
## Direction of CpG-Gene (eQTM) Associations by CpG Island Context
## Clean analysis + publication-quality plot + summary report
## =============================================================================

library(dplyr)
library(ggplot2)
library(scales)
library(openxlsx)

data=readRDS("Result_Combined/Combined_eQTM_Results_ALL_annotated.rds")
epi=read.xlsx("COVID_Manoj/10522_2025_10360_MOESM2_ESM.xlsx", sheet = 4)
intersect(epi$X2,data$CpG)
