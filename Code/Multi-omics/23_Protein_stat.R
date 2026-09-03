## =============================================================================
## Multi-Omic Directionality Report: Proteomic, Blood, and Viral Associations
## For candidate genes: NMU, NEK2, HORMAD1, PRAME, KIF18B, MMP1
## =============================================================================

library(dplyr)
library(tidyr)
library(openxlsx)


data1=readRDS("Result_Combined/Combined_Proteomic_Results_clean.rds")

head(data1)

library(dplyr)
library(tidyr)
library(stringr)

## -------------------------------------------------------------------
## 1. Extract antibody suffix (validation/species info) using regex
## -------------------------------------------------------------------
data1 <- data1 %>%
  mutate(
    AntibodySuffix = str_extract(Protein, "-[RM]-[CEV]$"),
    Species        = case_when(
      str_detect(AntibodySuffix, "^-R-") ~ "Rabbit",
      str_detect(AntibodySuffix, "^-M-") ~ "Mouse",
      TRUE ~ NA_character_
    ),
    ValidationType = case_when(
      str_detect(AntibodySuffix, "-C$") ~ "Validated-C",
      str_detect(AntibodySuffix, "-E$") ~ "Validated-E",
      str_detect(AntibodySuffix, "-V$") ~ "Validated-V",
      TRUE ~ NA_character_
    ),
    ProteinClean = str_remove(Protein, "-[RM]-[CEV]$")  # strip suffix for clean protein name
  )

## -------------------------------------------------------------------
## 2. Detect phosphoproteins vs total proteins using substr/grepl
## -------------------------------------------------------------------
data1 <- data1 %>%
  mutate(
    IsPhospho = grepl("_p[STY][0-9]", Protein),   # e.g., _pS79, _pT183, _pY204
    Phosphosite = str_extract(Protein, "_p[STY][0-9]+(_[STY][0-9]+)*")
  )

## -------------------------------------------------------------------
## 3. Pathway keyword tagging via grepl (customize keyword lists)
## -------------------------------------------------------------------
data1 <- data1 %>%
  mutate(
    PathwayGroup = case_when(
      grepl("MEK|MAPK|JNK|ERK|p38", ProteinClean, ignore.case = TRUE) ~ "MAPK/JNK signaling",
      grepl("mTOR|S6|4E-BP1|AKT|PI3K|ACC", ProteinClean, ignore.case = TRUE) ~ "PI3K/AKT/mTOR",
      grepl("Ku80|Mre11|MSH2|MSH6|BRCA|ATM|ATR|H2AX|Rad51", ProteinClean, ignore.case = TRUE) ~ "DNA damage/repair",
      grepl("N-Ras|K-Ras|H-Ras|Raf", ProteinClean, ignore.case = TRUE) ~ "RAS/RAF",
      grepl("Myosin|MYH11|N-Cadherin|Vimentin|E-Cadherin", ProteinClean, ignore.case = TRUE) ~ "Cytoskeleton/EMT",
      grepl("NF-kB|IkB", ProteinClean, ignore.case = TRUE) ~ "NF-kB signaling",
      grepl("Lck|LKB1", ProteinClean, ignore.case = TRUE) ~ "Kinase signaling (other)",
      grepl("14-3-3", ProteinClean, ignore.case = TRUE) ~ "14-3-3 scaffold",
      grepl("NDRG1|MIG-6", ProteinClean, ignore.case = TRUE) ~ "Stress response",
      TRUE ~ "Other/Unclassified"
    )
  )

## -------------------------------------------------------------------
## 4. Reshape to wide (Protein x Gene) to assign Panel I/II/III
## -------------------------------------------------------------------
wide_df <- data1 %>%
  select(Gene, Protein, ProteinClean, PathwayGroup, Species, ValidationType, IsPhospho, Estimate) %>%
  pivot_wider(names_from = Gene, values_from = Estimate)

gene_cols <- c("NMU", "NEK2", "HORMAD1", "PRAME", "KIF18B", "MMP1")

wide_df <- wide_df %>%
  rowwise() %>%
  mutate(
    n_pos = sum(c_across(all_of(gene_cols)) > 0, na.rm = TRUE),
    n_neg = sum(c_across(all_of(gene_cols)) < 0, na.rm = TRUE),
    Panel = case_when(
      n_neg == 6 ~ "Panel I (Negative)",
      n_pos == 6 ~ "Panel II (Positive)",
      TRUE ~ "Panel III (Mixed)"
    )
  ) %>%
  ungroup()

## -------------------------------------------------------------------
## 5. Summaries
## -------------------------------------------------------------------
table(wide_df$Panel)                          # should approximate 83/89/109
table(wide_df$Panel, wide_df$PathwayGroup)
table(wide_df$Panel, wide_df$Species)
table(wide_df$Panel, wide_df$ValidationType)
table(wide_df$Panel, wide_df$IsPhospho)

## Fisher's test for enrichment
fisher.test(table(wide_df$Panel, wide_df$PathwayGroup), simulate.p.value = TRUE)
fisher.test(table(wide_df$Panel, wide_df$IsPhospho))

fisher.test(table(wide_df$Panel, wide_df$Species))


# Test species enrichment
fisher.test(table(wide_df$Panel, wide_df$Species))

# Test validation type enrichment  
fisher.test(table(wide_df$Panel, wide_df$ValidationType), simulate.p.value = TRUE)

# Sanity check: where did N-Ras land?
wide_df %>% filter(grepl("N-Ras", Protein)) %>% select(Protein, all_of(gene_cols), Panel)

# Sanity check: DNA repair proteins
wide_df %>% filter(PathwayGroup == "DNA damage/repair") %>% 
  select(Protein, all_of(gene_cols), Panel)

