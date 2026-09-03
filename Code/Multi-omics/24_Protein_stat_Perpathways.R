## =============================================================================
## Multi-Omic Directionality Report: Proteomic, Blood, and Viral Associations
## For candidate genes: NMU, NEK2, HORMAD1, PRAME, KIF18B, MMP1
## =============================================================================

library(dplyr)
library(tidyr)
library(openxlsx)


data1=readRDS("Result_Combined/Combined_Proteomic_Results_clean.rds")

head(data1)
#############



library(openxlsx)
library(dplyr)
library(stringr)

## ============================================================
## 1. Read the RPPA antibody reference table
## ============================================================

prot_raw <- read.xlsx("Result_Combined/RPPA/41592_2013_BFnmeth2650_MOESM330_ESM.xlsx")

colnames(prot_raw) <- as.character(unlist(prot_raw[1, ]))
prot_clean <- prot_raw[-1, ] %>%
  select(
    Protein_Name = `Protein Name`,
    Gene_Name    = `Gene Name`,
    Validation   = `Antibody validation status`,
    Origin       = `Antibody Origin`,
    Pathway_Info = `Protein function and regulatory pathways`
  ) %>%
  filter(!is.na(Protein_Name))

## ============================================================
## 2. Clean, SHORT pathway name assignment based on keywords
##    (ordered by priority — first match wins)
## ============================================================



assign_pathway_short <- function(text) {
  t <- tolower(text)
  
  case_when(
    is.na(text) ~ "Unassigned",
    
    # --- mTOR complex components (very specific, check before generic PI3K/AKT) ---
    str_detect(t, "raptor|rictor|tsc1|tsc2|hamartin|tuberin|rheb|rapamycin-insensitive companion") ~ "mTOR Complex",
    
    # --- DNA Damage Repair ---
    str_detect(t, "dna repair|dna-dependent protein kinase|double-strand break|double strand break|dna damage checkpoint|dna damage sensor|dna-binding regulatory component|breast cancer susceptibility|brca") ~ "DNA Damage Repair",
    
    # --- p53 Pathway ---
    str_detect(t, "p53")                                                                            ~ "p53 Pathway",
    
    # --- Cell Cycle ---
    str_detect(t, "cell cycle|g1 phase|g1/s|late g1|early s phase|g2|mitosis|cyclin|retinoblastoma|e2f|cdks|restriction point|checkpoint signaling upon double strand breaks") ~ "Cell Cycle",
    
    # --- Apoptosis (strict: named apoptosis regulators) ---
    str_detect(t, "pro-apoptotic|anti-apoptotic|apoptosis suppressor|execution of apoptosis|bcl-2 interacting mediator|bim|bad|bax|bcl-xl")  ~ "Apoptosis",
    
    # --- Autophagy ---
    str_detect(t, "autophagy|beclin")                                                               ~ "Autophagy",
    
    # --- PI3K/AKT/mTOR (generic) ---
    str_detect(t, "pi3k|akt|mtor|p70 s6 kinase|s6 kinase|pten|pras40|pdk1|phosphoinositide-dependent|irs molecules|insulin signaling and play a central role") ~ "PI3K/AKT/mTOR",
    
    # --- Protein Translation ---
    str_detect(t, "elongation factor|initiation factor|eif|eef|peptidyl-trna|protein translation|helicase activity|ribosomal|translation initiation|eukaryotic initiation") ~ "Protein Translation",
    
    # --- MAPK/RAS ---
    str_detect(t, "mapk|erk|mek|raf|p38|ap-1|myc/max/mad|stress-activated protein kinase|mkk|jnk|extracellular-signal-related kinase|h-ras|k-ras|n-ras|ras oncogene|ras to activate the mek") ~ "RAS/MAPK",
    
    # --- RTK / growth factor receptors ---
    str_detect(t, "receptor tyrosine kinase|erbb|egfr|her2|vegf|tyrosine kinase receptor|growth factor recept|insulin receptor|igfbp|insulin-like growth factor|hepatocyte growth factor|transferrin receptor|met proto-oncogene|src family|src, akt and mapk|scaffold protein in signaling for a variety of receptor|adaptor protein.*receptor tyrosine kinase|recruited by.*receptor tyrosine") ~ "RTK Signaling",
    
    # --- Metabolism/Insulin ---
    str_detect(t, "insulin signaling|insulin-sensitive|glucose transport|energy homeostasis|fatty acid|glycolysis|glyceraldehyde|pentose phosphate|stearoyl-coa|asparagine synthesis|oleate|malonyl-coa|monounsaturated fatty") ~ "Metabolism",
    
    # --- Wnt signaling ---
    str_detect(t, "wnt|downstream of.*wnt")                                                        ~ "Wnt Signaling",
    
    # --- TGF-beta / TNF ---
    str_detect(t, "tgf beta|tgf-beta|tgfb|tnf beta superfamily|mediates tgf")                       ~ "TGF-beta/TNF Signaling",
    
    # --- Hormone/nuclear receptors ---
    str_detect(t, "androgen receptor|estrogen receptor|progesterone receptor|hormone regulation")   ~ "Hormone Signaling",
    
    # --- JAK/STAT ---
    str_detect(t, "signal transducer and activator of transcription|stat3|stat5|stat")              ~ "JAK/STAT Signaling",
    
    # --- NF-kB ---
    str_detect(t, "nf-κb|nf-kb|nuclear factor κ b|rel family")                                     ~ "NF-kB Signaling",
    
    # --- GTPase / small G proteins ---
    str_detect(t, "gtpase|ras superfamily|rab11|rab25|small ras-like|heterotrimric g protein|guanine nucleotide exchange") ~ "GTPase Signaling",
    
    # --- PKC signaling ---
    str_detect(t, "protein kinase c|pkc|diacylglycerol|phorbol ester")                             ~ "PKC Signaling",
    
    # --- Cell adhesion / cytoskeleton / ECM ---
    str_detect(t, "cell-cell adhesion|cadherin|claudin|tight junction|collagen|extracellular matrix|glycoprotein|adhesion|actin|myosin|integrin|cytoskeleton|cell motility|caveolae|tubulin|transglutaminase|focal adhesion|morphology") ~ "Adhesion/Cytoskeleton",
    
    # --- Development / cell fate receptors ---
    str_detect(t, "transmembrane receptor family|cell fate|osteoblast|development")                 ~ "Development/Cell Fate",
    
    # --- Chaperone/stress ---
    str_detect(t, "molecular chaperone|heat shock|protein homeostasis|environmental stress")       ~ "Chaperone/Stress Response",
    
    # --- Plasminogen Activation ---
    str_detect(t, "plasminogen|pai-1|tpa|upa")                                                     ~ "Plasminogen Activation",
    
    # --- Angiogenesis/endothelial ---
    str_detect(t, "endothelia|angiogenesis|vascular")                                              ~ "Angiogenesis",
    
    # --- Immune/Inflammation ---
    str_detect(t, "lymphocyte|immune|b cell|inflammat|tnf")                                        ~ "Immune/Inflammation",
    
    # --- 14-3-3 scaffolding ---
    str_detect(t, "14-3-3")                                                                        ~ "Scaffolding/Multi-functional",
    
    # --- Transcription/Chromatin (generic TFs, checked after specific families above) ---
    str_detect(t, "transcription factor|transfection factor|y-box|transcriptional regulator|transcriptional activity|transcriptional co-activator|x-box binding|repressor function in.*signaling|splicing factor|ets family|ets-1") ~ "Transcription/Chromatin",
    
    # --- Tumor suppressor (generic catch) ---
    str_detect(t, "tumor suppressor")                                                             ~ "Tumor Suppressor",
    
    # --- Growth/differentiation generic ---
    str_detect(t, "growth, differentiation, and cell survival|upregulated with cell growth|ubiquitous cytosolic phosphoprotein|poor prognosis") ~ "Growth/Differentiation",
    
    TRUE ~ "Other"
  )
}




# Apply to data
data2 <- data2 %>%
  mutate(PathwayShort = assign_pathway_short(Pathway))

# View distribution
table(data2$PathwayShort)

# Inspect any stragglers in "Other"
data2 %>% 
  filter(PathwayShort == "Other") %>% 
  select(Protein, Pathway) %>%
  distinct() %>%
  head(20)

data2 <- data2 %>%
  mutate(PathwayShort = assign_pathway_short(Pathway))

table(data2$PathwayShort)
prot_clean <- prot_clean %>%
  mutate(Pathway = assign_pathway_short(Pathway_Info))

## ============================================================
## 3. Check the distribution of assigned pathways
## ============================================================

print(table(prot_clean$Pathway))

## ============================================================
## 4. Build Protein name matching your data1$Protein convention
##    (Protein_Name-Origin_code-Validation_code)
## ============================================================

suffix_map <- c("Validated" = "V", "Use with Caution" = "C", "Under Evaluation" = "U")
origin_map <- c("Rabbit" = "R", "Mouse" = "M", "Goat" = "G")

prot_clean <- prot_clean %>%
  mutate(
    val_code    = suffix_map[Validation],
    origin_code = origin_map[Origin],
    Protein     = paste0(Protein_Name, "-", origin_code, "-", val_code)
  )

## ============================================================
## 5. Final clean mapping table
## ============================================================

pathway_map <- prot_clean %>%
  select(Protein, Protein_Name, Gene_Name, Validation, Pathway) %>%
  distinct()

head(pathway_map)

## ============================================================
## 6. Merge into your data1
## ============================================================

data1 <- data1 %>%
  left_join(pathway_map %>% select(Protein, Pathway), by = "Protein")

# Check unmatched
sum(is.na(data1$Pathway))
## ============================================================
## 7. Save
## ============================================================

write.csv(pathway_map, "Result_Combined/RPPA/RPPA_Protein_Pathway_Map_Clean.csv", row.names = FALSE)
write.csv(data1, "Result_Combined/RPPA/data1_with_Pathway.csv", row.names = FALSE)

