library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # for 450K array
# OR if you used EPIC array, use instead:
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## ------------------------------------------------------------
## Load annotation
## ------------------------------------------------------------
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # if EPIC

anno_df <- as.data.frame(anno) %>%
  tibble::rownames_to_column("CpG") %>%
  select(CpG, chr, pos, strand, Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Group)

cat("Annotation loaded. Dim:", dim(anno_df), "\n")
head(anno_df)

## ------------------------------------------------------------
## Load your eQTM results
## ------------------------------------------------------------
data <- readRDS("Result_Combined/Combined_eQTM_Results_ALL.rds")

cat("eQTM data loaded. Dim:", dim(data), "\n")

## ------------------------------------------------------------
## Merge CpG annotation info into eQTM results
## ------------------------------------------------------------
data_annotated <- data %>%
  left_join(anno_df, by = "CpG")

cat("Merged data dim:", dim(data_annotated), "\n")
cat("Number of CpGs matched to annotation:", sum(!is.na(data_annotated$chr)), "out of", nrow(data_annotated), "\n")

## ------------------------------------------------------------
## Check for any unmatched CpGs
## ------------------------------------------------------------
unmatched <- data_annotated %>% filter(is.na(chr)) %>% distinct(CpG)
if (nrow(unmatched) > 0) {
  cat("Warning: The following CpGs were NOT found in annotation:\n")
  print(unmatched)
}

## ------------------------------------------------------------
## Save annotated results
## ------------------------------------------------------------
saveRDS(data_annotated, "Result_Combined/Combined_eQTM_Results_ALL_annotated.rds")
write.csv(data_annotated, "Result_Combined/Combined_eQTM_Results_ALL_annotated.csv", row.names = FALSE)

cat("Done. Annotated file saved with columns:\n")
print(colnames(data_annotated))


table(data_annotated$chr)
table(data_annotated$Relation_to_Island)
