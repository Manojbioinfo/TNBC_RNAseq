############################################################
## 1. Load library
############################################################
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("sva", ask = FALSE)

library(sva)

############################################################
## 2. Input data
############################################################
# countm  = gene x sample matrix (counts or expression)
# coldata = sample metadata

countm  <- as.matrix(countm)
coldata <- as.data.frame(coldata)

# Ensure sample order matches
all(colnames(countm) == rownames(coldata))  # should be TRUE

############################################################
## 3. Format variables
############################################################
coldata$Disease   <- factor(coldata$Disease)
coldata$Ethnicity <- factor(coldata$Ethnicity)
coldata$Age       <- as.numeric(coldata$Age)

############################################################
## 4. (Optional but common) log transform counts
############################################################
# Recommended if using raw RNA-seq counts
expr <- log2(countm + 1)

# If already normalized/logged, use instead:
# expr <- countm

############################################################
## 5. Build model matrices
############################################################
mod  <- model.matrix(~ Disease + Age + Ethnicity, data = coldata)
mod0 <- model.matrix(~ Age + Ethnicity, data = coldata)

############################################################
## 6. Run SVA
############################################################
svobj <- sva(expr, mod, mod0)

# Number of surrogate variables
n_sv <- svobj$n.sv
cat("Number of SVs detected:", n_sv, "\n")

############################################################
## 7. Extract SVs
############################################################
SVs <- svobj$sv
colnames(SVs) <- paste0("SV", seq_len(n_sv))

head(SVs)

############################################################
## 8. Add SVs to metadata (for downstream models)
############################################################
coldata_sv <- cbind(coldata, SVs)

############################################################
## 9. Save outputs
############################################################
# Save SV table
write.csv(SVs, "SurrogateVariables.csv")

# Save full sva object as RDS
saveRDS(svobj, file = "svobj.rds")