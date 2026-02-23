
rm(list = ls())
gc()


# ==============================
# Load libraries
# ==============================
library(GEOquery)
library(DESeq2)
library(sva)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(WGCNA)

# ==============================
# Load DESeq2 object
# ==============================
dds <- readRDS("GSE38959_dds_final.rds")

# ==============================
# Variance stabilizing transformation
# ==============================
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Transpose: rows = samples, columns = genes
datExpr <- t(vsd_mat)

# ==============================
# Check data quality
# ==============================
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ==============================
# Prepare traits
# ==============================
traits <- data.frame(condition = as.factor(dds$condition))
rownames(traits) <- rownames(datExpr)

# ==============================
# Soft-threshold selection
# ==============================
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
saveRDS(sft,"sft.rds")


png(filename = "Pic_1_WGCNA_soft_threshold.png",
    width = 6, height = 5, units = "in", res = 300)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="tomato")
abline(h=0.9, col="tomato")
dev.off()
#### see the image and take the value that is very close to the line in this take 5 as it is lowest and have above 90% confidence

# ==============================
# Step 1: Build adjacency and TOM
# ==============================
softPower <- 5  # from your soft-threshold selection  from Pic1 anove
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# ==============================
# Step 2: Hierarchical clustering
# ==============================
geneTree <- hclust(as.dist(dissTOM), method = "average")

# ==============================
# Step 3: Module detection
# ==============================
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

# ==============================
# Step 4: Save dendrogram + module colors
# ==============================
png("Pic_2_WGCNA_gene_dendrogram.png", width = 10, height = 6, units = "in", res = 300)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# ==============================
# Step 5: Module eigengenes & trait correlation
# ==============================
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

saveRDS(MEs,"MEs.rds")


#https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
###Module Merging
ME.dissimilarity = 1-cor(MEList$eigengenes, use="complete") #Calculate eigengene dissimilarity
# ==============================
# Step: Cluster module eigengenes and save as PNG
# ==============================
# Hierarchical clustering of module eigengenes
METree <- hclust(as.dist(ME.dissimilarity), method = "average")

# Save dendrogram as 300 DPI PNG
png("Pic_3_WGCNA_module_eigengene_dendrogram.png",
    width = 10, height = 6, units = "in", res = 300)

par(mar = c(0, 4, 2, 0))   # set margins
par(cex = 0.6)             # scale text
plot(METree, main = "Clustering of Module Eigengenes")
abline(h = 0.25, col = "red")  # merge threshold (correlation 0.75)

dev.off()


#####

# ==============================
# Step: Dendrogram with original and merged module colors
# ==============================

# Merge similar modules (already done)
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors       # new module assignments
mergedMEs <- merge$newMEs          # eigengenes of merged modules

# Save dendrogram as PNG
png("Pic_4_WGCNA_merged_modules_dendrogram.png",
    width = 10, height = 6, units = "in", res = 300)

# Plot gene dendrogram with original and merged module colors
plotDendroAndColors(geneTree, 
                    cbind(dynamicColors, mergedColors),
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors: original vs merged")

dev.off()



#####

# ==============================
# 1. Dynamic Module-Trait and Hub Gene Analysis
# ==============================


# ==============================
# Module-Trait Correlation Heatmap
# ==============================
# Convert all traits to numeric for correlation
trait_numeric <- data.frame(lapply(traits, function(x) {
  if (is.factor(x)) as.numeric(x) else x
}))
rownames(trait_numeric) <- rownames(datExpr)

# Correlation and p-values
moduleTraitCor <- cor(MEs, trait_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Text matrix for heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Save heatmap as PNG
png("WGCNA_module_trait_heatmap_dynamic.png",
    width = 12, height = 8, units = "in", res = 300)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(trait_numeric),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

# ==============================
# Dynamic Hub Gene Calculation
# ==============================
hubGenesDF <- data.frame(Gene = colnames(datExpr), Module = dynamicColors)

# Get unique conditions/factors
conditions <- colnames(trait_numeric)

# Loop over each trait/condition and calculate GS and p-values
for (cond in conditions) {
  cat("Processing trait:", cond, "\n")
  
  trait_vector <- trait_numeric[[cond]]
  
  # Gene significance and p-values
  GS <- as.data.frame(cor(datExpr, trait_vector, use = "p"))
  GS_p <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(datExpr)))
  
  # Add columns dynamically
  GS_colname <- paste0("GS_", cond)
  GS_p_colname <- paste0("GS_pvalue_", cond)
  
  hubGenesDF[[GS_colname]] <- GS[,1]
  hubGenesDF[[GS_p_colname]] <- GS_p[,1]
}

# ==============================
# Module Membership (kME) Calculation
# ==============================
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# Extract module names without "ME" prefix
modNames <- substring(names(MEs), 3, nchar(names(MEs)))

# kME and p-value for each gene's assigned module
kME <- numeric(ncol(datExpr))
kME_p <- numeric(ncol(datExpr))
for (i in 1:ncol(datExpr)) {
  mod_idx <- match(dynamicColors[i], modNames)
  kME[i] <- geneModuleMembership[i, mod_idx]
  kME_p[i] <- MMPvalue[i, mod_idx]
}

hubGenesDF$kME <- kME
hubGenesDF$kME_pvalue <- kME_p

# Sort by module and kME
hubGenesDF <- hubGenesDF %>% arrange(Module, -kME)

# Save hub genes CSV
write.csv(hubGenesDF, "WGCNA_hub_genes_dynamic.csv", row.names = FALSE)
cat("Hub gene table saved: WGCNA_hub_genes_dynamic.csv\n")



dir.create("Images")

# ==============================
# Prepare geneTraitSignificance
# ==============================
# Assuming your trait is the first column in trait_numeric
geneTraitSignificance <- as.matrix(cor(datExpr, trait_numeric[,1], use = "p"))
rownames(geneTraitSignificance) <- colnames(datExpr)

# ==============================
# Create folder for module plots
# ==============================
if (!dir.exists("Images")) dir.create("Images")

# ==============================
# Loop through all modules
# ==============================
moduleResults <- data.frame(Module = character(),
                            Correlation = numeric(),
                            PValue = numeric(),
                            stringsAsFactors = FALSE)

for (module in modNames) {
  
  column <- match(module, modNames)
  moduleGenes <- mergedColors == module
  
  if (sum(moduleGenes) < 2) {
    warning(paste("Module", module, "has less than 2 genes. Skipping correlation."))
    next
  }
  
  MM <- abs(geneModuleMembership[moduleGenes, column])
  GS <- abs(geneTraitSignificance[moduleGenes, 1])
  
  # Remove NAs
  valid <- complete.cases(MM, GS)
  MM <- MM[valid]
  GS <- GS[valid]
  
  if (length(MM) < 2) {
    warning(paste("Module", module, "does not have enough valid genes. Skipping."))
    next
  }
  
  # Compute correlation and p-value
  corTest <- cor.test(MM, GS)
  corValue <- corTest$estimate
  pValue <- corTest$p.value
  
  # Save results
  moduleResults <- rbind(moduleResults, data.frame(Module = module,
                                                   Correlation = corValue,
                                                   PValue = pValue))
  
  # Plot
  png(filename = paste0("Images/ModuleMembership_vs_GeneSignificance_", module, ".png"),
      width = 6, height = 5, units = "in", res = 300)
  
  par(mar = c(1,1,1,1))
  verboseScatterplot(MM, GS,
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for trait",
                     main = paste("Module membership vs. gene significance\n", module),
                     cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5,
                     col = module)
  
  #legend("topleft", legend = paste0("cor = ", round(corValue, 2), ", p = ", signif(pValue, 3)), bty="n")
  dev.off()
}

# Save results table
write.csv(moduleResults, "WGCNA_ModuleMembership_vs_GeneSignificance.csv", row.names = FALSE)
print(moduleResults)