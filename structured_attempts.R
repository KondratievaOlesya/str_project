# Load necessary libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)

# -----------------------------------------------
# Try 1: Initial Data Loading and Matching
# -----------------------------------------------

# Load clinical data
clinical_data <- read.csv("coad_clinical.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", na.strings = c("NA", ""))

# Assign row names from Tumor_Sample_Barcode
rownames(clinical_data) <- clinical_data$Tumor_Sample_Barcode
clinical_data$Tumor_Sample_Barcode <- NULL

# Load count data
count_data <- read.csv("gene_expression/coad/count.csv", row.names = 1)

# Match clinical data to count data
clinical_data <- clinical_data[match(colnames(count_data), rownames(clinical_data)), ]

# Check alignment
stopifnot(identical(rownames(clinical_data), colnames(count_data)))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = clinical_data, design = ~ MSI_TMB)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results for MSS/TMB-H vs MSS/TMB-L
results_MSS_TMBH_vs_L <- results(dds, contrast = c("MSI_TMB", "MSS/TMB-H", "MSS/TMB-L"))

# Summary of results
summary(results_MSS_TMBH_vs_L)

# MA plot
plotMA(results_MSS_TMBH_vs_L, ylim = c(-5, 5))

# Save results
write.csv(as.data.frame(results_MSS_TMBH_vs_L), file = "DESeq2_results_MSS_TMBH_vs_L.csv")

# View top differentially expressed genes
head(results_MSS_TMBH_vs_L[order(results_MSS_TMBH_vs_L$padj), ])

# -----------------------------------------------
# Try 2: Checking Mismatched Samples and Data Issues
# -----------------------------------------------

# Check mismatches between count_data and clinical_data
setdiff(colnames(count_data), rownames(clinical_data))  # Samples in count_data but not in clinical_data
setdiff(rownames(clinical_data), colnames(count_data))  # Samples in clinical_data but not in count_data

# Reload clinical data
clinical_data <- read.csv("coad_clinical.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", na.strings = c("NA", ""))

# Set row names from Tumor_Sample_Barcode, handle missing
if ("Tumor_Sample_Barcode" %in% colnames(clinical_data)) {
  rownames(clinical_data) <- clinical_data$Tumor_Sample_Barcode
  clinical_data$Tumor_Sample_Barcode <- NULL
} else {
  stop("Tumor_Sample_Barcode column missing!")
}

# Remove rows with NA or empty row names
clinical_data <- clinical_data[!is.na(rownames(clinical_data)) & rownames(clinical_data) != "", ]

# Reload count data
count_data <- read.csv("gene_expression/coad/count.csv", row.names = 1)

# Check mismatches again
setdiff(colnames(count_data), rownames(clinical_data))
setdiff(rownames(clinical_data), colnames(count_data))

# Subset both datasets to matching samples
matching_samples <- intersect(colnames(count_data), rownames(clinical_data))
count_data <- count_data[, matching_samples]
clinical_data <- clinical_data[matching_samples, ]

# Re-check alignment
stopifnot(identical(rownames(clinical_data), colnames(count_data)))

# Create DESeqDataSet object again
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = clinical_data, design = ~ MSI_TMB)

# Run DESeq2 analysis
dds <- DESeq(dds)

# -----------------------------------------------
# Try 3: Ensure Data is Numeric and Handle NA Values in Count Data
# -----------------------------------------------

# Reload count data again
count_data <- read.csv("gene_expression/coad/count.csv", row.names = 1, stringsAsFactors = FALSE)

# Ensure all data in count_data is numeric
count_data <- as.data.frame(lapply(count_data, function(x) as.numeric(as.character(x))))

# Check for NA values and handle them
if (any(is.na(count_data))) {
  count_data[is.na(count_data)] <- 0  # Replace NA with 0
}

# Check mismatches between count_data and clinical_data again
setdiff(colnames(count_data), rownames(clinical_data))
setdiff(rownames(clinical_data), colnames(count_data))

# Subset datasets to matching samples
matching_samples <- intersect(colnames(count_data), rownames(clinical_data))
count_data <- count_data[, matching_samples]
clinical_data <- clinical_data[matching_samples, ]

# Re-check alignment
stopifnot(identical(rownames(clinical_data), colnames(count_data)))

# Create DESeqDataSet object again
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = clinical_data, design = ~ MSI_TMB)

# Run DESeq2 analysis
dds <- DESeq(dds)

# -----------------------------------------------
# Try 4: Plotting and More Advanced Steps (PCA, Heatmaps)
# -----------------------------------------------

# Extract results and apply independent filtering
res <- results(dds)
res <- lfcShrink(dds, coef = "MSI_TMB", type = "apeglm")

# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Save results to CSV
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

# MA plot
plotMA(res, main = "DESeq2 MA Plot", ylim = c(-2, 2))

# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "MSI_TMB")

# Heatmap of top 20 differentially expressed genes
top_genes <- head(order(res$padj), 20)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)[, c("MSI_TMB")]))

# -----------------------------------------------
# Try 5: Further Data Checks and Cleaning
# -----------------------------------------------

# Investigate mismatched samples
mismatched_in_count <- setdiff(colnames(count_data), rownames(clinical_data))
mismatched_in_clinical <- setdiff(rownames(clinical_data), colnames(count_data))

cat("Samples in count_data but not in clinical_data:\n")
print(mismatched_in_count)

cat("\nSamples in clinical_data but not in count_data:\n")
print(mismatched_in_clinical)

# Inspect sample identifiers
head(colnames(count_data))  # Samples in count_data
head(rownames(clinical_data))  # Samples in clinical_data

# Reload clinical and count data again and recheck alignment
clinical_data <- read.csv("coad_clinical.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8", na.strings = c("NA", ""))

if ("Tumor_Sample_Barcode" %in% colnames(clinical_data)) {
  rownames(clinical_data) <- clinical_data$Tumor_Sample_Barcode
  clinical_data$Tumor_Sample_Barcode <- NULL
}

count_data <- read.csv("gene_expression/coad/count.csv", row.names = 1)

matching_samples <- intersect(colnames(count_data), rownames(clinical_data))
count_data <- count_data[, matching_samples]
clinical_data <- clinical_data[matching_samples, ]

stopifnot(identical(rownames(clinical_data), colnames(count_data)))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = clinical_data, design = ~ MSI_TMB)

# Continue with DESeq2 analysis...

# -----------------------------------------------
# Try 6: Investigating Top Genes and Heatmap Plotting
# -----------------------------------------------

# Plot heatmap of top 20 differentially expressed genes
top_genes <- head(order(res$padj), 20)
heatmap_data <- assay(vsd)[top_genes, ]
heatmap_data <- heatmap_data - rowMeans(heatmap_data)

# Plot heatmap
pheatmap(heatmap_data, annotation_col = as.data.frame(colData(dds)[, "MSI_TMB"]))


