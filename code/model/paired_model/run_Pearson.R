rm(list=ls())

# Load required libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)

# Load multi-omics PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Data preprocessing - aggregation for SCEGHiC reference
agg.data <- process_data(pbmc, max_overlap = 0.5)
rna <- agg.data[["rna"]]
atac <- agg.data[["atac"]]
atac <- atac[rowSums(atac) != 0, ]  # remove peaks with zero counts
rownames(atac) <- gsub("-", "_", rownames(atac))  # standardize peak names

# Load transcription start sites (TSS) annotation for human hg38
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Select highly variable genes present in TSS annotation
gene <- pbmc@assays[["SCT"]]@var.features
gene1 <- intersect(gene, tssdata$TargetGene)
rna1 <- rna[which(rownames(rna) %in% gene1), ]
rna2 <- rna1[apply(rna1, 1, function(x) { sum(x != 0) > 2 }), ]  # keep genes expressed in >2 cells

peak <- rownames(atac)

# Calculate Pearson correlations between gene expression and peaks in ± 250 kb window
results <- list()
for (n in 1:nrow(rna2)) {
  G <- rownames(rna2)[n]
  y <- rna2[n, ]
  
  # Get chromosome and TSS position for gene G
  chr <- tssdata[which(G == tssdata$TargetGene), ]$chr
  start <- tssdata[which(G == tssdata$TargetGene), ]$TargetGeneTSS
  
  # Define promoter region ± 1 kb around TSS
  p1 <- paste0(chr, ":", max(start - 1000, 0), "-", start + 1000)
  # Define enhancer search window ±250 kb around TSS
  p2 <- paste0(chr, ":", max(start - 250000, 0), "-", start + 250000)
  
  # Find peaks overlapping promoter and enhancer windows
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  
  # Remove promoter peaks from enhancers
  enhancers <- setdiff(enhancers, promoters)
  
  if (length(enhancers) > 0) {
    x <- atac[which(rownames(atac) %in% enhancers), ]
    x2 <- rbind(y, x)
    rownames(x2) <- c("y", enhancers)
    dd <- cor(t(x2))  # compute Pearson correlation
    data <- data.frame(
      gene = G,
      peak2 = enhancers,
      score = dd[, 1][-1],  # correlations between gene and each enhancer
      method = "pearson"
    )
    results[[n]] <- data
  } else {
    results[[n]] <- NULL
  }
  
  print(n)
}

# Combine results for all genes
results <- do.call(rbind, results)

# Save results
save(results, file = "Pearson_PBMC_peak_gene.rda")

