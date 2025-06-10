rm(list=ls())

# Load required libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)

# Load PBMC multiomic dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Preprocess data for SCEGHiC (aggregating cells)
agg.data <- process_data(pbmc, max_overlap = 0.5)
rna <- agg.data[["rna"]]
atac <- agg.data[["atac"]]

# Remove peaks with zero counts across all cells
atac <- atac[rowSums(atac) != 0, ]

# Replace "-" with "_" in peak names to ensure consistent naming
rownames(atac) <- gsub("-", "_", rownames(atac))

# Load TSS annotation data for human hg38 genome
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Select highly variable genes present in TSS annotations
gene <- pbmc@assays[["SCT"]]@var.features
gene1 <- intersect(gene, tssdata$TargetGene)

# Subset RNA matrix for these genes
rna1 <- rna[which(rownames(rna) %in% gene1), ]

# Filter genes expressed in more than 2 cells
rna2 <- rna1[apply(rna1, 1, function(x) { sum(x != 0) > 2 }), ]

# Get peak names and replace "-" with "_"
peak <- pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
peak <- gsub("-", "_", peak)

# Initialize list to store results
results <- list()

# Loop through each gene
for (n in 1:dim(rna2)[1]) {
  G <- rownames(rna2)[n]
  y <- rna2[n, ]
  
  # Get chromosome and TSS position for the gene
  chr <- tssdata[which(G == tssdata$TargetGene), ]$chr
  start <- tssdata[which(G == tssdata$TargetGene), ]$TargetGeneTSS
  
  # Define promoter region ±1kb around TSS
  p1 <- paste(chr, ":", max(start - 1000, 0), "-", start + 1000, sep = "")
  
  # Define enhancer search region ±250kb around TSS
  p2 <- paste(chr, ":", max(start - 250000, 0), "-", start + 250000, sep = "")
  
  # Find overlapping peaks for promoters and enhancers
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  
  # Remove promoter peaks from enhancer peaks
  enhancers <- setdiff(enhancers, promoters)
  
  if (length(enhancers) > 0) {
    # Extract ATAC counts for enhancer peaks
    x <- atac[which(rownames(atac) %in% enhancers), ]
    
    # Combine gene expression and enhancer accessibility
    x2 <- rbind(y, x)
    rownames(x2) <- c("y", enhancers)
    
    # Calculate Spearman correlation between gene expression and enhancer accessibility
    dd <- cor(t(x2), method = "spearman")
    
    # Prepare results dataframe
    data <- data.frame(gene = G, peak2 = enhancers, score = dd[, 1][-1], method = "spearman")
    results[[n]] <- data
  } else {
    results[[n]] <- NULL
  }
  
  print(n)
}

# Combine all results into one dataframe
results <- do.call(rbind, results)

# Save to file
save(results, file = "Spearman_PBMC_peak_gene.rda")
