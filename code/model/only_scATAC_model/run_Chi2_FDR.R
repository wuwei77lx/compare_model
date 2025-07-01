rm(list=ls())  # Clear workspace

# Load libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)

# Load PBMC multi-omics data
pbmc <- readRDS("PBMC_multiomic.rds")

# Subset for specific cell type, e.g., CD8+ T cells
celltype <- "CD8 T"
dataset <- subset(pbmc, idents = celltype)

# Extract ATAC counts matrix
atac <- dataset@assays[["peaks"]]@counts

# Binarize ATAC counts (presence/absence)
atac@x[atac@x > 0] <- 1

# Load transcription start site (TSS) annotations
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Select highly variable genes that overlap with TSS annotations
focus_markers <- pbmc@assays[["SCT"]]@var.features
focus_markers <- intersect(focus_markers, tssdata$TargetGene)

# Format peak names to match TSS annotation format
peak <- rownames(atac)
peak <- gsub("-", "_", peak)
rownames(atac) <- peak

results <- list()

for(n in seq_along(focus_markers)) {
  G <- focus_markers[n]
  chr <- tssdata[which(G == tssdata$TargetGene),]$chr
  start <- tssdata[which(G == tssdata$TargetGene),]$TargetGeneTSS
  
  # Define promoter region (±1kb around TSS)
  p1 <- paste(chr, ":", max(start-1000, 0), "-", start+1000, sep = "")
  
  # Define enhancer region (±250kb around TSS)
  p2 <- paste(chr, ":", max(start-250000, 0), "-", start+250000, sep = "")
  
  # Find peaks overlapping promoters and enhancers
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  enhancers <- setdiff(enhancers, promoters)
  
  x <- atac[which(rownames(atac) %in% enhancers), ]
  
  if(length(enhancers) > 1 & length(promoters) > 0) {
    
    # Aggregate promoter accessibility if multiple promoters
    if(length(promoters) > 1) {
      y <- colSums(atac[which(rownames(atac) %in% promoters), ])
      y[y > 0] <- 1
    } else {
      y <- atac[which(rownames(atac) %in% promoters), ]
    }
    
    # Perform chi-square test only if promoter accessible in some cells
    if(sum(y) > 0) {
      
      if(length(enhancers) == 1) {
        observed <- table(y, x)
        result <- chisq.test(observed, simulate.p.value = TRUE)
        results[[n]] <- data.frame(gene = G,
                                   peak1 = promoters[1],
                                   peak2 = enhancers,
                                   p_value = result$p.value,
                                   fdr = result$p.value)
      } else {
        res <- data.frame(gene = G,
                          peak1 = promoters[1],
                          peak2 = enhancers,
                          p_value = 0)
        for(m in seq_along(enhancers)) {
          z <- x[m, ]
          observed <- table(y, z)
          result <- chisq.test(observed, simulate.p.value = TRUE)
          res[m, "p_value"] <- result$p.value
        }
        # Adjust p-values for multiple testing
        res$fdr <- p.adjust(res$p_value, method = "BH")
        results[[n]] <- res
      }
      
    } else {
      results[[n]] <- NULL
    }
    
  } else {
    results[[n]] <- NULL
  }
  
  print(n)
}

# Combine results into a single data frame
results <- do.call(rbind, results)

# Save results
save(results, file = "chi2_FDR_PBMC_peak_gene_scATAC.rda")
