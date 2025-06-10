rm(list=ls())

# Load required libraries
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load multiomic PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Set the default assay to ATAC peaks
DefaultAssay(pbmc) <- "peaks"

# Compute GC content for each peak region (needed by LinkPeaks)
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# Link peaks to genes using Signac's LinkPeaks function
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 250000,      # max distance between peak and gene TSS
  pvalue_cutoff = 1.01,  # use all p-values (cutoff slightly above 1 to keep all)
  score_cutoff = 0       # include all links with score >= 0
)

# Extract the peak-gene links dataframe
signac <- pbmc@assays[["peaks"]]@links

# Save results to file
save(signac, file = "Signac_PBMC_peak_gene.rda")
