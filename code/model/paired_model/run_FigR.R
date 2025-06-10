rm(list=ls())

# Load required libraries
library(Seurat)
library(Signac)
library(FigR)

# Load multiomic PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Function to convert Seurat ATAC assay to SummarizedExperiment for FigR
SEfromSignac <- function(obj, assayName="peaks", fetchRawCounts=TRUE) {
  message("Pulling info from assay container: ", assayName, "\n")
  
  peakRanges <- obj@assays[[assayName]]@ranges
  
  if (fetchRawCounts) {
    message("Using raw peak counts to store in SE ..\n")
    peakCounts <- obj@assays[[assayName]]@counts
  } else {
    message("WARNING: Using normalized peak counts. Set normalizeATACmat=FALSE when running FigR functions to avoid re-normalization.\n")
    peakCounts <- obj@assays[[assayName]]@data
  }
  
  cellMeta <- DataFrame(obj@meta.data)
  
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peakCounts),
    rowRanges = peakRanges,
    colData = cellMeta
  )
  
  return(SE)
}

# Convert ATAC assay to SummarizedExperiment
ATAC.se <- SEfromSignac(obj = pbmc, assayName = "peaks", fetchRawCounts = TRUE)

# Run peak-gene correlation analysis using FigR (supports parallelization)
cisCor <- runGenePeakcorr(
  ATAC.se,
  RNAmat = pbmc@assays[["SCT"]]@data,  # RNA assay matrix (normalized SCT counts)
  windowPadSize = 250000,               # +/- 250 kb window around genes
  genome = "hg38",                     # Genome build (hg38 or mm10)
  nCores = 48,                        # Number of CPU cores to use
  p.cut = NULL                        # No p-value cutoff here; will filter later
)

# Filter significant correlations by p-value (select)
cisCor.filt <- cisCor %>% dplyr::filter(pvalZ <= 0.05)

# Save results
save(cisCor, file = "FigR_PBMC_peak_gene.rda")
