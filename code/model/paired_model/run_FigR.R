rm(list=ls())
# load the required libraries
library(Seurat)
library(Signac)
library(FigR)

#  load data
pbmc<-readRDS("pbmc_multiomic.rds")

# construct the ATAC data format for FigR
SEfromSignac <- function(obj, # Seurat object containing assay that houses peak (ATAC) info
                         assayName="peaks", # Assay name for peak (ATAC) info
                         fetchRawCounts=TRUE){ # Whether to use raw counts (Default), otherwise fetch normalized counts if present 
  message("Pulling info from assay container: ",assayName,"\n")
  
  peakRanges <- obj@assays[[assayName]]@ranges
  
  if(fetchRawCounts){
    message("Using raw peak counts to store in SE ..\n")
    peakCounts <- obj@assays[[assayName]]@counts
  } else {
    message("WARNING: Pulling from normalized peak count slot.\n If running FigR's gene-peak association testing, make sure you set normalizeATACmat to FALSE in the runGenePeakcorr and getDORCscores functions to avoid renormalizing data internally")
    peakCounts <- obj@assays[[assayName]]@data
  }
  
  
  cellMeta <- DataFrame(obj@meta.data)
  
  
  SE <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=peakCounts),
                                                   rowRanges = peakRanges,
                                                   colData = cellMeta)
  return(SE)
}
ATAC.se <- SEfromSignac(obj=pbmc,assayName = "peaks",fetchRawCounts = TRUE)

# run using multiple cores if parallel support,predict peak-gene
cisCor <- runGenePeakcorr(ATAC.se,
                          RNAmat = pbmc@assays[["SCT"]]@data, #if Seurat V5 ,RNAmat = object@assays[[rna_assay]]@layers[["counts"]]
                          windowPadSize=250000,
                          genome = "hg38", # Also supports mm10 and hg38
                          nCores = 48,
                          p.cut=NULL)

# filter peak-gene correlations by p-value (select)
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.05)

save(cisCor,file="FigR_pbmc_peak_gene.rda")