rm(list=ls())  # Clear workspace

# Load required packages
library(DIRECTNET)
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(xgboost)
library(ggplot2)
library(cicero)

# Load data
pbmc <- readRDS("PBMC_multiomic.rds")

# Load promoter regions annotation, expect 4 columns: chrom, start, end, gene
genome.info <- read.table(file = "hg38.promoter.regions.txt")
colnames(genome.info) <- c("Chrom", "Starts", "Ends", "genes")

# Select highly variable genes from SCT assay
focus_markers <- pbmc@assays[["SCT"]]@var.features
# Note: make sure SCTransform was run before, else this will be empty

# DIRECT-NET expects ATAC assay; copy peaks assay to ATAC
pbmc[["ATAC"]] <- pbmc[["peaks"]]
# Confirm that "peaks" assay exists and is formatted correctly

# Rename UMAP reduction to wnn.umap as expected by DIRECT-NET
pbmc@reductions$wnn.umap <- pbmc@reductions$umap
# Make sure 'umap' is WNN-based UMAP; otherwise, this might cause issues

# Run model (only Seurat V4)
# pbmc <- Run_DIRECT_NET(pbmc, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap=0.5, size_factor_normalize = TRUE, genome.info = genome.info, focus_markers = focus_markers)

# Source the modified directnet.R (for Seurat v5 compatibility)
source('directnet.R')

# Run DIRECT-NET
pbmc <- Run_DIRECT_NET(pbmc, peakcalling=FALSE, k_neigh=50, atacbinary=TRUE,
                       max_overlap=0.5, size_factor_normalize=TRUE,
                       genome.info=genome.info, focus_markers=focus_markers)

# Extract direct.net results from pbmc Misc slot
direct.net_result <- Misc(pbmc, slot='direct.net')

# Combine list elements column-wise into one data.frame
direct.net_result <- as.data.frame(do.call(cbind, direct.net_result))
# Warning: ensure all elements have compatible row counts and types

# Replace strings in function_type column
if ("function_type" %in% colnames(direct.net_result)) {
  direct.net_result$function_type <- gsub("HF","HC", direct.net_result$function_type)
  direct.net_result$function_type <- gsub("Rest","MC", direct.net_result$function_type)
  direct.net_result$function_type <- gsub("LF","LC", direct.net_result$function_type)
}

# Save the result for future use
save(direct.net_result, file="DIRECTNET_PBMC_peak_gene.rda")
