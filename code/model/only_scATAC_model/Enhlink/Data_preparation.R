rm(list=ls())  # Clear workspace

# Load required libraries
library(Signac)
library(Seurat)
library(Matrix)

# Load multi-omics PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Select highly variable genes from SCT assay
focus_markers <- pbmc@assays[["SCT"]]@var.features
# Note: make sure SCTransform was run before, else this will be empty

# Write the sparse matrix in the MTX matrix format
writeMM(obj = pbmc@assays$peaks@counts, file = "pbmc.ATAC.mtx")

# Write the cell ID tags (the columns)
write.table(colnames(pbmc@assays$peaks@counts),
            file = "pbmc.ATAC.colnames", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write the chromosome peak coordinates (the rows)
# We need to format the peak coordinates to be tab separated (chrID, start, stop)
write.table(gsub("-", "\t", rownames(pbmc@assays$peaks@counts)),
            file = "pbmc.ATAC.rownames",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# Select target genes whose promoters are defined as ±1 kb around the transcription start site (TSS) and save them for downstream analysis
geneinfo<-read.table("Homo_sapiens.GRCh38.99.TSS.1K.bed")
gene<-intersect(focus_markers,geneinfo[,4])
write.table(gene,
            file = "target_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Extract cell IDs and cluster labels from the Seurat object, standardize cluster names, and save them as a CSV file for downstream analysis
clusters<-data.frame(celID=colnames(pbmc),clusterID=as.character(Idents(pbmc)))
clusters$clusterID <- gsub("[ -]", "_", clusters$clusterID)
write.csv(clusters,
          file = "clusters.csv",
          quote = FALSE, row.names = FALSE)

# If batch information is available, create a covariate table for cells to account for batch effects in downstream analyses
# cov<-data.frame(celID=colnames(pbmc),batch=as.character(pbmc$batch))
# write.csv(cov,file = "cov.csv",quote = FALSE, row.names = FALSE)