rm(list=ls())  # Clear the workspace

# Load necessary libraries
library(Seurat)
library(Signac)

# Load the Seurat object (assumed to contain multi-omic PBMC data)
pbmc <- readRDS("PBMC_multiomic.rds")

# Extract peak names and replace '-' with ':' for genomic coordinate formatting (e.g., chr1:1000:1500)
peak <- pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
count <- pbmc@assays[["peaks"]]@counts
peak2 <- sub("-", ":", peak)  # Convert to colon-separated format
count@Dimnames[[1]] <- peak2  # Update row names (peak names)
counts <- as.matrix(count)

# Export the peak count matrix to a text file (tsv)
write.table(counts, file = "PBMCcount.tsv", sep = "\t")

# Extract cell type annotation (assumes active.ident holds cell type labels)
celldata <- pbmc@active.ident
celldata1 <- as.matrix(celldata)
colnames(celldata1) <- c("celltype")

# Export cell type labels
write.table(celldata1, file = "celltype.tsv", sep = "\t")

# Switch default assay to RNA for downstream export
DefaultAssay(pbmc) <- "RNA"

# If using Seurat v5, use custom seurat2anndata function (assumes you've defined or downloaded seurat2anndata.R)
source('seurat2anndata.R')
seurat2anndata(pbmc, outFile = 'pbmc.h5ad', 
               assay = "RNA", 
               main_layer = "counts", 
               transfer_layers = NULL, 
               drop_single_values = TRUE)

# If using Seurat v4, the commented-out sceasy call would work:
# library(sceasy)
# sceasy::convertFormat(pbmc, from="seurat", to="anndata",
#                      outFile='PBMC.h5ad', assay = "RNA", 
#                      main_layer = "counts", transfer_layers = NULL, 
#                      drop_single_values = TRUE)

# Extract highly variable genes from SCT assay
gene = pbmc@assays[["SCT"]]@var.features
write.table(gene, file = "PBMC_gene.csv", sep = "\t", quote = FALSE, row.names = FALSE)
