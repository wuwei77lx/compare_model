rm(list=ls())  # Clear workspace

# Load required libraries
library(scMultiMap)
library(Signac)
library(Seurat) 
library(GenomicRanges)
library(rtracklayer)
library(SCEGHiC)

# Load PBMC multi-omics data
pbmc <- readRDS("PBMC_multiomic.rds")

# Select highly variable genes from SCT assay
focus_markers <- pbmc@assays[["SCT"]]@var.features
# Note: make sure SCTransform was run before, else this will be empty

# Annotate transcription start sites (TSS) for hg38
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Find highly variable genes TSS annotation
genes <- intersect(focus_markers, tssdata$TargetGene)

# Get peak names, replace "-" with "_" for consistency
peak <- pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
peak <- gsub("-", "_", peak)

# Initialize lists for promoters and enhancers per gene
promoter <- list()
enhancer <- list()

# For each gene, define promoter (± 1kb) and enhancer (± 250kb) regions and find overlapping peaks
for(n in seq_along(genes)) {
  G <- genes[n]
  chr <- tssdata[which(G == tssdata$TargetGene), ]$chr
  start <- tssdata[which(G == tssdata$TargetGene), ]$TargetGeneTSS
  p1 <- paste(chr, ":", max(start-1000, 0), "-", start+1000, sep = "")
  p2 <- paste(chr, ":", max(start-250000, 0), "-", start+250000, sep = "")
  
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  enhancers <- setdiff(enhancers, promoters)  # exclude promoters from enhancers
  
  if(length(promoters) > 0) promoter[[G]] <- data.frame(gene = G, promoters = promoters)
  if(length(enhancers) > 1) enhancer[[G]] <- data.frame(gene = G, enhancers = enhancers)
  
  print(n)
}

# Candidate enhancer-gene peaks
pair_df_all<-do.call(rbind,enhancer)
colnames(pairs_df)<-c("gene","peak")
pairs_df$peak<-gsub("_","-",pairs_df$peak)

# Run scMultiMap to identify enhancer-gene peaks for each cell type.
celltype<-as.character(unique(Idents(pbmc)))

scMultiMap_results<-list()
for( cell in celltype){
  dataset<-subset(pbmc, idents = cell)
  scMultiMap_res <- scMultiMap(dataset, pairs_df,
                               gene_assay = 'RNA', # name of the gene assay
                               peak_assay = 'peaks')
  scMultiMap_results[[cell]]<-scMultiMap_res
}

# Save results
save(scMultiMap_results,file="scMultiMap_peak_gene.rda"")