rm(list=ls())  # Clear workspace

# Load required libraries
require(GenomicFeatures)
library(ChIPseeker)
library(Signac)
library(Seurat)
library(SCEGHiC)
library(cicero)

# Load multi-omic PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Create transcript database using UCSC hg38 refGene table
hg38.refseq.db <- makeTxDbFromUCSC(genome="hg38", table="refGene") 
# Alternatively: makeTxDbFromGFF("hg38.refGene.gtf")

# Extract peak names from the ATAC assay
peak <- pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
peak1 <- matrix(peak, length(peak))

# Split peak names into chromosome, start, end
col1 <- apply(peak1, 1, function(x) strsplit(x, "-")[[1]][1])
col2 <- apply(peak1, 1, function(x) strsplit(x, "-")[[1]][2])
col3 <- apply(peak1, 1, function(x) strsplit(x, "-")[[1]][3])
peak2 <- cbind(col1, col2, col3)

# Convert to GRanges object for annotation
peak3 <- GRanges(
  seqnames = Rle(peak2[,1]),
  ranges = IRanges(start = as.numeric(peak2[,2]), end = as.numeric(peak2[,3]))
)

# Annotate peaks using ChIPseeker
txdb <- hg38.refseq.db
peaksAnno <- annotatePeak(peak3, tssRegion = c(-1000, 1000), TxDb = txdb)
peaksan <- as.data.frame(peaksAnno)

# Format peak names
peaksan$name <- paste(peaksan$seqnames, peaksan$start, peaksan$end, sep = "_")

# Simplify annotation: only retain "Promoter" and "Distal Intergenic", others are "others"
peaksan$annotation[!(peaksan$annotation %in% c("Distal Intergenic", "Promoter"))] <- "others"

# Annotate transcription start sites (TSS) for all genes
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Find genes common to RNA assay and TSS annotation
genes <- intersect(pbmc@assays[["SCT"]]@counts@Dimnames[[1]], tssdata$TargetGene)

# Replace "-" with "_" in peak names for consistency
peak <- gsub("-", "_", peak)

# Initialize empty lists to store promoter and enhancer peak mappings
promoter <- list()
enhancer <- list()

# Loop over genes to define promoter (±1 kb) and enhancer (±250 kb) windows
for(n in seq_along(genes)) {
  G <- genes[n]
  chr <- tssdata[which(G == tssdata$TargetGene), ]$chr
  start <- tssdata[which(G == tssdata$TargetGene), ]$TargetGeneTSS
  
  # Define promoter and enhancer regions as genomic intervals
  p1 <- paste(chr, ":", max(start - 1000, 0), "-", start + 1000, sep = "")
  p2 <- paste(chr, ":", max(start - 250000, 0), "-", start + 250000, sep = "")
  
  # Find overlapping peaks within these regions
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  
  # Remove any promoter peaks from the enhancer set
  enhancers <- setdiff(enhancers, promoters)
  
  # Store results only if there are overlaps
  if (length(promoters) > 0) promoter[[G]] <- data.frame(gene = G, promoters = promoters)
  if (length(enhancers) > 1) enhancer[[G]] <- data.frame(gene = G, enhancers = enhancers)
  
  print(n)  # Progress tracker
}

# Combine all promoter entries into one dataframe
promoter <- do.call(rbind, promoter)
colnames(promoter) <- c("gene", "name")

# Compute pseudobulk expression for RNA and ATAC assays grouped by cell identity
pseudobulk_RNA <- AverageExpression(
  pbmc,
  assays = "SCT",
  return.seurat = FALSE,
  group.by = "ident",
  slot = "data"  # assumes log-normalized data; converts back to raw scale
)

pseudobulk_ATAC <- AverageExpression(
  pbmc,
  assays = "peaks",
  return.seurat = FALSE,
  group.by = "ident",
  slot = "data"
)

# Convert to data frames
atac <- as.data.frame(pseudobulk_ATAC[["peaks"]])
rna <- as.data.frame(pseudobulk_RNA[["SCT"]])
rna$gene <- rownames(rna)
atac$name <- gsub("-", "_", rownames(atac))

# Merge promoter regions with ATAC signals
promoter1 <- merge(promoter, atac, by = "name", all.x = TRUE)
promoter1 <- promoter1[, -1]  # Drop 'name' column

# Aggregate multiple promoter peaks per gene by averaging signal
promoter2 <- promoter1 %>%
  group_by(gene) %>%
  summarise_all(mean)

# Prepare TSS annotation table
colnames(tssdata) <- c("chr", "gene", "tss")
tss <- tssdata

# Merge TSS, promoter activity, and gene expression
tss1 <- merge(tss, promoter2, by = "gene")
tss1[is.na(tss1)] <- 0  # Replace NA with 0
tss2 <- merge(tss1, rna, by = "gene", all.x = TRUE)

# Merge peak activity with functional annotations
atac1 <- merge(atac, peaksan, by = "name", all.x = TRUE)

# Extract CD4 T cell-specific promoter and enhancer information
Genelist <- data.frame(
  chr = tss2$chr,
  symbol = tss2$gene,
  tss = tss2$tss,
  Expression = tss2$`CD4 T.y`,
  PromoterActivityQuantile = tss2$`CD4 T.x`
)

enhancerlist <- data.frame(
  chr = atac1$seqnames,
  start = atac1$start,
  end = atac1$end,
  name = atac1$name,
  class = atac1$annotation,
  activity_base = atac1$`CD4 T`
)

# Save output files for downstream use
write.table(enhancerlist, file = "EnhancerList_CD4.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Genelist, file = "GeneList_CD4.txt", sep = "\t", row.names = FALSE, quote = FALSE)