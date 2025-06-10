rm(list=ls())  # Clear workspace

# Load required libraries
library(monocle3)
library(cicero)
library(Signac)
library(Seurat)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(magrittr)

# Load multi-omic PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Extract peak count matrix from "peaks" assay
counts <- pbmc@assays[["peaks"]]@counts

# Filter peaks with total counts > 100 to reduce noise and computational burden
s <- rowSums(counts)
index <- which(s > 100)
counts <- counts[index,]

peaks <- rownames(counts)
peaks <- gsub("-", ":", peaks)

# Write peaks to text for separate chromosome/start/end extraction (external shell commands commented)
write.table(peaks, file = "peaks_for_cicero.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# cat peaks_for_cicero.txt | cut -d':' -f1 > cicero_peaks_chr.txt
# cat peaks_for_cicero.txt | cut -d':' -f2 > cicero_peaks_starts.txt
# cat peaks_for_cicero.txt | cut -d':' -f3 > cicero_peaks_ends.txt

chrs <- read.table(file = "cicero_peaks_chr.txt", sep = "\t")
starts <- read.table(file = "cicero_peaks_starts.txt", sep = "\t")
ends <- read.table(file = "cicero_peaks_ends.txt", sep = "\t")

# Construct data frame of peak coordinates
chep <- data.frame(R1.chrom = chrs$V1, R1.start = starts$V1, R1.end = ends$V1)
head(chep)

# Convert to hg19 need a liftover
# wget "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
# gunzip hg38ToHg19.over.chain

# Liftover chain file preparation
# External download of hg38 to hg19 chain file is required and assumed done beforehand
f <- "hg38ToHg19.over.chain.gz"
f <- sub(".gz", "", f)
liftover.chain<-basename(f) %>% import.chain() # Load liftover chain file for coordinate conversion
basename(f) %>% unlink() # Remove unzipped chain file
# Add unique ID to each peak for tracking after liftover
chep$line <- 1:nrow(chep)

# Create GRanges object with hg38 coordinates for liftover
chepr1 <- with(chep, GRanges(seqnames=Rle(R1.chrom), ranges=IRanges(start=R1.start, end=R1.end), id=line))

# Perform liftover to hg19
chepr1.19 <- unlist(liftOver(chepr1, liftover.chain))

# Extract liftover results and remove duplicates (keep unique peaks)
df1 <- data.frame(iranges = chepr1.19)
index <- which(!duplicated(df1$iranges.id))
hg38_id <- df1$iranges.id[index]

# Subset counts matrix to peaks with successful liftover
counts_hg19 <- counts[hg38_id,]

# Construct new peak names based on lifted coordinates (hg19)
new_peaks <- paste0(df1$iranges.seqnames, "-", df1$iranges.start, "-", df1$iranges.end)
rownames(counts_hg19) <- new_peaks[index]

# Binarize the counts matrix (convert counts >0 to 1)
counts_hg19@x[counts_hg19@x > 0] <- 1

# Prepare peak metadata dataframe for Monocle3
peakinfo <- data.frame(chr = df1$iranges.seqnames, bp1 = df1$iranges.start, bp2 = df1$iranges.end)[index,]
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep = "_")
row.names(peakinfo) <- peakinfo$site_name

# Prepare cell metadata dataframe
cellinfo <- data.frame(cells = colnames(counts_hg19))
row.names(cellinfo) <- cellinfo$cells

# Fix counts matrix row and column names to match metadata
row.names(counts_hg19) <- row.names(peakinfo)
colnames(counts_hg19) <- row.names(cellinfo)

# Create CellDataSet (CDS) object for Cicero (Monocle3)
input_cds <- suppressWarnings(new_cell_data_set(counts_hg19, cell_metadata = cellinfo, gene_metadata = peakinfo))

# Preprocessing: detect genes (peaks), estimate size factors, and run LSI normalization
set.seed(2017)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensionality reduction with UMAP on LSI embeddings
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

# Create Cicero CDS using UMAP coordinates for co-accessibility analysis
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Load hg19 genome data from Cicero package
data("human.hg19.genome")

# Run Cicero to compute co-accessibility scores
conns_cicero <- run_cicero(cicero_cds, human.hg19.genome)

# Replace NA co-accessibility scores with zero
conns_cicero$coaccess[is.na(conns_cicero$coaccess)] <- 0

# Ensure Peak2 column is character type
conns_cicero$Peak2 <- as.character(conns_cicero$Peak2)

# Convert peak coordinates back to hg38
peaks_hg38 <- peaks[hg38_id]  # original hg38 peak names
conns_cicero_hg38 <- conns_cicero
peaks_hg19 <- rownames(counts_hg19)

# Match Cicero peaks with hg19 peak names to replace with hg38 peaks
id1 <- match(conns_cicero$Peak1, peaks_hg19)
id2 <- match(conns_cicero$Peak2, peaks_hg19)

# Check how many NA mappings (should be minimal)
length(which(is.na(id1)))
length(which(is.na(id2)))

# Replace hg19 peaks with corresponding hg38 peak names
conns_cicero_hg38$Peak1 <- peaks_hg38[id1]
conns_cicero_hg38$Peak2 <- peaks_hg38[id2]

# Replace ":" with "_" for compatibility with downstream tools
conns_cicero_hg38$Peak1 <- gsub(":", "_", conns_cicero_hg38$Peak1)
conns_cicero_hg38$Peak2 <- gsub(":", "_", conns_cicero_hg38$Peak2)

# Save Cicero peak-gene co-accessibility results
save(conns_cicero_hg38,file = "Cicero_PBMC_peak_gene_scATAC.rda")
