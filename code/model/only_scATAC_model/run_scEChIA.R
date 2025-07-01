rm(list=ls())  # Clear workspace

# Load required libraries
library(glassoFast)
library(nor1mix)
library(roxygen2)
library(scEChIA)
library(varbvs)
library(pracma)
library(cicero)
library(GenomicRanges)
library(rtracklayer)
library(SCEGHiC)
library(Seurat)
library(dplyr)

# Load PBMC multi-omics data
pbmc <- readRDS("PBMC_multiomic.rds")

# Annotate transcription start sites (TSS) for hg38
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Find genes present both in SCT assay and TSS annotation
genes <- intersect(pbmc@assays[["SCT"]]@counts@Dimnames[[1]], tssdata$TargetGene)

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

# Extract peak count matrix
data <- pbmc@assays[["peaks"]]@counts

# For genes with multiple promoters, sum counts into the first promoter peak and remove duplicates
for(i in seq_along(promoter)) {
  if(length(promoter[[i]]$promoters) > 1) {
    id <- which(rownames(data) %in% gsub("_", "-", promoter[[i]]$promoters))
    data[id[1], ] <- colSums(data[id, ])
    print(i)
  }
}

# Remove redundant promoter peaks (all but first promoter per gene)
idp1 <- c()
idp2 <- c()
for(i in seq_along(promoter)) {
  idp1 <- c(idp1, gsub("_", "-", promoter[[i]]$promoters)[1])
  idp2 <- c(idp2, gsub("_", "-", promoter[[i]]$promoters)[-1])
}
idp <- setdiff(idp2, idp1)
data <- data[-which(rownames(data) %in% idp), ]

# Select data for a specific cell type (e.g., "CD8 T")
celltype <- "CD8 T"
data1 <- data[, Idents(pbmc) == celltype]
data1 <- as.matrix(data1)

# Prepare peak names for liftover from hg38 to hg19
peaks <- rownames(data1)
peaks <- gsub("-", ":", peaks)
write.table(peaks, file = "peaks_for_scEChIA.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# External steps (done outside R) to extract chromosome, start, and end positions:
# cat peaks_for_scEChIA.txt | cut -d':' -f1 > scEChIA_peaks_chr.txt
# cat peaks_for_scEChIA.txt | cut -d':' -f2 > scEChIA_peaks_starts.txt
# cat peaks_for_scEChIA.txt | cut -d':' -f3 > scEChIA_peaks_ends.txt

# Read peak coordinates
chrs <- read.table(file = "scEChIA_peaks_chr.txt", sep = "\t")
starts <- read.table(file = "scEChIA_peaks_starts.txt", sep = "\t")
ends <- read.table(file = "scEChIA_peaks_ends.txt", sep = "\t")
chep <- data.frame(R1.chrom = chrs$V1, R1.start = starts$V1, R1.end = ends$V1)
head(chep)

# Liftover peaks from hg38 to hg19
f <- "hg38ToHg19.over.chain.gz"
f <- sub(".gz", "", f)
liftover.chain <- basename(f) %>% import.chain() 
basename(f) %>% unlink()
chep$line <- 1:nrow(chep)
chepr1 <- with(chep, GRanges(seqnames = Rle(R1.chrom), ranges = IRanges(start = R1.start, end = R1.end), id = line))
chepr1.19 <- unlist(liftOver(chepr1, liftover.chain))  

df1 <- data.frame(iranges = chepr1.19)
colnames(chep) <- c("R1.chrom", "R1.start", "R1.end", "iranges.id")
peakinfo <- merge(df1, chep, by = "iranges.id", all.x = TRUE)

# Filter peaks to keep only those successfully mapped and unique
peakinfo <- peakinfo[peakinfo$iranges.seqnames == peakinfo$R1.chrom, ]
peakinfo <- peakinfo[!duplicated(peakinfo$iranges.id), ]

peak_liftover <- data.frame(
  chr = peakinfo$iranges.seqnames,
  start = peakinfo$iranges.start,
  end = peakinfo$iranges.end,
  peak = paste(peakinfo$R1.chrom, peakinfo$R1.start, peakinfo$R1.end, sep = "_")
)
peak_liftover$midpoint <- (as.numeric(peak_liftover$start) + as.numeric(peak_liftover$end)) / 2

# Subset data to liftover peaks
data2 <- data1[which(rownames(data1) %in% gsub("_", "-", peak_liftover$peak)), ]
peak_liftover <- peak_liftover[which(gsub("_", "-", peak_liftover$peak) %in% rownames(data2)), ]
data3 <- cbind(peak_liftover[, 1:3], data2)

# Load chromosome sizes for hg19
chrom <- read.table("hg19.chrom.sizes.txt")

# Initialize list to hold predicted interactions per chromosome
predicted_interaction <- list()

# Loop through chromosomes 1-22 to predict peak-gene interactions
for(i in 1:22) {
  chr <- paste0("chr", i)
  
  # Load rho matrices for two reference cell types (GM12878 and IMR90)
  genomic_region1 <- read.table(paste0("scEChIA/GM12878/", chr, ".txt"))
  genomic_region1[is.na(genomic_region1)] <- 0
  
  genomic_region2 <- read.table(paste0("scEChIA/IMR90/", chr, ".txt"))
  genomic_region2[is.na(genomic_region2)] <- 0
  
  gap <- 25000  # Bin size 25kb, adjustable
  
  data4 <- data3[data3[,1] == chr, ]
  chrinfo <- data4[, 1:3]
  
  startCell <- 1
  endCell <- ncol(data4) - 3
  chromSize <- chrom[i, 2]
  
  # Calculate combined rho matrix from two cell types
  rhomatrix <- rhomatAvg(genomic_region1, genomic_region2, gap, i, data4, chrinfo)
  
  # Predict interactions using scEChIA model
  predicted_interaction[[chr]] <- Interaction_Prediction_1(chrinfo, data4, rhomatrix, i, startCell, endCell, chromSize)
  
  print(i)
}

# Save predicted interactions
save(predicted_interaction, file = "scEchIA_PBMC_peak_gene_scATAC.rda")
