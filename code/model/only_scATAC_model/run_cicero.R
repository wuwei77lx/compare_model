rm(list=ls())
# load the required libraries
library(monocle3)
library(cicero)
library(Signac)
library(Seurat)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(magrittr)
# load data
pbmc<-readRDS("pbmc_multiomic.rds")
counts<-pbmc@assays[["peaks"]]@counts
s <- rowSums(counts)
index <- which(s > 100)
counts <- counts[index,]

# convert to hg19 first
peaks <- rownames(counts)
peaks <- gsub("-",":",peaks)
write.table(peaks, file = "peaks_for_cicero.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
# cat peaks_for_cicero.txt | cut -d':' -f1 > cicero_peaks_chr.txt
# cat peaks_for_cicero.txt | cut -d':' -f2 > cicero_peaks_starts.txt
# cat peaks_for_cicero.txt | cut -d':' -f3 > cicero_peaks_ends.txt

chrs <- read.table(file = "cicero_peaks_chr.txt",sep = "\t")
starts <- read.table(file = "cicero_peaks_starts.txt",sep = "\t")
ends <- read.table(file = "cicero_peaks_ends.txt",sep = "\t")
chep <- data.frame(R1.chrom = chrs$V1,R1.start = starts$V1,R1.end = ends$V1)
head(chep)

# convert to hg19
# will need a liftover
# wget "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
# gunzip hg38ToHg19.over.chain

f <- "hg38ToHg19.over.chain.gz"
f <- sub(".gz","",f)
liftover.chain<-basename(f) %>% import.chain() 
basename(f) %>% unlink()
chep$line<-1:nrow(chep)
chepr1<-with(chep,GRanges(seqnames=Rle(R1.chrom),ranges=IRanges(start=R1.start,end=R1.end),id=line))
chepr1.19<-unlist(liftOver(chepr1,liftover.chain))  
df1 <-  data.frame(iranges = chepr1.19)
index <- which(!duplicated(df1$iranges.id))
hg38_id <- df1$iranges.id[index]
counts_hg19 <- counts[hg38_id,]
new_peaks <- paste0(df1$iranges.seqnames,"-",df1$iranges.start,"-",df1$iranges.end)
rownames(counts_hg19) <- new_peaks[index]

# binarize the matrix
counts_hg19@x[counts_hg19@x > 0] <- 1

# format peak info
peakinfo<-data.frame(chr=df1$iranges.seqnames,bp1=df1$iranges.start,bp2=df1$iranges.end)[index,]
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

# format cell info
cellinfo<-data.frame(colnames(counts_hg19))
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"
row.names(counts_hg19) <- row.names(peakinfo)
colnames(counts_hg19) <- row.names(cellinfo)

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(counts_hg19,cell_metadata = cellinfo,gene_metadata = peakinfo))

# data preprocessing
set.seed(2017)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# load reference genome information
data("human.hg19.genome")

# running Cicero
conns_cicero <- run_cicero(cicero_cds,human.hg19.genome)
conns_cicero$coaccess[is.na(conns_cicero$coaccess)] <- 0
conns_cicero$Peak2 <- as.character(conns_cicero$Peak2)

# convert to hg38
peaks_hg38 <- peaks[hg38_id]
conns_cicero_hg38 <- conns_cicero
peaks_hg19 <- rownames(counts_hg19)
id1 <- match(conns_cicero$Peak1,peaks_hg19)
f1 <- which(is.na(id1))
length(f1)
id2 <- match(conns_cicero$Peak2,peaks_hg19)
f2 <- which(is.na(id2))
length(f2)
conns_cicero_hg38$Peak1 <- peaks_hg38[id1]
conns_cicero_hg38$Peak2 <- peaks_hg38[id2]
conns_cicero_hg38$Peak1<-gsub(":","_",conns_cicero_hg38$Peak1)
conns_cicero_hg38$Peak2<-gsub(":","_",conns_cicero_hg38$Peak2)
save(conns_cicero_hg38,file = "cicero_pbmc_peak_gene_scATAC.rda")
