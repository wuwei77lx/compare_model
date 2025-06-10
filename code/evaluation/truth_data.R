rm(list=ls())
# Load necessary libraries for downstream analysis
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)
library(reader)

# Load PBMC multiomic dataset containing scRNA-seq and scATAC-seq data
pbmc <- readRDS("PBMC_multiomic.rds")

# Annotate transcription start sites (TSS) for human genome (hg38)
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Extract peak information (chromosome, start, end) from PBMC peak assay
peak <- pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
split_peaks <- do.call("rbind", strsplit(as.character(peak), "-", fixed = TRUE))
peakinfo <- data.frame(chr = split_peaks[, 1],
                       start = as.numeric(split_peaks[, 2]),
                       end = as.numeric(split_peaks[, 3]))
peakinfo$name <- peak  # Assign peak names

# -----------------------------------
# Prepare Hi-C contact-based truth for CD8 T cells
# Hi-C data source: [ENCODE database](https://www.encodeproject.org/files/ENCFF009ONH/@@download/ENCFF009ONH.hic)
# Juicer tools command example provided in comments to extract 5kb resolution contacts for each chromosome
# wget https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar
# chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X") 
# for chromosome in "${chromosomes[@]}"; do input_hic="contact/CD8/ENCFF009ONH.hic"; output_prefix="contact/CD8/chr${chromosome}.SCALEobserved"; java -jar juicer_tools_2.13.06.jar dump observed SCALE "$input_hic" "chr$chromosome" "chr$chromosome" BP 5000 "$output_prefix"; done
# -----------------------------------

contact <- list()  # Initialize list to store contacts by chromosome

for(i in c(1:22, "X")){
  # Load contact matrix for chromosome i at 5kb resolution
  file_path <- paste0("contact/CD8/", "chr", i, ".SCALEobserved")
  data <- read.table(file_path)
  
  chr <- paste0("chr", i)
  # Filter TSS for current chromosome
  tss <- tssdata[tssdata$chr == chr, ]
  # Filter peaks for current chromosome
  enhancer <- peakinfo[peakinfo$chr == chr, ]
  
  # Calculate midpoint of each peak for bin assignment
  enhancer$midpoint <- (as.numeric(enhancer$start) + as.numeric(enhancer$end)) / 2
  enhancer$num2 <- floor(enhancer$midpoint / 5000) * 5000  # Assign peak to 5kb bin
  
  # Assign TSS to 5kb bins
  tss$num1 <- floor(tss$TargetGeneTSS / 5000) * 5000
  
  # Prepare symmetrical contact data by swapping columns
  data1 <- data[, c(2, 1, 3)]
  colnames(data) <- c("num1", "num2", "score")
  colnames(data1) <- c("num1", "num2", "score")
  
  # Combine original and swapped contact data
  data2 <- rbind(data, data1)
  
  # Merge TSS bins with contact scores
  re <- merge(tss, data2, by = "num1")
  # Merge enhancer bins with above result by enhancer bin
  re1 <- merge(enhancer, re, by = "num2")
  
  # Extract relevant columns: peak name (5), gene name (9), contact score (11)
  re2 <- re1[, c(5, 9, 11)]
  
  # Store chromosome-specific contact data
  contact[[chr]] <- re2
  print(i)
}

# Save contact list as RDS for downstream use as ground truth
saveRDS(contact, file = "truth_CD8.rds")

# -----------------------------------
# Prepare eQTL-based ground truth from GTEx Whole Blood data
# -----------------------------------

# Read significant variant-gene pairs and eGene annotation from GTEx
gtex <- read.table("GTEx/Whole_Blood.signifpairs.txt", header = TRUE)
anno <- read.table("GTEx/Whole_Blood.v8.egenes.txt", header = TRUE, fill = TRUE)
length(unique(gtex$gene_id))
length(unique(anno$gene_id))

# Merge variant-gene pairs with gene annotations by gene_id
gtex1 <- merge(gtex, anno[, c(1, 2)], by = c("gene_id"), all.x = TRUE)

# Parse variant chromosome and position from variant_id (format: chr_pos_ref_alt)
gtex2 <- data.frame(
  chr = apply(gtex1, 1, function(x) { unlist(strsplit(x[2], "_"))[1] }),
  start = as.numeric(apply(gtex1, 1, function(x) { unlist(strsplit(x[2], "_"))[2] })),
  end = as.numeric(apply(gtex1, 1, function(x) { unlist(strsplit(x[2], "_"))[2] })) + 1,
  name = gtex1$variant_id,
  gene = gtex1$gene_name
)

# Write eQTL variants with gene names to a TSV file for downstream intersection with peaks
write.table(gtex2, file = "GTEx.tsv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# -----------------------------------
# Prepare peak BED file for intersection with GTEx variants
# -----------------------------------

peak <- pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
split_peaks <- do.call("rbind", strsplit(as.character(peak), "-", fixed = TRUE))
peakinfo <- data.frame(chr = split_peaks[, 1],
                       start = as.numeric(split_peaks[, 2]),
                       end = as.numeric(split_peaks[, 3]))
peakinfo$name <- paste0("peak", 1:dim(peakinfo)[1])

# Save peaks as BED-like file without headers for bedtools usage
write.table(peakinfo, file = "macs2peak.tsv", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# -----------------------------------
# Command line notes for processing:
# mv GTEx.tsv GTEx.bed
# perl -p -i -e 's/ //g' GTEx.bed
# mv macs2peak.tsv macs2peak.bed
# perl -p -i -e 's/ //g' macs2peak.bed
# bedtools intersect -a macs2peak.bed -b GTEx.bed -wo | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$NF"\t"$9}' > GTEx1.bed

