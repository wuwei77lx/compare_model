rm(list=ls())  # Clear workspace

# Load required libraries
library(Seurat)
library(Signac)
library(ArchR)

# Load multi-omics PBMC dataset
pbmc <- readRDS("PBMC_multiomic.rds")

# Set seed for reproducibility, configure ArchR threads, and specify reference genome
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")

# Create Arrow files (binary format) with a non-binarized TileMatrix for downstream analysis
ArrowFiles <- createArrowFiles(
  inputFiles = 'pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz',
  sampleNames = 'PBMC',
  minTSS = 0, #Dont set this too high because you can always increase later
  maxFrags = Inf,
  minFrags = 0, 
  minFragSize = 0,
  maxFragSize = Inf,
  addTileMat = TRUE,
  TileMatParams=list(binarize = FALSE),
  addGeneScoreMat = TRUE,
  force=T
)

# Create an ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "SCARlink/pbmc_archR",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

# Prefix Seurat object cell names to match ArchR cell names
colnames(pbmc) <- paste0("PBMC#", colnames(pbmc))
colnames(pbmc)[1:10]

# Prepare Seurat object for integration: copy SCT assay to RNA slot, assign cell type identities, set default assay, and save object
pbmc[["RNA"]]=pbmc[["SCT"]]
pbmc$celltype<-as.character(Idents(pbmc))
DefaultAssay(pbmc) <- "RNA"
saveRDS(pbmc,file="SCARlink/pbmc.rds")

# Find cells present in both ArchR and Seurat objects
common_cells <- intersect(proj$cellNames, colnames(pbmc))
length(common_cells)

# Subset ArchR project to include only cells present in Seurat object
proj <- subsetArchRProject(proj, cells = common_cells, outputDirectory = "SCARlink/pbmc_Seurat_to_ArchR", dropCells = TRUE)

# Perform iterative latent semantic indexing (LSI) on scATAC-seq TileMatrix
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# Save the final ArchR project
proj <- saveArchRProject(ArchRProj = proj)