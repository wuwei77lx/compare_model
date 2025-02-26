rm(list=ls())
# load the required libraries
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
#  load data
pbmc<-readRDS("pbmc_multiomic.rds")

DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance=250000,
  pvalue_cutoff = 1.01,
  score_cutoff = 0
)

signac<-pbmc@assays[["peaks"]]@links

save(signac,file="signac_pbmc_peak_gene.rda")