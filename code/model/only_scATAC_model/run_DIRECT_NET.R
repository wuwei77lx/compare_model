rm(list=ls())
# load the required libraries
library(DIRECTNET)
library(Signac)
library(Seurat)
library(Matrix)
library(data.table)
library(xgboost)
library(ggplot2)
library(cicero)

# load data(remain ATAC data)
pbmc<-readRDS("pbmc_multiomic.rds")
DefaultAssay(pbmc) <- "peaks"
pbmc@assays$RNA <- NULL



# format genome.info
genome.info <- read.table(file = "hg38.promoter.regions.txt")
colnames(genome.info)<-c("Chrom","Starts","Ends","genes")

# select highly variable features
focus_markers=pbmc@assays[["SCT"]]@var.features

# note DIRECT-NET process 'ATAC' assay ,need 'peaks' assay
pbmc[["ATAC"]] <- pbmc[["peaks"]]

# name of the model dimension is 'wnn.umap'
pbmc@reductions$wnn.umap<-pbmc@reductions$umap 

#run model(only Seurat V4)
#pbmc <- Run_DIRECT_NET(pbmc, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap=0.5, size_factor_normalize = FALSE, genome.info = genome.info, focus_markers = focus_markers)

#if Seurat V5 and we need promoter:1000bp upstream and downstream of gene transcription start site(TSS),modify code(directnet.R')
source('directnet.R')
pbmc <- Run_DIRECT_NET(pbmc, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap=0.5, size_factor_normalize = FALSE, genome.info = genome.info, focus_markers = focus_markers)

direct.net_result <- Misc(pbmc, slot = 'direct.net')
direct.net_result <- as.data.frame(do.call(cbind,direct.net_result)) # links for markers

# check the function type name
direct.net_result$function_type <- gsub("HF","HC",direct.net_result$function_type)
direct.net_result$function_type <- gsub("Rest","MC",direct.net_result$function_type)
direct.net_result$function_type <- gsub("LF","LC",direct.net_result$function_type)

save(direct.net_result,file="DIRECTNET_pbmc_peak_gene_scATAC.rda")