rm(list=ls())
# load the required libraries
library(Seurat)
library(Signac)

#  load data
pbmc<-readRDS("pbmc_multiomic.rds")

# format peak count
peak<-pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
count<-pbmc@assays[["peaks"]]@counts
peak2<-sub("-",":",peak)
count@Dimnames[[1]]=peak2
counts<-as.matrix(count)
write.table(counts, file = "pbmccount.tsv",sep="\t")

# format celltype info
celldata<-pbmc@active.ident
celldata1<-as.matrix(celldata)
colnames(celldata1)<-c("celltype")
write.table(celldata1, file = "celltype.tsv",sep="\t")

# convert a Seurat object to an Anndata object(scRNA-seq)
DefaultAssay(pbmc) <- "RNA"

# use R sceasy ,but only Seurat(V4)
#sceasy::convertFormat(pbmc, from="seurat", to="anndata",
#                      outFile='pbmc.h5ad',assay = "RNA", main_layer = "counts", transfer_layers = NULL, drop_single_values = TRUE)

# if Seurat(V5),modify code("seurat2anndata")
source('seurat2anndata.R')
seurat2anndata(pbmc, outFile = 'pbmc.h5ad', assay = "RNA", main_layer = "counts", transfer_layers = NULL, drop_single_values = TRUE)




# select highly variable features
gene=pbmc@assays[["SCT"]]@var.features
write.table(gene,file="pbmc_gene.csv",sep="\t",quote=F,row.names=F)