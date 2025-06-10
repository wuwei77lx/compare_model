rm(list=ls())

# Load required libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)

# Load multiomic PBMC dataset (remain ATAC data)
pbmc <- readRDS("PBMC_multiomic.rds")
DefaultAssay(pbmc) <- "peaks"
pbmc@assays$RNA <- NULL
pbmc<-RunUMAP(pbmc, reduction = 'lsi', dims = 2:40)

# Data preprocessing - aggregation for SCEGHiC reference
agg.data <- process_data(pbmc, max_overlap = 0.5)

# Remain ATAC data
atac<-agg.data[["atac"]]
atac<-atac[rowSums(atac)!=0,]
rownames(atac) <- gsub("-", "_", rownames(atac))

# Load transcription start sites (TSS) annotation for human hg38
tssdata <- annotateTSS("Homo sapiens", "hg38")

# Select highly variable features
focus_markers=pbmc@assays[["SCT"]]@var.features
focus_markers<-intersect(focus_markers,tssdata$TargetGene)

# Peak info: convert "-" to "_" to match naming convention 
peak<-pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
peak<-gsub("-","_",peak)

# Calculate Pearson correlations between gene expression and peaks in ±250 kb window
results<-list()
for(n in 1:length(focus_markers)){
  G=focus_markers[n]
  chr<-tssdata[which(G==tssdata$TargetGene),]$chr
  start<-tssdata[which(G==tssdata$TargetGene),]$TargetGeneTSS
  p1 <- paste(chr,":",max(start-1000,0),"-",start+1000,sep = "")
  p2 <- paste(chr,":",max(start-250000,0),"-",start+250000,sep = "")
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  enhancers <- setdiff(enhancers,promoters)
  x=atac[which(rownames(atac)%in%enhancers),]
  if(length(enhancers)>0&length(promoters)>0){
    if(length(promoters)>1){
      y=colSums(atac[which(rownames(atac)%in%promoters),])
    }else{
      y=atac[which(rownames(atac)%in%promoters),]
    }
    x2<-rbind(y,x)
    rownames(x2)<-c(promoters[1],enhancers)
    dd<-cor(t(x2))
    data<-data.frame(gene=G,peak1=promoters[1],peak2=enhancers,score=dd[,1][-1],method="pearson")
    results[[n]]=data
  }
  else{
    results[[n]]=NULL
    
  }
  print(n)
}
# Combine results for all genes
results <- do.call(rbind, results)

# Save results
save(results,file="Pearson_PBMC_peak_gene_scATAC.rda")
