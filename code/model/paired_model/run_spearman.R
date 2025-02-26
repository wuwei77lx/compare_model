rm(list=ls())
# load the required libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)

# load data
pbmc<-readRDS("pbmc_multiomic.rds")

# data preprocessing(aggregate)-reference SCEGHiC
agg.data <- process_data(pbmc, max_overlap = 0.5)
rna<-agg.data[["rna"]]
atac<-agg.data[["atac"]]
atac<-atac[rowSums(atac)!=0,]
rownames(atac) <- gsub("-", "_", rownames(atac))

# load annotate transcription start sites(TSS)
tssdata <- annotateTSS("Homo sapiens", "hg38")

# select highly variable features
gene<-pbmc@assays[["SCT"]]@var.features
gene1<-intersect(gene,tssdata$TargetGene)
rna1<-rna[which(rownames(rna)%in%gene1),]
rna2<-rna1[apply(rna1,1,function(x){sum(x!=0)>2}),]

# peak info 
peak<-pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
peak<-gsub("-","_",peak)

#spearman
results<-list()
for(n in 1:dim(rna2)[1]){
  G=rownames(rna2)[n]
  y=rna2[n,]
  chr<-tssdata[which(G==tssdata$TargetGene),]$chr
  start<-tssdata[which(G==tssdata$TargetGene),]$TargetGeneTSS
  p1 <- paste(chr,":",max(start-1000,0),"-",start+1000,sep = "")
  p2 <- paste(chr,":",max(start-250000,0),"-",start+250000,sep = "")
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  enhancers <- setdiff(enhancers,promoters)
  if(length(enhancers)>0){
    x=atac[which(rownames(atac)%in%enhancers),]
    x2<-rbind(y,x)
    rownames(x2)<-c("y",enhancers)
    dd<-cor(t(x2),method="spearman")
    data<-data.frame(gene=G,peak2=enhancers,score=dd[,1][-1],method="spearman")
    results[[n]]=data
  }else{
    results[[n]]=NULL
  }
  print(n)
}
results<-do.call(rbind,results)

save(results,file="spearman_pbmc_peak_gene.rda")