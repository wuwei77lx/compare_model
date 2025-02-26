rm(list=ls())
# load the required libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)

# load data
pbmc<-readRDS("pbmc_multiomic.rds")
celltype="CD8 T"
dataset<-subset(pbmc, idents = celltype)

#remain ATAC data
atac=dataset@assays[["peaks"]]@counts

# binarize the matrix
atac@x[atac@x > 0] <- 1

# load annotate transcription start sites(TSS)
tssdata <- annotateTSS("Homo sapiens", "hg38")

# select highly variable features
focus_markers=pbmc@assays[["SCT"]]@var.features
focus_markers<-intersect(focus_markers,tssdata$TargetGene)

# peak info
peak=rownames(atac)
peak=gsub("-","_",peak)
rownames(atac)=peak

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
  if(length(enhancers)>1&length(promoters)>0){
    if(length(promoters)>1){
      y=colSums(atac[which(rownames(atac)%in%promoters),])
      y[y>0]=1
    }else{
      y=atac[which(rownames(atac)%in%promoters),]
    }
    if(sum(y)>0){
      if(length(enhancers)==1){
        observed=table(y,x)
        result <- chisq.test(observed,simulate.p.value = TRUE)
        results[[n]]=data.frame(gene=G,peak1=promoters[1],peak2=enhancers,p_value=result$p.value,fdr=result$p.value)
      }else{
        res=data.frame(gene=G,peak1=promoters[1],peak2=enhancers,p_value=0)
        for(m in 1:length(enhancers)){
          z=x[m,]
          observed=table(y,z)
          result <- chisq.test(observed,simulate.p.value = TRUE)
          res[m,4]<-result$p.value
        }
        res$fdr<-p.adjust(res$p_value, method = "BH")
        results[[n]]=res
      } 
    }else{
      results[[n]]=NULL
    }
  }
  else{
    results[[n]]=NULL
    
  }
  print(n)
}
results<-do.call(rbind,results)

save(results,file="chi2_fdr_pbmc_peak_gene_scATAC.rda")