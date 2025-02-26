rm(list=ls())
# load the required libraries
library(glassoFast)
library(nor1mix)
library(roxygen2)
library(scEChIA)
library(varbvs)
library(pracma)
library(scEChIA)
library(cicero)
library(GenomicRanges)
library(rtracklayer)
library(SCEGHiC)
library(Seurat)
library(dplyr)

#  load data
pbmc<-readRDS("pbmc_multiomic.rds")

# load annotate transcription start sites(TSS)
tssdata <- annotateTSS("Homo sapiens", "hg38")

# enhancer and promoter of each gene
genes<-intersect(pbmc@assays[["SCT"]]@counts@Dimnames[[1]],tssdata$TargetGene)
peak<-pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
peak<-gsub("-","_",peak)
promoter<-list()
enhancer<-list()
for(n in 1:length(genes)){
  G=genes[n]
  chr<-tssdata[which(G==tssdata$TargetGene),]$chr
  start<-tssdata[which(G==tssdata$TargetGene),]$TargetGeneTSS
  p1 <- paste(chr,":",max(start-1000,0),"-",start+1000,sep = "")
  p2 <- paste(chr,":",max(start-250000,0),"-",start+250000,sep = "")
  promoters <- find_overlapping_coordinates(peak, p1)
  enhancers <- find_overlapping_coordinates(peak, p2)
  enhancers <- setdiff(enhancers,promoters)
  if(length(promoters)>0){
    promoter[[G]]=data.frame(gene=G,promoters=promoters)
  }
  if(length(enhancers)>1){
    enhancer[[G]]=data.frame(gene=G,enhancers=enhancers)
  }
  print(n)
}

# create scEChIA data
data=pbmc@assays[["peaks"]]@counts

# retain only one promoter for each gene, summing the expression levels of multiple promoters
for(i in 1:length(promoter)){
  if(length(promoter[[i]]$promoters)>1){
    id<-which(rownames(data)%in%gsub("_","-",promoter[[i]]$promoters))
    data[id[1],]<-colSums(data[id,])
    print(i)
  }
}
idp1<-c()
idp2<-c()
for(i in 1:length(promoter)){
  idp1<-c(idp1,gsub("_","-",promoter[[i]]$promoters)[1])
  idp2<-c(idp2,gsub("_","-",promoter[[i]]$promoters)[-1])
}
idp<-setdiff(idp2,idp1)
data<-data[-which(rownames(data)%in%idp),]

# select celltype
celltype="CD8 T"
data1=data[,Idents(pbmc)==celltype]
data1<-as.matrix(data1)

# convert to hg19 first
peaks <- rownames(data1)
peaks <- gsub("-",":",peaks)
write.table(peaks, file = "peaks_for_scEChIA.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
# cat peaks_for_cicero.txt | cut -d':' -f1 > scEChIA_peaks_chr.txt
# cat peaks_for_cicero.txt | cut -d':' -f2 > scEChIA_peaks_starts.txt
# cat peaks_for_cicero.txt | cut -d':' -f3 > scEChIA_peaks_ends.txt

chrs <- read.table(file = "scEChIA_peaks_chr.txt",sep = "\t")
starts <- read.table(file = "scEChIA_peaks_starts.txt",sep = "\t")
ends <- read.table(file = "scEChIA_peaks_ends.txt",sep = "\t")
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
colnames(chep) <-c("R1.chrom","R1.start","R1.end","iranges.id")
peakinfo<-merge(df1,chep,by="iranges.id",all.x=T)
peakinfo<-peakinfo[peakinfo$iranges.seqnames==peakinfo$R1.chrom,]
peakinfo<-peakinfo[which(!duplicated(peakinfo$iranges.id)),]
peak_liftover<-data.frame(chr=peakinfo$iranges.seqnames,start=peakinfo$iranges.start,end=peakinfo$iranges.end,peak=paste(peakinfo$R1.chrom,peakinfo$R1.start,peakinfo$R1.end,sep="_"))
peak_liftover$midpoint<-(as.numeric(peak_liftover$start)+as.numeric(peak_liftover$end))/2


data2<-data1[which(rownames(data1)%in%gsub("_","-",peak_liftover$peak)),]
peak_liftover<-peak_liftover[which(gsub("_","-",peak_liftover$peak)%in%rownames(data2)),]
data3<-cbind(peak_liftover[,1:3],data2)

# run model(predict interaction with the help of two different cell types' rho value)
chrom<-read.table("hg19.chrom.sizes.txt")
predicted_interaction=list()
for(i in c(1:22)){
  chr=paste0("chr",i)
  genomic_region1 = read.table(paste0("scEChIA/GM12878/",chr,".txt"))
  genomic_region1[is.na(genomic_region1)] = 0
  genomic_region2 = read.table(paste0("scEChIA/IMR90/",chr,".txt"))
  genomic_region2[is.na(genomic_region2)] = 0
  gap = 25000 # bin size 25kb (User can change it as per his choice)
  patternf = i # chromosome number 
  chrNo = patternf
  data4<-data3[data3[,1]==chr,]
  chrinfo = data4[, 1:3] # chromosome location (chr, start, end)
  startCell = 1 # start sample column
  endCell = dim(data4)[2]-3# end sample
  chromSize = chrom[i,2]  # size of chromosome 
  rhomatrix = rhomatAvg(genomic_region1, genomic_region2, gap, patternf, data4, chrinfo) # calculate rho value
  predicted_interaction[[chr]] = Interaction_Prediction_1(chrinfo, data4, rhomatrix, chrNo, startCell, endCell, chromSize)  #predict interaction
  print(i)
}

save(predicted_interaction,file="scEchIA_pbmc_peak_gene_scATAC.rda")