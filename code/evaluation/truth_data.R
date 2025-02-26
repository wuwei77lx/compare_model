rm(list=ls())
# real Hi-C data generation
# load the required libraries
library(cicero)
library(Signac)
library(Seurat)
library(SCEGHiC)
library(reader)

# load data
pbmc<-readRDS("pbmc_multiomic.rds")

# load annotate transcription start sites(TSS)
tssdata <- annotateTSS("Homo sapiens", "hg38")

#peak info
peak<-pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
split_peaks <- do.call("rbind", strsplit(as.character(peak), "-", fixed = TRUE))
peakinfo <- data.frame(chr = split_peaks[, 1], start = as.numeric(split_peaks[, 2]), end = as.numeric(split_peaks[, 3]))
peakinfo$name<-peak

# the truth of CD8 T is sourced from the [ENCODE database](https://www.encodeproject.org/files/ENCFF009ONH/@@download/ENCFF009ONH.hic)
# wget https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar
# chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X") 
# for chromosome in "${chromosomes[@]}"; do input_hic="contact/CD8/ENCFF009ONH.hic"; output_prefix="contact/CD8/chr${chromosome}.SCALEobserved"; java -jar juicer_tools_2.13.06.jar dump observed SCALE "$input_hic" "chr$chromosome" "chr$chromosome" BP 5000 "$output_prefix"; done

# truth peak-gene by Hi-C
contact<-list()
for(i in c(1:22,"X")){
  file_path=paste0("contact/CD8/","chr",i,".SCALEobserved")
  data <- read.table(file_path)
  chr=paste0("chr",i)
  tss<-tssdata[tssdata$chr==chr,]
  enhancer=peakinfo[peakinfo$chr==chr,]
  enhancer$midpoint<-(as.numeric(enhancer$start)+as.numeric(enhancer$end))/2
  enhancer$num2<-floor(enhancer$midpoint/5000)*5000
  tss$num1<-floor(tss$TargetGeneTSS/5000)*5000
  data1<-data[,c(2,1,3)]
  colnames(data)<-c("num1","num2","score")
  colnames(data1)<-c("num1","num2","score")
  data2<-rbind(data,data1)
  re<-merge(tss,data2,by="num1")
  re1<-merge(enhancer,re,by="num2")
  re2<-re1[,c(5,9,11)]
  contact[[chr]]=re2
  print(i)
}

saveRDS(contact,file="truth_CD8.rds")

# real eQTL data generation
# the truth of PBMC is sourced from the [GTEx database](https://www.gtexportal.org)
gtex<-read.table("GTEx/Whole_Blood.signifpairs.txt",header=T)
anno<-read.table("GTEx/Whole_Blood.v8.egenes.txt",header=T,fill=T)
length(unique(gtex$gene_id))
length(unique(anno$gene_id))
gtex1<-merge(gtex,anno[,c(1,2)],by=c("gene_id"),all.x=T)
gtex2<-data.frame(chr=apply(gtex1,1,function(x){unlist(strsplit(x[2],"_"))[1]}),start=as.numeric(apply(gtex1,1,function(x){unlist(strsplit(x[2],"_"))[2]})),end=as.numeric(apply(gtex1,1,function(x){unlist(strsplit(x[2],"_"))[2]}))+1,name=gtex1$variant_id,gene=gtex1$gene_name)
write.table(gtex2,file = "GTEx.tsv",quote=F,col.names = F,row.names = F,sep="\t")

# peak info
peak<-pbmc@assays[["peaks"]]@counts@Dimnames[[1]]
split_peaks <- do.call("rbind", strsplit(as.character(peak), "-", fixed = TRUE))
peakinfo <- data.frame(chr = split_peaks[, 1], start = as.numeric(split_peaks[, 2]), end = as.numeric(split_peaks[, 3]))
peakinfo$name<-paste0("peak",1:dim(peakinfo)[1])
write.table(peakinfo, file = "macs2peak.tsv",row.names = FALSE ,col.names=FALSE, sep="\t",quote=FALSE)

# mv GTEx.tsv GTEx.bed
# perl -p -i -e 's/ //g' GTEx.bed
# mv macs2peak.tsv macs2peak.bed
# perl -p -i -e 's/ //g' macs2peak.bed
# bedtools intersect -a macs2peak.bed -b GTEx.bed -wo |awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$NF"\t"$9}'> GTEx1.bed

