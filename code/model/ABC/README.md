# Activity-by-contact (ABC) model

## Introduction
This directory primarily contains two components:

1. **Running the ABC model**: Scripts and resources for applying the ABC model to predict enhancer–gene interactions.

2. **Mouse average Hi-C data generation**: Step-by-step scripts for processing multiple mouse Hi-C datasets and generating a bulk average Hi-C data.


## Running the ABC model

**Running the ABC model consists of the following steps**:

 1. Generate required input files

 2. Computing the ABC Score

### Step 1. Generate required input files

Prepare ```EnhancerList.txt``` and ```GeneList.txt``` from paired scRNA-seq and scATAC-seq data. We first summarized gene expression and chromatin accessibility for each cell type. Peaks were annotated using ChIPseeker and categorized as “Promoter,” “Distal Intergenic,” or “Others” based on their genomic locations relative to the TSS. See the code in ```|-- code/model/ABC/Data_preparation```.

Main output files:

  * **EnhancerList.txt**: Candidate enhancer regions with summarized ATAC read counts
  * **GeneList.txt**: Summarized RNA read counts on gene bodies and gene promoter regions

### Step 2. Computing the ABC Score
Use the provided scripts to compute the ABC score for each enhancer–gene pair by integrating enhancer activity (such as accessibility or expression) with bulk average Hi-C data.

**Sample Command**:

```
python predict.py \
--enhancers /picb/bigdata/project/liangxuan/data/human/pbmc/ABC/EnhancerList_CD4.txt \
--genes /picb/bigdata/project/liangxuan/data/human/pbmc/ABC/GeneList_CD4.txt \
--hic_file /picb/bigdata/project/liangxuan/data/human_contact/AvgHiC \
--chrom_sizes /picb/bigdata/project/liangxuan/data/human/hg38.chrom.sizes.txt \
--hic_resolution 5000 \
--window 250000 \
--expression_cutoff 0.0  \
--tss_slop 1000 \
--promoter_activity_quantile_cutoff 0.0 \
--cellType average \
--outdir /picb/bigdata/project/liangxuan/data/human/pbmc/ABC/CD4/ \
--hic_type avg \
--hic_gamma 0.87 \
--hic_scale 5.41 \
--make_all_putative

```

The main output file is:
* **EnhancerPredictionsAllPutative.tsv.gz**: all element-gene pairs with scores above the provided threshold.

## Mouse average Hi-C data generation

**Generating the average Hi-C consists of the following steps**:

1. Download the Hi-C matrix files from Juicebox (in .hic format)

2. Fit each Hi-C dataset to a power-law model to extract the decay parameters

3. Generate the average Hi-C

### Step 1. Download .hic matrix files for all selected cell types from Juicebox
  
```
wget https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar

# hic_file can use local files to save time
python code/model/ABC/juicebox_dump.py \
--hic_file https://www.encodeproject.org/files/ENCFF409ZRS/@@download/ENCFF409ZRS.hic \
--juicebox "java -jar juicer_tools_2.13.06.jar" \
--resolution 5000 \
--outdir mouse_contact/ESC1
```

### Step 2. Fit Hi-C contact frequency data to a power-law model and extract the corresponding parameters

```
python code/model/ABC/compute_powerlaw_fit_from_hic.py \
--hic_dir mouse_contact/ESC1 \
--outDir mouse_contact/ESC1/powerlaw \
--hic_resolution 5000
```


### Step 3. Generate the average Hi-C

```
# chromosome = ['chr' + str(x) for x in range(1,20)] + ['chrX']
python code/model/ABC/makeAverageHiC.py \
--celltypes ESC1,ESC2,CH12F3,B,CH12LX,fiber,epi \
--chromosome chr1 \
--basedir mouse_contact \
--outDir mouse_contact/average_HiC
```

**The mouse average Hi-C data (~13.9 GB) can be downloaded from**: [10.5281/zenodo.14849886](https://zenodo.org/record/14849886).

## Reference

Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

Nasser, J., Bergman, D. T., Fulco, C. P., Guckelberger, P., Doughty, B. R., Patwardhan, T. A., Jones, T. R., Nguyen, T. H., Ulirsch, J. C., Lekschas, F., Mualim, K., Natri, H. M., Weeks, E. M., Munson, G., Kane, M., Kang, H. Y., Cui, A., Ray, J. P., Eisenhaure, T. M., Collins, R. L., … Engreitz, J. M. (2021). Genome-wide enhancer maps link risk variants to disease genes. Nature, 593(7858), 238–243. https://doi.org/10.1038/s41586-021-03446-x

Yu, G., Wang, L. G., & He, Q. Y. (2015). ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics (Oxford, England), 31(14), 2382–2383. https://doi.org/10.1093/bioinformatics/btv145
