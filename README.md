# Original code to predict gene-peak model from paired scRNA-seq and scATAC-seq or only scATAC-seq and evaluate models

## Introduction

We presented a novel model, predicting single-cell enhancer-gene interactions by integrating priori Hi-C information (SCEG-HiC), which predicts interactions between genes and enhancers using multiomics (paired scRNA-seq and scATAC-seq or only scATAC-seq) sequencing data and three-dimensional omics data (average Hi-C) . The approach uses weighted graphical lasso (wglasso) model, taking the average Hi-C data as priori knowledge, to predict cell specific interactions between genes and enhancers based on the correlation between genes and enhancers.

To validate the quality of peak-gene relationships predicted by SCEG-HiC, we benchmarked its predictions to against other tools that predict peak-gene using multiomics data.For paired scRNA-seq and scATAC-seq data, compared SCEG-HiC to pearson, spearman, cellOracle, DIRECT-NET, SCENIC+, FigR, and Signac. For only scATAC-seq data, we compared it to Cicero, DIRECT-NET, chi2+fdr, scEChIA, and pearson. We assessed the quality of predicted peak-gene associations making use of  Hi-C data and eQTL data.

## Download data 

We evaluated and applied the model using a total of 5 human paired scRNA-seq and scATAC-seq data, 5 mouse paired scRNA-seq and scATAC-seq data, and COVID-19 scATAC-seq data. We used Seurat and Signac to process datas.
Data can be downloaded from: [10.5281/zenodo.14849886](https://zenodo.org/record/14849886).We used publicly available 10x Genomic Multiome dataset for human PBMCs  ```pbmc_multiomic.rds```.

## Predict gene-peak model

We provided code for predicting peak-gene models using multiomics data.```|-- code/model/paired_model``` included predict gene-peak model from paired scRNA-seq and scATAC-seq, namely pearson, spearman, cellOracle, DIRECT-NET, SCENIC+, FigR, and Signac.```|-- code/model/only_scATAC_model``` included predict gene-peak model from only scATAC-seq, namely Cicero, DIRECT-NET, chi2+fdr, scEChIA, and pearson.

`NOTE`:
- promoter:1000bp upstream and downstream of gene transcription start site(TSS).
- enhancer:250kb upstream and downstream of gene transcription start site(TSS).

Table 1 below shows **the information predict gene-peak model from paired scRNA-seq and scATAC-seq**.

| Model       | Platform  | method          | URL                                                   |
| ----------- | --------- | --------------- | ----------------------------------------------------- |
| pearson     | R         | pearson         |                                                       |
| spearman    | R         | spearman        |                                                       |
| cellOracle  | R;Python  | cicero(glasso)  | https://github.com/morris-lab/CellOracle              |
| DIRECT-NET  | R         | XGBooset        | https://github.com/zhanglhbioinfor/DIRECT-NET         |
| SCENIC+     | R;Python  | spearman;GBM    | https://github.com/aertslab/scenicplus                |
| FigR        | R         | spearman        | https://github.com/buenrostrolab/FigR                 |
| Siganc      | R         | pearson         | https://github.com/stuart-lab/signac                  |
| SCEG-HiC    | R         | wglasso         | https://github.com/wuwei77lx/SCEGHiC                  |


Table 2 below shows **the information predict gene-peak model from only scATAC-seq**.


| Model       | Platform  | method          | URL                                                   |
| ----------- | --------- | --------------- | ----------------------------------------------------- |
| cicero      | R         | glasso          | https://github.com/cole-trapnell-lab/cicero-release   |
| DIRECT-NET  | R         | XGBoost         | https://github.com/zhanglhbioinfor/DIRECT-NET         |
| chi2+fdr    | R         | chi-square      |                                                       |
| scEChIA     | R         | glasso          | https://github.com/reggenlab/scEChiA                  |
| pearson     | R         | pearson         |                                                       |
| SCEG-HiC    | R         | wglasso         | https://github.com/wuwei77lx/SCEGHiC                  |

## Model evaluation
We evaluated model  making use of  Hi-C data and eQTL data.```|-- code/evaluation``` included evaluation metrics(AUPRC and Early Precision) and data(Hi-C and eQTL) processing.

### data(Hi-C and eQTL) processing:
- Hi-C:normalized scores(SCALCE or KR) across bins of 5 kb using Juicer Tools were extracted, keeping only links with appropriate scores and involving a bin that overlaps at least one of the consensus peaks and a TSS.
- eQTL:keeping only links with the peak overlapp an eQTL locus and the locus is associated with the expression of the gene.

### evaluation metrics:
- AUPRC:area under the precision-recall curve.
- Early Precision:the fraction of true positives in the top-k edges(k equaled the number of edges in the ground-truth).

## Average Hi-C
In [Activity by Contact (ABC)](https://www.nature.com/articles/s41588-019-0538-0) model, proved using an average Hi-C profile gives approximately equally good performance as using a cell-type specific Hi-C profile.To facilitate making model predictions in a large panel of cell types, including those without cell type-specific Hi-C data, we have provided an average Hi-C matrix .

### human average Hi-C(provided by ABC model)
The celltypes used by human for averaging are: GM12878, NHEK, HMEC, RPE1, THP1, IMR90, HUVEC, HCT116, K562, KBM7.

Average Hi-C data can be downloaded from: 

**The human avg HiC file here**: https://www.encodeproject.org/files/ENCFF134PUN/@@download/ENCFF134PUN.bed.gz (54.1 GB)

Extract the human bulk average Hi-C using  ```|-- code/Hi_C/extract_avg_hic.py```.

```{sh eval=FALSE}
python code/Hi_C/extract_avg_hic.py --avg_hic_bed_file ../ENCFF134PUN.bed.gz --output_dir ../
```

### mouse average Hi-C
The celltypes used by mouse for averaging are: mESC1, mESC2, CH12LX, CH12F3, fiber, epithelium, B.

Table 3 below shows **the source of mouse Hi-C for producing average Hi-C**.


| celltype                                   | source      | URL                                                                                                                             |
| ------------------------------------------ | ---------   | ------------------------------------------------------------------------------------------------------------------------------- |
| mouse embryonic stem cell(C57BL6/129s4)    | GSE124193   | https://www.encodeproject.org/files/ENCFF409ZRS/@@download/ENCFF409ZRS.hic                                                      |
| mouse embryonic stem cell(strain 129/Ola)  | GSE82144    | https://www.encodeproject.org/files/ENCFF584EDJ/@@download/ENCFF584EDJ.hic                                                      |
| CH12F                                      | GSE98119    | https://www.encodeproject.org/files/ENCFF909ODS/@@download/ENCFF909ODS.hic                                                      |
| mature B cell                              | GSE82144    | https://www.encodeproject.org/files/ENCFF026SNO/@@download/ENCFF026SNO.hic                                                      |
| CH12.LX                                    | GSE63525    | https://www.encodeproject.org/files/ENCFF076GYR/@@download/ENCFF076GYR.hic                                                      |
| epithelium                                 | GSE243851   | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243851&format=file&file=GSE243851%5Fepi%5Fboth%5Fsamples%5Finter%5F30%2Ehic   |
| fiber                                      | GSE243851   | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243851&format=file&file=GSE243851%5Ffiber%5Fboth%5Fsamples%5Finter%5F30%2Ehic |



**producing average Hi-C consists of the following steps(```|-- code/Hi_C```)**:

1.Download hic matrix file from juicebox.

2.Fit HiC data to powerlaw model and extract parameters.

3.Make average Hi-C.

- Step 1. Download hic matrix file from juicebox(all celltypes)
  
```
wget https://github.com/aidenlab/Juicebox/releases/download/v2.13.06/juicer_tools_2.13.06.jar
python code/Hi_C/juicebox_dump.py \
--hic_file https://www.encodeproject.org/files/ENCFF409ZRS/@@download/ENCFF409ZRS.hic \
--juicebox "java -jar juicer_tools_2.13.06.jar" \
--resolution 5000 \
--outdir mouse_contact/ESC1
```

- Step 2. Fit HiC data to powerlaw model and extract parameters

```
python code/Hi_C/compute_powerlaw_fit_from_hic.py \
--hic_dir mouse_contact/ESC1 \
--outDir mouse_contact/ESC1/powerlaw \
--hic_resolution 5000
```


- Step 3. Make average Hi-C

```
python code/Hi_C/makeAverageHiC.py \
--celltypes ESC1,ESC2,CH12F3,B,CH12LX,fiber,epi \
--chromosome chr1 \
--basedir mouse_contact \
--outDir mouse_contact/average_HiC
```

**The mouse avg HiC file here**: [10.5281/zenodo.14849886](https://zenodo.org/record/14849886) (13.9 GB)

## Reference

cellOracle:Kamimoto, K., Stringa, B., Hoffmann, C. M., Jindal, K., Solnica-Krezel, L., & Morris, S. A. (2023). Dissecting cell identity via network inference and in silico gene perturbation. Nature, 614(7949), 742–751. https://doi.org/10.1038/s41586-022-05688-9

DIRECT-NET:Zhang, L., Zhang, J., & Nie, Q. (2022). DIRECT-NET: An efficient method to discover cis-regulatory elements and construct regulatory networks from single-cell multiomics data. Science advances, 8(22), eabl7393. https://doi.org/10.1126/sciadv.abl7393

SCENIC+:Bravo González-Blas, C., De Winter, S., Hulselmans, G., Hecker, N., Matetovici, I., Christiaens, V., Poovathingal, S., Wouters, J., Aibar, S., & Aerts, S. (2023). SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. Nature methods, 20(9), 1355–1367. https://doi.org/10.1038/s41592-023-01938-4

FigR:Kartha, V. K., Duarte, F. M., Hu, Y., Ma, S., Chew, J. G., Lareau, C. A., Earl, A., Burkett, Z. D., Kohlway, A. S., Lebofsky, R., & Buenrostro, J. D. (2022). Functional inference of gene regulation using single-cell multi-omics. Cell genomics, 2(9), 100166. https://doi.org/10.1016/j.xgen.2022.100166

Signac:Stuart, T., Srivastava, A., Madad, S., Lareau, C. A., & Satija, R. (2021). Single-cell chromatin state analysis with Signac. Nature methods, 18(11), 1333–1341. https://doi.org/10.1038/s41592-021-01282-5

cicero:Pliner, H. A., Packer, J. S., McFaline-Figueroa, J. L., Cusanovich, D. A., Daza, R. M., Aghamirzaie, D., Srivatsan, S., Qiu, X., Jackson, D., Minkina, A., Adey, A. C., Steemers, F. J., Shendure, J., & Trapnell, C. (2018). Cicero Predicts cis-Regulatory DNA Interactions from Single-Cell Chromatin Accessibility Data. Molecular cell, 71(5), 858–871.e8. https://doi.org/10.1016/j.molcel.2018.06.044

scEChIA:Pandey, N., Omkar Chandra, Mishra, S., & Kumar, V. (2021). Improving Chromatin-Interaction Prediction Using Single-Cell Open-Chromatin Profiles and Making Insight Into the Cis-Regulatory Landscape of the Human Brain. Frontiers in genetics, 12, 738194. https://doi.org/10.3389/fgene.2021.738194


## Contact
If you have any problems, comments or suggestions, please contact us at XuanLiang (<liangxuan2022@sinh.ac.cn>).
