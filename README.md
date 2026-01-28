# Original code for predicting peak-gene models from paired scATAC-seq/RNA-seq or scATAC-seq data alone, and model evaluation

## Introduction

We present a novel model, SCEG-HiC, which predicts links between genes and enhancers by integrating multi-omics sequencing data (paired scATAC-seq/RNA-seq or scATAC-seq alone) with three-dimensional chromatin conformation data (bulk average Hi-C). Our approach employs the weighted graphical lasso (wglasso) model to incorporate average bulk Hi-C data, effectively regularizing the correlation matrix with the prior Hi-C contact matrix as a penalty term.

To validate the accuracy of peak-gene links predicted by SCEG-HiC, we benchmarked its predictions against other tools that infer peak-gene links using single-cell omics data. For paired scATAC-seq/RNA-seq data, we compared SCEG-HiC with Pearson, Spearman, CellOracle, DIRECT-NET, SCENIC+, FigR, Signac, SCARlink, Enhlink, and scMultiMap. For scATAC-seq data alone, we compared it with Cicero, DIRECT-NET, chi-squared test with FDR correction (Chi2+FDR), scEChIA, Pearson and Enhlink. Since SCEG-HiC employs a strategy similar to the activity-by-contact (ABC) model by using bulk average Hi-C data instead of cell type-specific Hi-C, we also included ABC in the comparison. The accuracy of predicted peak-gene links was evaluated using both cell type-specific Hi-C and eQTL data.

## Download data 

We evaluated and applied the model using a total of five human paired scATAC-seq/RNA-seq datasets, five mouse paired scATAC-seq/RNA-seq datasets, and COVID-19 scATAC-seq datasets. Data processing was performed using Seurat and Signac. The processed datasets are publicly available and can be downloaded from [10.5281/zenodo.14849886](https://zenodo.org/record/14849886). Here, we use the publicly available 10x Genomics Multiome dataset, ```PBMC_multiomic.rds```, as an example to demonstrate how to run the code for predicting peak–gene links models.

## Predict peak-gene model

We provided code for predicting peak-gene models using single-cell omics data. The directory ```|-- code/model/paired_model``` includes methods for predicting gene-peak models from paired scATAC-seq/RNA-seq data, such as Pearson, Spearman, CellOracle, DIRECT-NET, SCENIC+, FigR, Signac, SCARlink, Enhlink, and scMultiMap.. The directory ```|-- code/model/only_scATAC_model``` includes methods for predicting peak-gene models using scATAC-seq data alone, including Cicero, DIRECT-NET, Chi2+FDR, scEChIA, Pearson and Enhlink. The directory ```|-- code/model/ABC``` contains a method for predicting peak-gene models using the ABC model.

`NOTE`:
- **Promoter**: 1 kp upstream and downstream of the gene transcription start site (TSS).
- **Enhancer**: 250 kb upstream and downstream of the gene TSS.

Table 1 below shows **information on predicting peak-gene models from scATAC-seq/RNA-seq**.

| Model       | Platform  | method                | URL                                                   |
| ----------- | --------- | --------------------- | ----------------------------------------------------- |
| Pearson     | R         | Pearson               |                                                       |
| Spearman    | R         | Spearman              |                                                       |
| CellOracle  | R; Python | Cicero (glasso)       | https://github.com/morris-lab/CellOracle              |
| DIRECT-NET  | R         | XGBooset              | https://github.com/zhanglhbioinfor/DIRECT-NET         |
| SCENIC+     | R; Python | Spearman; GBM         | https://github.com/aertslab/scenicplus                |
| FigR        | R         | Spearman              | https://github.com/buenrostrolab/FigR                 |
| Siganc      | R         | Pearson               | https://github.com/stuart-lab/signac                  |
| SCARlink    | R; Python | Poisson regression    | https://github.com/snehamitra/SCARlink                |
| Enhlink     | R; Python | Random forest-like    | https://gitlab.com/Grouumf/enhlinktools               |
| scMultiMap  | R         | Joint latent variable | https://github.com/ChangSuBiostats/scMultiMap         |
| SCEG-HiC    | R         | wglasso               | https://github.com/wuwei77lx/SCEGHiC                  |


Table 2 below shows **information on predicting peak-gene models from scATAC-seq data alone**.


| Model       | Platform  | method             | URL                                                   |
| ----------- | --------- | ------------------ | ----------------------------------------------------- |
| Cicero      | R         | glasso             | https://github.com/cole-trapnell-lab/cicero-release   |
| DIRECT-NET  | R         | XGBoost            | https://github.com/zhanglhbioinfor/DIRECT-NET         |
| Chi2+FDR    | R         | chi-square         |                                                       |
| scEChIA     | R         | glasso             | https://github.com/reggenlab/scEChiA                  |
| Pearson     | R         | Pearson            |                                                       |
| Enhlink     | R; Python | Random forest-like | https://gitlab.com/Grouumf/enhlinktools               |
| SCEG-HiC    | R         | wglasso            | https://github.com/wuwei77lx/SCEGHiC                  |


Table 3 below shows **information on predicting peak-gene model using the ABC model**.


| Model       | Platform  | method          | URL                                                                                      |
| ----------- | --------- | --------------- | ---------------------------------------------------------------------------------------- |
| ABC         | R; Python | score-based     | https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019 |


## Model evaluation
We evaluated the models using cell type-specific Hi-C data and eQTL data. The directory ```|-- code/evaluation``` contains evaluation metrics (AUPRC and early precision) as well as scripts for processing the cell type-specific Hi-C and eQTL data.

**Note:** The cell type-specific Hi-C data used for validation are independent and were not included in the bulk average Hi-C data.

### Data (Hi-C and eQTL) Processing:
- **Hi-C**: Normalized scores (SCALCE or KR) across 5 kb bins are extracted using Juicer Tools. Only interactions with appropriate scores are kept, where one bin overlaps at least one consensus peak and the other overlaps a TSS.
- **eQTL**: Only links where the peak overlaps an eQTL locus, and the locus is associated with the expression of the corresponding gene, are retained.

### Evaluation metrics:
- **AUPRC**: Area under the precision-recall curve.
- **Early precision**: Fraction of true positives in the top-k edges (where k equals the number of edges in the ground truth).
- **AUPRC ratio and early precision ratio**: Calculated by dividing the AUPRC and early precision values by the precision of a random predictor, which is determined based on the ground-truth network.

## The bulk average Hi-C data
In the [ABC](https://www.nature.com/articles/s41588-019-0538-0) model, it has been demonstrated that using an average Hi-C profile yields approximately the same performance as using cell type-specific Hi-C data. To enable model predictions across a broad range of cell types, including those lacking cell type-specific Hi-C data, we provide an average Hi-C contact matrix for use in model inference.

### Human average Hi-C (provided by ABC model)

The human cell types used for generating the average Hi-C profile include: [34 Hi-C datasets](https://www.encodeproject.org/annotations/ENCSR382HAW/) from the ENCODE project.

**The human average Hi-C data (~54.1 GB) can be downloaded from**:
https://www.encodeproject.org/files/ENCFF134PUN/@@download/ENCFF134PUN.bed.gz

To extract the human bulk average Hi-C data, use the script located at ```|-- code/model/ABC/extract_avg_hic.py```:

```{sh eval=FALSE}
python code/model/ABC/extract_avg_hic.py --avg_hic_bed_file ../ENCFF134PUN.bed.gz --output_dir ../
```

### Mouse bulk average Hi-C
The mouse cell types used for generating the average Hi-C profile include: two embryonic stem cell types (mESC1, mESC2) CH12LX, CH12F3, fiber, epithelium, B cells.

Table 3 below shows **the sources of mouse Hi-C datasets used to generate the average Hi-C data**.


| Cell type                                  | Source      | URL                                                                                                                             |
| ------------------------------------------ | ---------   | ------------------------------------------------------------------------------------------------------------------------------- |
| Mouse embryonic stem cell (C57BL6/129s4)   | GSE124193   | https://www.encodeproject.org/files/ENCFF409ZRS/@@download/ENCFF409ZRS.hic                                                      |
| Mouse embryonic stem cell (strain 129/Ola) | GSE82144    | https://www.encodeproject.org/files/ENCFF584EDJ/@@download/ENCFF584EDJ.hic                                                      |
| CH12F                                      | GSE98119    | https://www.encodeproject.org/files/ENCFF909ODS/@@download/ENCFF909ODS.hic                                                      |
| Mature B cells                             | GSE82144    | https://www.encodeproject.org/files/ENCFF026SNO/@@download/ENCFF026SNO.hic                                                      |
| CH12.LX                                    | GSE63525    | https://www.encodeproject.org/files/ENCFF076GYR/@@download/ENCFF076GYR.hic                                                      |
| Epithelium                                 | GSE243851   | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243851&format=file&file=GSE243851%5Fepi%5Fboth%5Fsamples%5Finter%5F30%2Ehic   |
| Fiber                                      | GSE243851   | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243851&format=file&file=GSE243851%5Ffiber%5Fboth%5Fsamples%5Finter%5F30%2Ehic |


For detailed steps on generating the average Hi-C, please refer to the scripts and documentation at:
https://github.com/wuwei77lx/compare_model/blob/main/code/model/ABC/README.md

**The mouse average Hi-C data (~13.9 GB) can be downloaded from**: [10.5281/zenodo.14849886](https://zenodo.org/record/14849886).

## Reference

CellOracle: Kamimoto, K., Stringa, B., Hoffmann, C. M., Jindal, K., Solnica-Krezel, L., & Morris, S. A. (2023). Dissecting cell identity via network inference and in silico gene perturbation. Nature, 614(7949), 742–751. https://doi.org/10.1038/s41586-022-05688-9

DIRECT-NET: Zhang, L., Zhang, J., & Nie, Q. (2022). DIRECT-NET: An efficient method to discover cis-regulatory elements and construct regulatory networks from single-cell multiomics data. Science advances, 8(22), eabl7393. https://doi.org/10.1126/sciadv.abl7393

SCENIC+: Bravo González-Blas, C., De Winter, S., Hulselmans, G., Hecker, N., Matetovici, I., Christiaens, V., Poovathingal, S., Wouters, J., Aibar, S., & Aerts, S. (2023). SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. Nature methods, 20(9), 1355–1367. https://doi.org/10.1038/s41592-023-01938-4

FigR: Kartha, V. K., Duarte, F. M., Hu, Y., Ma, S., Chew, J. G., Lareau, C. A., Earl, A., Burkett, Z. D., Kohlway, A. S., Lebofsky, R., & Buenrostro, J. D. (2022). Functional inference of gene regulation using single-cell multi-omics. Cell genomics, 2(9), 100166. https://doi.org/10.1016/j.xgen.2022.100166

Signac:Stuart, T., Srivastava, A., Madad, S., Lareau, C. A., & Satija, R. (2021). Single-cell chromatin state analysis with Signac. Nature methods, 18(11), 1333–1341. https://doi.org/10.1038/s41592-021-01282-5

Cicero: Pliner, H. A., Packer, J. S., McFaline-Figueroa, J. L., Cusanovich, D. A., Daza, R. M., Aghamirzaie, D., Srivatsan, S., Qiu, X., Jackson, D., Minkina, A., Adey, A. C., Steemers, F. J., Shendure, J., & Trapnell, C. (2018). Cicero Predicts cis-Regulatory DNA Interactions from Single-Cell Chromatin Accessibility Data. Molecular cell, 71(5), 858–871.e8. https://doi.org/10.1016/j.molcel.2018.06.044

scEChIA: Pandey, N., Omkar Chandra, Mishra, S., & Kumar, V. (2021). Improving Chromatin-Interaction Prediction Using Single-Cell Open-Chromatin Profiles and Making Insight Into the Cis-Regulatory Landscape of the Human Brain. Frontiers in genetics, 12, 738194. https://doi.org/10.3389/fgene.2021.738194

ABC: Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0


## Contact
If you have any problems, comments or suggestions, please contact us at XuanLiang (<liangxuan2022@sinh.ac.cn>).
