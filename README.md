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
