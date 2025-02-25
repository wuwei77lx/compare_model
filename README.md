# Original code to predict gene-peak model from paired scRNA-seq and scATAC-seq or only scATAC-seq and evaluate models

## Introduction

We present a novel model, predicting single-cell enhancer-gene interactions by integrating priori Hi-C information (SCEG-HiC), which predicts interactions between genes and enhancers using multiomics (paired scRNA-seq and scATAC-seq or only scATAC-seq) sequencing data and three-dimensional omics data (average Hi-C) . The approach uses weighted graphical lasso (wglasso) model, taking the average Hi-C data as priori knowledge, to predict cell specific interactions between genes and enhancers based on the correlation between genes and enhancers.
To validate the quality of peak-gene relationships predicted by SCEG-HiC, we benchmarked its predictions to against other tools that predict peak-gene using multiomics data.For paired scRNA-seq and scATAC-seq data, compared SCEG-HiC to Pearson, Spearman, cellOracle, DIRECT-NET, SCENIC+, FigR, and Signac. For only scATAC-seq data, we compared it to Cicero, DIRECT-NET, chi2+fdr, scEChIA, and Pearson. We assessed the quality of predicted peak-gene associations making use of  Hi-C data and eQTL data.

## Predict gene-peak model
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
