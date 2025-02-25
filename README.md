# Original code to predict gene-peak model from paired scRNA-seq and scATAC-seq or only scATAC-seq and evaluate models

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
