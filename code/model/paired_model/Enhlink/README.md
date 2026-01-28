# Enhlink

## Running Enhlink

### Step 1. Processing a multi-omics dataset from Seurat and Signac R objects

scRNA-seq and scATAC-seq count data were exported from a Seurat object in MTX matrix format using the writeMM function from the Matrix R package. Cell type annotations, gene and peak information, and batch effect covariates (if available) were also extracted from the Seurat object. See the code in ```|-- code/model/Enhlink/Data_preparation```.

### Step 2. Processing the data with Enhlink

To process these data, we will use a custom GTF file in which promoters are defined as regions spanning 1 kb upstream and downstream of the transcription start site (TSS). See the data in ```|-- data/Homo_sapiens.GRCh38.99.TSS.1K.bed```.

```enhlink``` employs an ensemble strategy that integrates cell covariates and produces robust p-values for any link and covariate-specific link. Enhlink is run with ```-cellRanger``` to specify Cell Ranger–formatted matrices, ```-isExpr``` to treat the expression matrix as continuous, and ```-clusters``` to incorporate cell cluster annotations.

**Sample Command**:

```
enhlink -mat pbmc.multi.ATAC.mtx \
    -xgi pbmc.multi.ATAC.colnames \
    -ygi pbmc.multi.ATAC.rownames \
    -mat2 pbmc.multi.RNA.mtx \
    -xgi2 pbmc.multi.RNA.colnames \
    -ygi2 pbmc.multi.RNA.rownames \
    -gtf Homo_sapiens.GRCh38.99.TSS.1K.bed \
    -clusters clusters.tsv \
    -out output_multi/ \
    -threads 6 \
    -isExpr \
    -targets target_genes.txt \
    -format cellRanger \
    -threshold 1
```

Optionally, additional covariates can be included using the ```-covariates``` option, which incorporates binarized covariates into the model. The ```-secondOrder``` flag enables the analysis of second-order interactions to identify covariate-specific links, at the cost of increased computational time.

**Sample Command**:

```
enhlink -mat pbmc.multi.ATAC.mtx \
    -xgi pbmc.multi.ATAC.colnames \
    -ygi pbmc.multi.ATAC.rownames \
    -mat2 pbmc.multi.RNA.mtx \
    -xgi2 pbmc.multi.RNA.colnames \
    -ygi2 pbmc.multi.RNA.rownames \
    -gtf Homo_sapiens.GRCh38.99.TSS.1K.bed \
    -clusters clusters.tsv \
    -covariates cov.tsv \
    -secondOrder \
    -out output_multi/ \
    -threads 6 \
    -isExpr \
    -targets target_genes.txt \
    -format cellRanger \
    -threshold 1
```

