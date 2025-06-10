# Helper function: removes columns in a data.frame that contain only a single unique value
.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[["name"]] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0) {
      warning(
        paste("Dropping single category variables:"),
        paste(colnames(df)[k_singular], collapse = ", ")
      )
    }
    df <- df[, !k_singular, drop = F]
    if (ncol(df) == 0) df[["name"]] <- rownames(df)
  }
  return(df)
}

# Main function: converts a Seurat object to an AnnData object
seurat2anndata<- function(obj, outFile = NULL, assay = "RNA", main_layer = "data", transfer_layers = NULL, drop_single_values = TRUE) {
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  main_layer <- match.arg(main_layer, c("data", "counts", "scale.data"))
  transfer_layers <- transfer_layers[
    transfer_layers %in% c("data", "counts", "scale.data")
  ]
  transfer_layers <- transfer_layers[transfer_layers != main_layer]
  
  if (compareVersion(as.character(obj@version), "3.0.0") < 0) {
    obj <- Seurat::UpdateSeuratObject(object = obj)
  }
  
  X <- Seurat::GetAssayData(object = obj, layer = assay, slot = main_layer)
  
  obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)
  dx<-Seurat::GetAssay(obj, assay = assay)@meta.data #Seurat(V5)
  rownames(dx)=rownames(Seurat::GetAssay(obj, assay = assay))
  var <- .regularise_df(dx, drop_single_values = drop_single_values)
  
  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction = name)),
      simplify = FALSE
    )
    names(obsm) <- paste0("X_", tolower(names(obj@reductions)))
  }
  
  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = assay, slot = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }
  
  anndata <- reticulate::import("anndata", convert = FALSE)
  
  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )
  
  if (!is.null(outFile)) {
    adata$write(outFile, compression = "gzip")
  }
  
  adata
}
