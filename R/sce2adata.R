#' Coverts SingleCellExperiment object from R to anndata.AnnData object in Python
#'
#' The AnnData object here can be saved to .h5ad file and read into Python
#' interactive console.
#' @param SCE A SingleCellExperiment object.
#' @param assay Character, optional, default "counts". The name of assay of
#' interest that will be set as the primary matrix of the output AnnData.
#' Available options can be listed by `names(SCE@assays$data)`
#' @param save_h5ad Boolean, optional, default FALSE. Whether to save the
#' converted AnnData object to file.
#' @param h5ad_filename character, optional, default NULL.
#' @return A Python anndata.AnnData object
#' @export
#' @examples
#' adata <- adata2sce(sce)
sce2adata <- function(SCE, assay = 'counts', save_h5ad = FALSE, h5ad_filename = NULL) {
    # Transfer SCE object back to AnnData
    # Argument check first
    stopifnot(class(SCE) == "SingleCellExperiment")
    if (save_h5ad) {
        stopifnot(!is.null(h5ad_filename))
    }
    # r-reticulate required, python anndata required
    anndata <- import('anndata')
    # Extract information that correspond to AnnData structure
    X <- t(SCE@assays$data[[assay]])
    obs <- as.data.frame(SCE@colData)
    var <- as.data.frame(SCE@elementMetadata)
    AnnData <- anndata$AnnData(X = X, obs = obs, var = var)
    # For uns
    AnnData$uns <- sce@metadata
    # For obsm
    obsm_names <- names(SCE@reducedDims)
    for (i in 1:length(obsm_names)) {
        AnnData$obsm$'__setitem__'(obsm_names[i],
                                   sce@reducedDims[[obsm_names[i]]])
    }
    # Furthermore, the other assays will for now also be saved to .obsm
    assayNames <- names(SCE@assays$data)
    for (i in 1:length(assayNames)) {
        name <- assayNames[i]
        if (!name == assay) {
            AnnData$obsm$'__setitem__'(name, t(SCE@assays$data[[name]]))
        }
    }
    return(AnnData)
}
