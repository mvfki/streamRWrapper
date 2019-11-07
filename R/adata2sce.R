#' Coverts AnnData object from Python to SingleCellExperiment object in R. 
#' 
#' The AnnData object here can be read from .h5ad file with the help of 
#' `reticulate`.
#' @param AnnData A Python adata.AnnData object. 
#' @return A SingleCellExperiment object
#' @export
#' @examples
#' sce <- adata2sce(adata)
adata2sce <- function(AnnData) {
    # Wrap a Python anndata.AnnData object to SingleCellExperiment object. 
    # I know there is a module called anndata2ri but ...
    
    # First check if the input data is parsable. 
    AnnDataClass <- c("anndata.core.anndata.AnnData", "python.builtin.object")
    if (!(FALSE %in% (class(AnnData) == AnnDataClass))) {
        write("Parsing the input AnnData...", stdout())
    } else {
        stop("Please input a Python AnnData object. ")
    }
    
    # The matrix. Since SCE requires the dimension should match at many places, 
    # only the filtered information will be added rather than the raw count of
    # all cells and genes before filtration
    calculatedMatrix <- data.frame(t(AnnData$X), 
                                   row.names = AnnData$var_names$to_list())
    if (!is.null(AnnData$raw)) {
        colnames(calculatedMatrix) <- AnnData$obs_names$to_list()
        rawCount <- data.frame(t(AnnData$raw$X), 
                               row.names = AnnData$raw$var_names$to_list())
        colnames(rawCount) <- AnnData$raw$obs_names$to_list()
        filteredRowNames <- row.names(calculatedMatrix)
        filteredColNames <- colnames(calculatedMatrix)
        filteredRawCount <- rawCount[filteredRowNames, filteredColNames]
        sce <- SingleCellExperiment(assays = 
                                      list(counts = as.matrix(filteredRawCount), 
                                           stream_matrix = 
                                             as.matrix(calculatedMatrix)), 
                                    reducedDims = AnnData$obsm$as_dict())
    } else {
        sce <- SingleCellExperiment(assays = 
                                        list(stream_matrix = 
                                                 as.matrix(calculatedMatrix)), 
                                    reducedDims = AnnData$obsm$as_dict())
    }
    

    # For gene information
    genes <- AnnData$var_names$to_list()
    var <- as.list(AnnData$var)
    sce@int_elementMetadata@rownames <- genes
    sce@int_elementMetadata@listData <- var
    # The following two line might not be allowed in R 3.6.1
    sce@elementMetadata@rownames <- genes
    sce@elementMetadata@listData <- var
    
    rowData(sce) <- var
    
    # For cell information
    cells <- AnnData$obs_names$to_list()
    obs <- as.list(AnnData$obs)
    sce@int_colData@rownames <- cells
    sce@colData@rownames <- cells
    for (i in 1:length(obs)) {
        sce@int_colData@listData[[names(obs)[i]]] <- obs[[names(obs)[i]]]
        sce@colData@listData[[names(obs)[i]]] <- obs[[names(obs)[i]]]
    }

    return(sce)
}
