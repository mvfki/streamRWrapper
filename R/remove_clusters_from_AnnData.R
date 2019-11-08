#' Inplace Remove Small Outlier Groups from AnnData
#'
#' @param AnnData A Python anndata.AnnData object.
#' @param clusters Array of int. The array that `plot_AnnData_UMAP_2D(adata, label_by_cluster = TRUE)` returns
#' @param ToRemove numeric. Can be an array of numbers or just one number that 
#' specify the small groups in the plot. `c(1, 2)`, `15L`, `c('1', 2)` are all
#' acceptable
#' @param plotCheck Boolean, default FALSE. Whether to directly plot the 
#' filtered AnnData to have a check. The coloring follows the originial labels 
#' from the AnnData. 
#' @param plotCheckMethod Character, default "umap", can only be chosen from 
#' "umap" and "tsne". This would serve as the `method` argument when this 
#' function internally calls `plot_AnnData_UMAP_2D(adata, method = plotCheckMethod)`.
#' @examples
#' clusters <- plot_AnnData_UMAP_2D(adata, label_by_cluster = TRUE, 
#'                                  n_cluster = 30)
#' remove_clusters_from_AnnData(adata, clusters, c(18, 29), plotCheck = TRUE)
remove_clusters_from_AnnData <- function(AnnData, clusters, ToRemove,
                                         plotCheck = FALSE,
                                         plotCheckMethod = 'umap') {
    AnnDataClass <- c("anndata.core.anndata.AnnData", "python.builtin.object")
    if (!(FALSE %in% (class(AnnData) == AnnDataClass))) {
        stopifnot(AnnData$n_obs == length(clusters))
        cellToRemove <- c()
        for (i in 1:length(ToRemove)) {
            cellToRemove <- c(cellToRemove, which(clusters == ToRemove[i]))
        }
        # Now in cellToRemove, they are the 1-based index of the cell to remove,
        # be careful when operating the AnnData since python index is 0-based.
        cellToRemove <- unique(sort(cellToRemove))
        remain <- setdiff(1:AnnData$n_obs, cellToRemove)
        AnnData$'_inplace_subset_obs'(remain - 1L)
        if (plotCheck) {
            plot_AnnData_UMAP_2D(AnnData, method = plotCheckMethod)
        }
    } else {
        stop("Please give a correct AnnData Object.")
    }
}