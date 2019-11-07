#' Plot a 2D UMAP in the same way that STREAM plots.
#' @param AnnData A Python anndata.AnnData object that has 
#' anndata.AnnData.obsm['X_vis_umap'] value
#' @return NULL
#' @export
#' @examples
#' plot_AnnData_UMAP_2D(adata)
plot_AnnData_UMAP_2D <- function(AnnData) {
	# This function extract the precomputed UMAP visualization coordinates and 
	# plot in R space. 
    AnnDataClass <- c("anndata.core.anndata.AnnData", "python.builtin.object")
    if (!(FALSE %in% (class(AnnData) == AnnDataClass))) {
        # Condition that the input annData is the Python AnnData Object. 
        if (is.null(AnnData$obsm$get('X_vis_umap'))) {
            # Since environmental issue not solved yet, I can't support invoking
            # STREAM to calculate it automatically.
            stop('UMAP not calculated yet.')
            #write('Calculating...', stdout())
            # TODO support other methods later
            #st$plot_visualization_2D(AnnData)
        } 
        write('Importing calculated UMAP visualization', stdout())
        Vis <- AnnData$obsm$get("X_vis_umap")
        uniqLabels <- unique(AnnData$obs$label)
        allColors <- AnnData$obs$label_color
    } else {
        stop("Please give a correct AnnData Object.")
    }
    if (is.factor(allColors)) {
        # With dev-ver reticulate, the allColors extracted is a Factor object. 
        c <- array()
        colorNames <- levels(allColors)
        for (i in 1:length(colorNames)) {
            c[which(allColors == colorNames[i])] <- colorNames[i]
        }
        allColors <- c
    }
    par(mar = c(3, 3, 1, 1))
    plot(Vis[,1], Vis[,2], pch = 16, cex=0.8, col = allColors, xlab = '', 
         ylab = '')
    legend("topright", legend = uniqLabels, pch = 16, col = unique(allColors), 
           bty = 'n', cex=0.8)
}

