#' Plot 2D visualization in STREAM style
#'
#' This function imports the STREAM style UMAP/tSNE visualization coordinates
#' and plot them in R space. Alternatively, the plot can be labeled by a naive
#' K-Means clustering in order to remove small outlier cell groups.
#' @param AnnData A Python anndata.AnnData object.
#' @param method Character, default "umap". The calculation method for the 2D
#' visualization
#' @param label_by_cluster Boolean, default FALSE. Whether to plot with
#' clustering labels.
#' @param n_cluster int, default 15L. Number of clusters when plot with labels
#' of clusters.
#' @param save_fig Boolean, default FALSE. Whether to save the plot
#' automatically.
#' @param fig_name character, default "STREAM_<method>_visualization_2D.pdf".
#' @return clusters. Will be `NULL`` unless `label_by_cluster` is set to `TRUE`.
#' An array of int, each of which specifies the cluster a cell belongs to.
#' @examples
#' plot_AnnData_UMAP_2D(adata)
plot_AnnData_UMAP_2D <- function(AnnData, method = "umap",
                                 label_by_cluster = FALSE, n_cluster = 15L,
                                 save_fig = FALSE,
                                 fig_name = "visualization_2D.pdf") {
    # This function imports the STREAM style UMAP/tSNE visualization
    # coordinates and plot them in R space.
    AnnDataClass <- c("anndata.core.anndata.AnnData", "python.builtin.object")
    if (!(FALSE %in% (class(AnnData) == AnnDataClass))) {
        # Condition that the input annData is the Python AnnData Object.
        if (method == "umap") {
            obsmKey <- "X_vis_umap"
        } else if (method == "tsne") {
            obsmKey <- "X_vis_tsne"
        } else {
            stop("Only 'umap' and 'tsne' supported for the plotting method.")
        }
        st$plot_visualization_2D(AnnData, method = method,
                                 save_fig = save_fig,
                                 fig_name = paste("STREAM", method,
                                                  fig_name, sep = "_"))
        write(paste('Importing calculated', method,
                    'visualization', sep = ' '),
              stdout())
        Vis <- AnnData$obsm$get(obsmKey)
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
    if (label_by_cluster) {
        # If we want tot label the umap by clustering for removing the outliers
        cl <- kmeans(Vis, n_cluster)
        palette <- rainbow(n_cluster)
        allColors <- palette[cl$cluster]
    }
    par(mar = c(3, 3, 1, 1))
    plot(Vis[,1], Vis[,2], pch = 16, cex=0.8, col = allColors, xlab = '',
         ylab = '')
    if (label_by_cluster) {
        text(cl$centers[,1], cl$centers[,2], labels = 1:n_cluster)
        return(cl$cluster)
    } else {
        legend("topright", legend = uniqLabels, pch = 16, col = unique(allColors),
               bty = 'n', cex=0.8)
    }
}
