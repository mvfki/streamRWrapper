# THIS SCRIPT IS STILL A TEST NOTE

# Basically with package "reticulate", I can call most of the Python functions
# from R. I guess now what I want to do is to just run the same thing here and
# reproduce the results that I can get from directly running in Python.
# Next thing to do is to make things more R, including converting Python Objects
# saved visualizable in R, callable with R functions.

#####Load basic environment and packages########################################

## Working directry
#setwd('~/camplab/work/git/streamRWrapper/')
#setwd('/mnt/d/BU_MS_BF/camplab/work/git/streamRWrapper/')

## Loading library
library(devtools)
devtools::install_github("rstudio/reticulate", quiet = TRUE)
library(reticulate)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if(!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    BiocManager::install('SingleCellExperiment')
}
library(SingleCellExperiment)

## Import package
st <- import('stream')

#####Functions that work and might be useful####################################
plot_AnnData_UMAP_2D <- function(AnnData, method = "umap",
                                 label_by_cluster = FALSE, n_cluster = 15,
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
#TODO: Check for plotting method that works from SingleCellExperiment object.

#####Pipeline###################################################################
adata <- st$read('testData_real/matrix.mtx', file_format = 'mtx')
st$add_cell_labels(adata, file_name = 'testData_real/cell_label.tsv')
st$add_cell_colors(adata, file_name = 'testData_real/cell_label_color.tsv')
adata$obs_names_make_unique()
adata$var_names_make_unique()
st$remove_mt_genes(adata)
st$normalize_per_cell(adata)
st$log_transform(adata)
st$filter_cells(adata)
# In small test case use 1L, in real test case use 5L.
st$filter_genes(adata, min_num_cells = max(5L, adata$n_obs * 0.001))
st$select_variable_genes(adata, n_genes = min(2000L, adata$n_vars),
                         save_fig = TRUE)
# By specifying save_fig to TRUE, a file called 'stream_result/std_vs_mean.pdf'
# is saved at your working directory.
st$dimension_reduction(adata, nb_pct = 0.01)
plot_AnnData_UMAP_2D(adata)
clusters <- plot_AnnData_UMAP_2D(adata_backup, label_by_cluster = TRUE,
                                 n_cluster = 30)
remove_clusters_from_AnnData(adata, clusters, c(10, 8), plotCheck = TRUE)
sce <- adata2sce(adata)
#####DEBUGGING AREA#############################################################
