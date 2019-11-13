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
if(!require("devtools")) {
    install.packages("devtools")
}
devtools::install_github("rstudio/reticulate", quiet = TRUE)
devtools::install_github("Albluca/distutils", quiet = TRUE)
devtools::install_github("Albluca/ElPiGraph.R", quiet = TRUE)
devtools::install_github("mvfki/streamRWapper", quiet = TRUE)
library(reticulate)
library(png)
library(scales)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if(!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    BiocManager::install('SingleCellExperiment')
}
library(SingleCellExperiment)

## Import package
st <- import('stream')

#####Simple functions###########################################################

viewPNG <- function(path) {
    img <- readPNG(path)
    plot.new()
    grid::grid.raster(img)
}

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
st$select_variable_genes(adata, save_fig = TRUE, loess_frac = 0.01,
                         fig_name = 'std_vs_means.png')
viewPNG('stream_result/std_vs_means.png')

# Currently, for all the plot that STREAM makes, I use the strategy of save it
# to file, and import the saved figure to RStudio plot panel.
# Planning to find a better way by either making it possible to directly allow
# the display of the plot like Jupyter Notebook does, or write R based plotting
# functions that draw equivalent figure.

st$dimension_reduction(adata, nb_pct = 0.01)
plot_AnnData_UMAP_2D(adata)
clusters <- plot_AnnData_UMAP_2D(adata, label_by_cluster = TRUE,
                                 n_cluster = 30)
remove_clusters_from_AnnData(adata, clusters, c(5, 13), plotCheck = TRUE)

st$seed_elastic_principal_graph(adata)

st$plot_branches(adata, save_fig = TRUE, fig_name = 'branches.png')
viewPNG('stream_result/branches.png')

st$plot_branches_with_cells(adata, save_fig = TRUE,
                            fig_name = 'branches_with_cells.png')
viewPNG('stream_result/branches_with_cells.png')

st$elastic_principal_graph(adata)

st$plot_branches(adata, save_fig = TRUE, fig_name = 'branches.png')
viewPNG('stream_result/branches.png')

st$plot_branches_with_cells(adata, save_fig = TRUE,
                            fig_name = 'branches_with_cells.png')
viewPNG('stream_result/branches_with_cells.png')

st$optimize_branching(adata)

st$extend_elastic_principal_graph(adata)

st$plot_branches(adata, save_fig = TRUE, fig_name = 'branches.png')
viewPNG('stream_result/branches.png')

st$plot_branches_with_cells(adata, save_fig = TRUE,
                            fig_name = 'branches_with_cells.png')
viewPNG('stream_result/branches_with_cells.png')

st$plot_flat_tree(adata, save_fig = TRUE, fig_name = 'flat_tree.png')
viewPNG('stream_result/flat_tree.png')

root <- 'S0'

st$subwaymap_plot(adata, root = root, save_fig = TRUE, percentile_dist = 100,
                  fig_name = 'subway_map.png')
viewPNG(paste('stream_result', root, 'subway_map.png', sep = '/'))

st$stream_plot(adata, root = root, save_fig = TRUE,
               fig_name = 'stream_plot.png')
viewPNG(paste('stream_result', root, 'stream_plot.png', sep = '/'))

# This step was said to be time consuming. Make sure previous stuffs are
# reasonable before you run the following.
st$detect_leaf_genes(adata,root = root)
leaf_genes <- adata$uns$leaf_genes

st$detect_transistion_genes(adata, root = root)
st$plot_transition_genes(adata, save_fig = TRUE)
#####Wrapping test##############################################################
sce <- adata2sce(adata)

#####DEBUGGING AREA#############################################################
anndata <- import('anndata')
sce2adata <- function(SCE, assay = 'counts', save_h5ad = FALSE, h5ad_filename = NULL) {
    # Transfer SCE object back to AnnData
    # Argument check first
    stopifnot(class(SCE) == "SingleCellExperiment")
    if (save_h5ad) {
        stopifnot(!is.null(h5ad_filename))
    }
    # Extract information that correspond to AnnData structure
    X <- t(SCE@assays$data[[assay]])
    obs <- as.data.frame(SCE@colData)
    var <- as.data.frame(SCE@elementMetadata)
    AnnData <- anndata$AnnData(X = X, obs = obs, var = var)
    # For uns

    # For obsm
    obsm_names <- names(SCE@reducedDims)
    for (i in 1:length(obsm_names)) {
        AnnData$obsm$'__setitem__'(obsm_names[i], sce@reducedDims[[obsm_names[i]]])
    }
    return(AnnData)
}
