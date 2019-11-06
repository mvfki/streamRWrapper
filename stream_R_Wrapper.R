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
if("SingleCellExperiment" %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install('SingleCellExperiment')
} 
library(SingleCellExperiment)

## Import package
st <- import('stream')
anndata <- import('anndata')

#####Functions that work and might be useful####################################
plot_AnnData_UMAP_2D <- function(AnnData) {
	# This function extract the precomputed UMAP visualization coordinates and 
	# plot in R space. 
    AnnDataClass <- c("anndata.core.anndata.AnnData", "python.builtin.object")
    if (!(FALSE %in% (class(AnnData) == AnnDataClass))) {
        # Condition that the input annData is the Python AnnData Object. 
        if (is.null(AnnData$obsm$get('X_vis_umap'))) {
            write('Calculating...', stdout())
            # TODO support other methods later
            st$plot_visualization_2D(AnnData)
        } 
        write('Importing calculated UMAP visualization', stdout())
        Vis <- AnnData$obsm$get("X_vis_umap")
        uniqLabels <- unique(AnnData$obs$label)
        allColors <- AnnData$obs$label_color
    } else {
        stop("Please give a correct AnnData Object.")
    }
    par(mar = c(3, 3, 1, 1))
    plot(Vis[,1], Vis[,2], pch = 16, cex=0.8, col = allColors, xlab = '', 
         ylab = '')
    legend("topright", legend = uniqLabels, pch = 16, col = unique(allColors), 
           bty = 'n', cex=0.8)
}

adata2sce <- function(AnnData) {
    # Wrap a Python anndata.AnnData object to SingleCellExperiment object. 
    # I know there is a module called anndata2ri but ...
    
    # First check if the input data is parsable. 
    AnnDataClass <- c("anndata.core.anndata.AnnData", "python.builtin.object")
    if (!(FALSE %in% (class(adata) == AnnDataClass))) {
        write("Parsing the input AnnData...", stdout())
    } else {
        stop("Please input a Python AnnData object. ")
    }
    
    # The matrix. Since SCE requires the dimension should match at many places, 
    # only the filtered information will be added rather than the raw count of
    # all cells and genes before filtration
    calculatedMatrix <- data.frame(t(AnnData$X), 
                                   row.names = AnnData$var_names$to_list())
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
                                reducedDims = adata$obsm$as_dict())

    # For gene information
    genes <- AnnData$var_names$to_list()
    var <- as.list(AnnData$var)
    sce@int_elementMetadata@rownames <- genes
    sce@int_elementMetadata@listData <- var
    rowData(sce) <- list(gene_ids = gene_ids, n_counts = n_counts, 
                         n_cells = n_cells)
    
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
st$select_variable_genes(adata, n_genes = min(2000L, adata$n_vars))
st$dimension_reduction(adata, nb_pct = 0.01)
plot_AnnData_UMAP_2D(adata)
sce <- adata2sce(adata)
#####DEBUGGING AREA#############################################################
