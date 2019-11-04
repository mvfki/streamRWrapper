# THIS SCRIPT IS STILL A TEST NOTE

# Basically with package "reticulate", I can call most of the Python functions 
# from R. I guess now what I want to do is to just run the same thing here and 
# reproduce the results that I can get from directly running in Python. 
# Next thing to do is to make things more R, including converting Python Objects
# saved visualizable in R, callable with R functions. 

#####Load basic environment and packages########################################

#setwd('/mnt/d/BU_MS_BF/camplab/work/git/streamRWrapper/')
if("reticulate" %in% rownames(installed.packages()) == FALSE) {
    install.packages('reticulate')
    }
library(reticulate)
use_condaenv(condaenv = 'STREAM', required = TRUE)
Sys.setenv(LD_LIBRARY_PATH="~/.conda/envs/STREAM/lib/")
st <- import('stream')

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if("SingleCellExperiment" %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install('SingleCellExperiment')
} 
library('SingleCellExperiment')

#####Functions that work and might be useful####################################

convertAdata <- function(PyAnnData) {
    # This function creates an RAnnData Object from the given AnnData from 
    # Python. The structure imitates the Python style. 
    RAnnData <- new("RAnnData", X = PyAnnData$X, obs = PyAnnData$obs, 
                    var = PyAnnData$var, uns = PyAnnData$uns, 
                    obsm = PyAnnData$obsm$as_dict())
    return(RAnnData)
}
setClass("RAnnData", slots = c(X = "matrix", obs = "data.frame", 
                               var = 'data.frame', uns = 'list', obsm = 'list'))

plotUMAP2D <- function(Adata) {
    # This function plots the STREAM style UMAP visualization in R
    if ((class(Adata) == c("anndata.core.anndata.AnnData", 
                           "python.builtin.object"))[1] && 
        (class(Adata) == c("anndata.core.anndata.AnnData", 
                           "python.builtin.object"))[2]) {
        # Condition that the input adata is the Python AnnData Object. 
        if (is.null(adata$obsm$get('X_vis_umap'))) {
            write('Calculating...', stdout())
            # TODO support other methods later
            st$plot_visualization_2D(Adata, method = method)
        } else {
        }
        write('Importing calculated UMAP visualization', stdout())
        Vis <- Adata$obsm$get("X_vis_umap")
        uniqLabels <- unique(Adata$obs$label)
        allColors <- Adata$obs$label_color
    } else if (class(Adata) == 'RAnnData') {
        if (is.null(Adata@obsm$X_vis_umap)) {
            stop("The UMAP visualization arrays are not found in the input 
RAnnData. Try input the Python adata.AnnData Object and then use 
convertAdata(adata) to get the visualizable RAnnData")
        } else {
            Vis <- Adata@obsm$X_vis_umap
            uniqLabels <- unique(Adata@obs$label)
            allColors <- Adata@obs$label_color
        }
    } else {
        stop("Please give a correct AnnData Object (either an anndata.AnnData 
Object in Python space, or an RAnnData Object converted from Python). ")
    }
    par(mar = c(3, 3, 1, 1))
    plot(Vis[,1], Vis[,2], pch = 16, cex=0.8, col = allColors, xlab = '', 
         ylab = '')
    legend("topright", legend = uniqLabels, pch = 16, col = unique(allColors), 
           bty = 'n', cex=0.8)
}

adata2sce <- function(PyAnnData) {
    ## Wrap a Python anndata.AnnData object to SingleCellExperiment object
    # First check if the input data is parsable. 
    if ((class(PyAnnData) == c("anndata.core.anndata.AnnData", 
                               "python.builtin.object"))[1] && 
        (class(PyAnnData) == c("anndata.core.anndata.AnnData", 
                               "python.builtin.object"))[2]) {
        write("Parsing the input AnnData...", stdout())
    } else {
        stop("Please input a Python AnnData object. ")
    }
    # The matrix. Since SCE requires the dimension should match at many places, 
    # only the filtered information will be added rather than the raw count of 
    # all cells and genes before filtration. 
    calculatedMatrix <- data.frame(t(PyAnnData$X), 
                                   row.names = PyAnnData$var_names$to_list())
    colnames(calculatedMatrix) <- PyAnnData$obs_names$to_list()
    rawCount <- data.frame(t(PyAnnData$raw$X), 
                           row.names = PyAnnData$raw$var_names$to_list())
    colnames(rawCount) <- PyAnnData$raw$obs_names$to_list()
    filteredRowNames <- row.names(calculatedMatrix)
    filteredColNames <- colnames(calculatedMatrix)
    filteredRawCount <- rawCount[filteredRowNames, filteredColNames]
    sce <- SingleCellExperiment(assays = 
                                    list(counts = as.matrix(filteredRawCount), 
                                         stream_matrix = 
                                             as.matrix(calculatedMatrix)))
    # For gene information
    genes <- PyAnnData$var_names$to_list()
    gene_ids <- PyAnnData$var$gene_ids
    n_counts <- PyAnnData$var$n_counts
    n_cells <- PyAnnData$var$n_cells
    sce@int_elementMetadata@rownames <- genes
    sce@int_elementMetadata@listData$gene_ids <- gene_ids
    sce@int_elementMetadata@listData$n_counts <- n_counts
    sce@int_elementMetadata@listData$n_cells <- n_cells
    sce@elementMetadata@rownames <- genes
    sce@elementMetadata@listData$gene_ids <- gene_ids
    sce@elementMetadata@listData$n_counts <- n_counts
    sce@elementMetadata@listData$n_cells <- n_cells
    rowData(sce) <- list(geneSymbol = genes, gene_ids = gene_ids, 
                         n_counts = n_counts, n_cells = n_cells)
    # For cell information
    cells <- PyAnnData$obs_names$to_list()
    label <- PyAnnData$obs$label
    label_color <- PyAnnData$obs$label_color
    n_counts <- PyAnnData$obs$n_counts
    n_genes <- PyAnnData$obs$n_genes
    sce@int_colData@rownames <- cells
    sce@int_colData@listData$barcodes <- cells
    sce@int_colData@listData$label <- label
    sce@int_colData@listData$label_color <- label_color
    sce@int_colData@listData$n_counts <- n_counts
    sce@int_colData@listData$n_genes <- n_genes
    sce@colData@rownames <- cells
    sce@colData@listData$barcodes <- cells
    sce@colData@listData$label <- label
    sce@colData@listData$label_color <- label_color
    sce@colData@listData$n_counts <- n_counts
    sce@colData@listData$n_genes <- n_genes
    #colData(sce) <- list(barcodes = cells, label = label, 
    #                     label_color = label_color, n_counts = n_counts, 
    #                     n_genes = n_genes)
    # Other calculation results
    # TODO: add things in PyAnnData$obsm; Be careful if SingleCellExperiment 
    # objectrequires dimension match. 
    
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
# In small test case I say "1" there, but in real data remember to say "5". 
st$filter_genes(adata, min_num_cells = max(1L, adata$n_obs * 0.001))
st$select_variable_genes(adata, n_genes = min(2000L, adata$n_vars))
st$dimension_reduction(adata, nb_pct = 0.01)
plotUMAP2D(adata)
Radata <- convertAdata(adata)
sce <- adata2sce(adata)
#####DEBUGGING AREA#############################################################
