# THIS SCRIPT IS STILL A TEST NOTE

# Basically with package "reticulate", I can call most of the Python functions 
# from R. I guess now what I want to do is to just run the same thing here and 
# reproduce the results that I can get from directly running in Python. 
# Next thing to do is to make things more R, including converting Python Objects
# saved visualizable in R, callable with R functions. 

# This first section would save the commands that worked everywhere so far.
library(reticulate)
use_condaenv(condaenv = 'STREAM', required = TRUE)
Sys.setenv(LD_LIBRARY_PATH="~/.conda/envs/STREAM/lib/")

#####These commands worked on my local but not yet on SCC4######################
#####Functions that work and might be useful####################################
convertAdata <- function(PyAnnData) {
    # This function creates/updates an RAnnData Object from the given AnnData
    # from Python. 
    #TODO Not finished yet! Solve the PyAnnData$obsm parsing!
    RAnnData <- new("RAnnData", X = PyAnnData$X, obs = PyAnnData$obs, 
                    var = PyAnnData$var, uns = PyAnnData$uns, 
                    obsm = PyAnnData$obsm$as_dict())
    return(RAnnData)
}

setClass("RAnnData", slots = c(X = "matrix", obs = "data.frame", 
                               var = 'data.frame', uns = 'list', obsm = 'list'))
plotUMAP2D <- function(Adata, method = 'mlle') {
    if ((class(Adata) == c("anndata.core.anndata.AnnData", 
                           "python.builtin.object"))[1] && 
        (class(Adata) == c("anndata.core.anndata.AnnData", 
                           "python.builtin.object"))[2]) {
        # Condition that the input adata is the Python AnnData Object. 
        if (is.null(adata$obsm$get('X_vis_umap'))) {
            write('Calculating with method: MLLE', stdout())
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
################################################################################
st <- import('stream')
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
plotUMAP2D(Radata)
Radata <- convertAdata(adata)

#####DEBUGGING AREA#############################################################
