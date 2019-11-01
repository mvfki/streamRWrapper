# THIS SCRIPT IS STILL A TEST NOTE

# Basically with package "reticulate", I can call most of the Python functions 
# from R. I guess now what I want to do is to just run the same thing here and 
# reproduce the results that I can get from directly running in Python. 
# Next thing to do is to make things more R, including converting Python Objects
# saved visualable in R, callable with R functions. 

# This first section would save the commands that worked everywhere so far.
library(reticulate)
use_condaenv(condaenv = 'STREAM', required = TRUE)
Sys.setenv(LD_LIBRARY_PATH="~/.conda/envs/STREAM/lib/")

#####These commands worked on my local but not yet on SCC4######################
#####Functions that work and might be useful####################################
defaultMinNumCells <- function(adata) {
    shapeFactor <- adata$n_obs
    toCompare <- c(5, shapeFactor)
    ans <- max(toCompare)
    return(ans)
}
################################################################################
st <- import('stream')
adata <- st$read('testData/matrix.mtx', file_format = 'mtx')
st$add_cell_labels(adata, file_name = 'testData/cell_label.tsv')
st$add_cell_colors(adata, file_name = 'testData/cell_label_color.tsv')
adata$obs_names_make_unique()
adata$var_names_make_unique()
st$remove_mt_genes(adata)
st$normalize_per_cell(adata)
st$log_transform(adata)
st$filter_cells(adata)
# In small test case I say "1" there, but in real data remember to say "5". 
st$filter_genes(adata, min_num_cells = max(1L, adata$n_obs * 0.001))
# The "L" after the integer in the following line is important as it specifies
# the datatype as integer rather than a float which might not be accepbtable in 
# Python. 
st$select_variable_genes(adata, n_genes = min(2000L, adata$n_vars))
st$dimension_reduction(adata, nb_pct = 0.01)
#####DEBUGGING AREA#############################################################

library(methods)
setClass("AnnData", representation(a = "character", b = "numeric"))
readMTX <- function(fileName = "character") {
    # This function would run the Python stream.read(), and then from the
    # Python side access some basic informations and save them in the AnnData
    # class instance defined here.
    adata <- new("AnnData")
    return(adata)
}
x <- readMTX()
