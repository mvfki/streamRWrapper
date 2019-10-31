# THIS SCRIPT IS STILL A TEST NOTE


# Here I want to try to implement some instances or functions, so that we can run
# commands in R, which actually invoke functions from Python, but meanwhile the
# contents (i.e. eadable numbers, matrices) can also be called from R.

library(reticulate)
main <- import_main()
sys <- import("sys")



# Load Everything needed from Python scripts
use_virtualenv("STREAM") #TODO: This one seems still not working. For now run `conda activate STREAM` first before doing anything

source_python('streamUtil.py')

adata2 <- st.read('testData/matrix.mtx')

adata <- streamSingleCellSamples('testData/matrix.mtx', 'testData/cell_label.tsv', 'testData/cell_label_color.tsv')
# Get a dataframe structure with R commands
exprDF <- data.frame(adata$adata$X, row.names = adata$allCells)
colnames(exprDF) <- adata$allGenes


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
