# Here I want to try to implement some instances or functions, so that we can run
# commands in R, which actually invoke functions from Python, but meanwhile the
# contents (i.e. eadable numbers, matrices) can also be called from R.

library(reticulate)
# Load Everything needed from Python scripts
use_virtualenv("STREAM")
source_python('streamUtil.py')
# Note that after trying a bit, a Python function is actually easy to invoke.
# But nothing visible is saved in R workspace except the physical memory
# address of the functions and variables.

# Basic idea of how this wrapper works:
# In R                                  In Python
# 1. read raw data as AnnData  ___
#                                 |
#                                 ----> 2. Convert data to Python and
#                                 _____    calculate with STREAM
# 3. Call back the calculated <---|
#    AnnData, converted.

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
