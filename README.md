# streamRWrapper
Wraps Python3 functions into R commands so that people can run STREAM based single-cell analysis in R and do anything else they like. 

## Prerequisites  
### STREAM installation  
reference link: https://github.com/pinellolab/STREAM
```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda create -n STREAM python=3.6 stream jupyter
$ conda activate STREAM
```

### AnnData installation
reference link: https://github.com/theislab/anndata
While creating the `STREAM` environment, it is already installed as a dependency. If it is necessary to work outside of conda environment, install with:  
```
pip install anndata
```
Make sure the `pip` is **not** a conda pip. 

## Features  
In `stream_R_Wrapper.R`, there is a pipeline that follows the [STREAM suggestions](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/1.STREAM_scRNA-seq.ipynb?flush_cache=true), and functions such as one for converting an Python [AnnData](https://github.com/theislab/anndata) object to R's [SingleCellExperiment](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) object which is widely used. 

## Note
There is a conflict between environment when running the script using R from `module load R` and importing modules from a conda based Python. 
Here are two alternative solution:
- Use R from the conda environment.
In this way, people who need a bunch of packages from another R version will complain; while people who have been working with a conda R will have a better experience. 
- Save the AnnData finished in Python to `.h5ad` format, and read to your R. 
You will need to install `anndata` as mentioned above. In R, do: 
```
library(reticulate)
anndata <- import('anndata')
adata <- anndata$read_h5ad(PATH_TO_THE_H5AD_ADATA)
sce <- adata2sce(adata)
```