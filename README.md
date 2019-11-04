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

## Features  
In `stream_R_Wrapper.R`, there is a pipeline that follows the [STREAM suggestions](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/1.STREAM_scRNA-seq.ipynb?flush_cache=true), and functions such as one for converting an Python [AnnData](https://github.com/theislab/anndata) object to R's [SingleCellExperiment](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) object which is widely used. 
