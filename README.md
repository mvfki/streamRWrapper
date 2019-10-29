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
