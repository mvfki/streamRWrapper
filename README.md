# streamRWrapper
Wraps Python3 functions into R commands so that people can run STREAM based single-cell analysis in R and do anything else they like. 

## Prerequisites  
### STREAM installation  
reference link: https://github.com/pinellolab/STREAM  
Since there might be a conflict between Anaconda3 based environment VS computation cluster pre-installed modules, here two approaches to the installation are provided. If you work with R from Anaconda3, install STREAM with `conda`, elif you work with pre-installed R module, you might need to go for the manual pip installation.
#### Installation with Anaconda3
```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda create -n STREAM python=3.6 stream jupyter
$ conda activate STREAM
```
#### Conda free installation
Before this, **make sure** you check `which python` and `which pip` which are supposed to show you the path of the pre-installed modules. 
```
git clone https://github.com/pinellolab/STREAM.git
cd STREAM
python setup.py install --user
pip install python-slugify --user
pip install rpy2 --user
pip install anndata --user
# Above are the most possible missing packages that you need to manually install. 
# There might be other dependencies you need to manually insntall with pip.
```
After these, try import STREAM from Python
```
$ python
>>> import stream as st
```
If there is R shared library not found issue, try `echo $LD_LIBRARY_PATH` to see if the R shared library path is appended in your system. If not, append it.
```
LD_LIBRARY_PATH=PATH/TO/R/lib:$LD_LIBRARY_PATH
```

## Features  
In `stream_R_Wrapper.R`, there is a pipeline that follows the [STREAM suggestions](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/1.STREAM_scRNA-seq.ipynb?flush_cache=true), and functions such as one for converting an Python [AnnData](https://github.com/theislab/anndata) object to R's [SingleCellExperiment](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) object which is widely used. 
