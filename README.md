# streamRWrapper
Wraps Python3 functions into R commands so that people can run STREAM based single-cell analysis in R and do anything else they like. 

## Prerequisites  
### STREAM installation  
reference link: https://github.com/pinellolab/STREAM  
Since there might be a conflict between Anaconda3 based environment VS computation cluster pre-installed modules, here two approaches to the installation are provided. If you work with R from Anaconda3, install STREAM with `conda`, elif you work with pre-installed R module, you might need to go for the manual pip installation.
#### Installation with Anaconda3
```console
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda create -n STREAM python=3.6 stream jupyter
$ conda activate STREAM
```
#### Conda free installation
Before this, **make sure** you check `which python` and `which pip` which are supposed to show you the path of the pre-installed modules. 
```{bash}
git clone https://github.com/pinellolab/STREAM.git
cd STREAM
python setup.py install --user
pip install -r pip_requirements.txt --user
```
One thing that was terribly annoying is the changing argument naming in the dependency `networkx`. By default, STREAM requires networkx 2.1 (You can see this when installing from conda). However, in my preinstalled environment, networkx is under version 2.3, where there are different argument naming for equivalent values. 
The version of your networkx can be checked by running in bash command
```{bash}
pip list | grep networkx
```
If so, I've got the solution. Run bash commands in the stream git repo directory you just cloned before you install STREAM. 
```{bash}
cp PATH/TO/THIS/REPO/STREAM_diff.patch ./
git apply STREAM_diff.patch
```
The most possible missing packages are listed in `pip_requirements.txt`, but there can also be other missing ones for different users. You might need to install them by yourself.  
After these, try import STREAM from Python to check if the installation is done.
```{python}
import stream as st
```
If there is R shared library not found issue, try `echo $LD_LIBRARY_PATH` to see if the R shared library path is appended in your system. If not, append it.
```{bash}
LD_LIBRARY_PATH=PATH/TO/R/lib:$LD_LIBRARY_PATH
```
Besides, there are also some additional R packages needed for STREAM to work. 
```{r}
if(!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("Albluca/distutils")
devtools::install_github("Albluca/ElPiGraph.R")
```
## Installation
In R, install the package from this repository with `devtools`:
```{r}
library(devtools)
install_github('mvfki/streamRWrapper')
```
Note that this would only work when I set the repo to public. 
## Features  
In `stream_R_Wrapper.R`, there is the minimum steps to load the necessary environment, a pipeline for `testData_real/` that follows the [STREAM suggestions](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/1.STREAM_scRNA-seq.ipynb?flush_cache=true), and functions such as one for converting an Python [AnnData](https://github.com/theislab/anndata) object to R's [SingleCellExperiment](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) object which is widely used. All the functions can be installed with packaging. 
