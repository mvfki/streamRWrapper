import sys
#import stream as st
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans
from collections import OrderedDict, defaultdict
import copy
#import anndata

def rgb216(color):
    '''Mainly for converting seaborn RGB color code (from 0 to 1) to 
    hexadecimal, which stream plotting methods support.'''
    c = '#%02x%02x%02x' % (int(round(color[0] * 256)), 
                           int(round(color[1] * 256)), 
                           int(round(color[2] * 256)))
    return c

class streamSingleCellSamples(object):
    """A streamSamples object stores information read and processed by stream 
    package and wraps standard pipeline commands as method functions"""
    def __init__(self, fileName, rawCount = True, cellLabel = None, 
                 cellLabelColor = None):
        '''Initialization
        fileName       - Required. Only support cellranger MTX format for now. 
                         The "barcodes.tsv" and "genes.tsv" file should be 
                         located at the exactly same path and named exactly as 
                         mentioned here.
        rawCount       - Boolean, default True. If True, stream.normalize_per_
                         cell() and stream.log_transform() will be performed. 
                         Note that this normalization produces RPM values 
                         rather than z-score. 
        cellLabel      - A TSV file with a label each line, and each line 
                         corresponds to a cell in "barcodes.tsv".
        cellLabelColor - A TSV file with each line written in "<label>\\t<hex-
                         color-code>" format.'''
        self.originalAdata = st.read(file_name = fileName, file_format = 'mtx')
        st.add_cell_labels(self.originalAdata, file_name = cellLabel)
        st.add_cell_colors(self.originalAdata, file_name = cellLabelColor)
        self.originalAdata.var_names_make_unique()
        self.originalAdata.obs_names_make_unique()
        print('Raw input parsed...')
        print(self.originalAdata)
        self.nCells = self.originalAdata.n_obs
        self.nGenes = self.originalAdata.n_vars
        self._keepCurrentRecords()
        st.remove_mt_genes(self.originalAdata)
        if rawCount:
            st.normalize_per_cell(self.originalAdata)
            st.log_transform(self.originalAdata)
        self.backup = {}
        self.backupKey = 0
        self.backup[0] = copy.deepcopy(self.originalAdata)
        print('Initial backup saved with key: 0')
        print('Restore with self.restoreFromBackup()')

    def backup(self, key = None):
        sys.stderr.write('WARNING: Maintaining backup takes memory. \n')
        if key == None:
            self.backupKey += 1
        tmp = copy.deepcopy(self.originalAdata)
        self.backup[self.backupKey] = tmp
        if type(key) == str:
            key = '"' + key + '"'
        print('Current results backuped with key:', key)
        print('Restore with self.restoreFromBackup(%s)' % str(key))

    def restoreFromBackup(self, key = 0):
        self.originalAdata = copy.deepcopy(self.backup)

    def preprocess(self, min_num_genes = None, min_num_cells = None, 
                   expr_cutoff = 1, loess_frac = 0.01, nb_pct = 0.1):
        '''Do the filtration to the dataset. 
        Arguments:
        min_num_genes - int, optional (default: None). 
                        Minimum number of genes expressed. For filtering 
                        cells. 
        min_num_cells - int, optional (default: max(5, (nCells * 0.001)) ). 
                        Minimum number of cells expressing one gene. For 
                        filtering genes. 
        expr_cutoff   - float, optional (default: 1). 
                        Expression cutoff. If greater than expr_cutoff, the 
                        gene is considered 'expressed'. 
        loess_frac    - float, optional (default: 0.01). Between 0 and 1.
                        The fraction of the data used when estimating each 
                        y-value in LOWESS function. For selecting variable 
                        genes. 
        nb_pct        - float, optional (default: 0.1). 
                        The percentage neighbor cells. For doing MLLE 
                        dimension reduction. Smaller number such as 0.01 is 
                        recommended when dealing with large dataset with more 
                        than 10,000 cells. 
        '''
        st.filter_cells(self.originalAdata, min_num_genes = min_num_genes)
        if min_num_cells == None:
            min_num_cells = int(round(self.originalAdata.shape[0] * 0.001))
        st.filter_genes(self.originalAdata, min_num_cells = min_num_cells, expr_cutoff = expr_cutoff) # updated from default # Then back to default
        st.select_variable_genes(self.originalAdata, loess_frac = loess_frac) # n_genes parameter updated
        st.dimension_reduction(self.originalAdata, nb_pct = nb_pct)
        print('Initial UMAP visualization')
        st.plot_visualization_2D(self.originalAdata, fig_legend_ncol = 4)
        #TODO KMeans

    def removeSmallCluster(self):
        pass
    
    def confirmRemoval(self):
        pass
    
    def _keepCurrentRecords(self):
        self.allCells = self.originalAdata.obs.index.to_list()
        self.label2cells = defaultdict(set)
        self.label2color = {}
        for cell, info in self.originalAdata.obs.iterrows():
            self.label2cells[info['label']].add(cell)
            if info['label'] in self.label2color:
                if self.label2color[info['label']] != info['label_color']:
                    sys.stderr.write('In cell %s - label %s - color %s, does not match with %s. Probably it is not a unique cell. \n' % (cell, info['label'], info['label_color'], self.label2color[info['label']]))
                    sys.stderr.write('Ingnoring cell %s. \n' % (cell))
            else:
                self.label2color[info['label']] = info['label_color']

    def _recoverRecords(self):
        self.originalAdata.uns['label_color'] = self.label2color
        existingCells = set(self.originalAdata.obs.index.to_list())
        for label, cellSet in self.label2cells.items():
            cellSet = cellSet.intersection(existingCells)
            for cell in cellSet:
                self.originalAdata.obs.at[cell, 'label'] = label
                self.originalAdata.obs.at[cell, 'label_color'] = self.label2color[label]

def labelByCluster(adata, kmeans):
    '''Change the sample labels in place. This function supports labeling the 
    cell's from an sklearn.cluster.KMeans object. Probably will also work with
    other sklearn.cluster methods but not tested yet. 
    Input
        adata  - anndata.AnnData object
        kmeans - sklearn.cluster.KMeans object, where you should have already 
                 run "KMeans.fit(adata.obsm['X_vis_umap'])"

    As a result, the label in the input anndata.AnnData object will be changed 
    automatically to 'c0', 'c1', 'c2', ... pattern. 
    '''
    nC = kmeans.n_clusters
    palettes = sns.color_palette("husl", nC)
    palettes = [rgb216(c) for c in palettes]
    adata.obs['label'] = 'unknown'
    adata.obs['label_color'] = 'gray'
    adata.uns['label_color'] = {}
    adata.uns['label_color']['unknown'] = 'gray'
    for i in range(nC):
        c0 = np.where(kmeans.labels_ == i)[0]
        label = 'c' + str(i)
        adata.obs['label'][c0] = label
        adata.obs['label_color'][c0] = palettes[i]
        adata.uns['label_color'][label] = palettes[i]

def removeCellsByLabel(adata, labels):
    '''In place remove all the cell with the given labels from the input 
    anndata.AnnData object. Usually used when having done a clustering and you 
    want to remove some clusters. 

    Input
        adata    - anndata.AnnData object
        labels - a str or a list. If a str, it should be exactly the label 
                   for that cluster, usually 'c0', 'c1', ... If a list, it 
                   should be a list of labels that exist.
    '''
    cellToRemove = []
    if type(labels) == str:
        cellToRemove = list(adata.obs['label'][adata.obs['label'] == labels].index)
    elif type(labels) == list:
        for c in labels:
            cellToRemove.extend(list(adata.obs['label'][adata.obs['label'] == c].index))
    cellToRemove = set(cellToRemove)
    remain = set(adata.obs.index).difference(cellToRemove)
    Idx = []
    for i in remain:
        Idx.append(list(adata.obs.index).index(i))
    adata._inplace_subset_obs(Idx)

class testClass(object):
    """docstring for testClass"""
    def __init__(self, arg):
        super(testClass, self).__init__()
        self.arg = arg

    def func1(self):
        return self.arg

    def func2(self, arg2):
        self.arg += arg2
        return self.arg
        

def testFunc(integer):
    integer = int(integer)
    return (type(integer), integer)